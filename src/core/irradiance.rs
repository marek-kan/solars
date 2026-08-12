/// Parameters for the uncoated-glass physical incidence angle modifier.
///
/// `refractive_index` is unitless, `extinction_coefficient` is in 1/m, and
/// `thickness` is in meters.
#[derive(Clone, Copy, Debug)]
pub struct OpticalLossParameters {
    pub refractive_index: f64,
    pub extinction_coefficient: f64,
    pub thickness: f64,
}

impl Default for OpticalLossParameters {
    fn default() -> Self {
        Self {
            refractive_index: 1.526,
            extinction_coefficient: 4.0,
            thickness: 0.002,
        }
    }
}

/// Calculate geometric AOI in degrees, or physical IAM when optical parameters are provided.
///
/// All geometric inputs are degrees. Without `optical_loss`, this returns the
/// geometric angle of incidence in degrees. With `optical_loss`, it returns
/// the dimensionless physical incidence angle modifier for uncoated glass.
pub fn aoi(
    zenith_angle: f64,
    topocentric_azimuth_from_north: f64,
    panel_tilt: f64,
    panel_azimuth: f64,
    optical_loss: Option<OpticalLossParameters>,
) -> f64 {
    let zenith_rad = zenith_angle.to_radians();
    let tilt_rad = panel_tilt.to_radians();
    let azimuth_diff_rad = (topocentric_azimuth_from_north - panel_azimuth).to_radians();
    let cos_aoi = zenith_rad.cos() * tilt_rad.cos()
        + zenith_rad.sin() * tilt_rad.sin() * azimuth_diff_rad.cos();
    let geometric_aoi = cos_aoi.clamp(-1.0, 1.0).acos().to_degrees();

    optical_loss.map_or(geometric_aoi, |parameters| {
        physical(geometric_aoi, parameters)
    })
}

/// Calculate pvlib's uncoated-glass physical incidence angle modifier.
pub fn physical(aoi: f64, parameters: OpticalLossParameters) -> f64 {
    let n = parameters.refractive_index;
    let cos_incidence = aoi.to_radians().cos().max(0.0);
    if cos_incidence == 0.0 {
        return 0.0;
    }

    let sin_incidence = (1.0 - cos_incidence.powi(2)).sqrt();
    let sin_refraction = sin_incidence / n;
    let cos_refraction = (1.0 - sin_refraction.powi(2)).sqrt();

    let rho_s =
        ((cos_incidence - n * cos_refraction) / (cos_incidence + n * cos_refraction)).powi(2);
    let rho_p =
        ((cos_refraction - n * cos_incidence) / (cos_refraction + n * cos_incidence)).powi(2);
    let rho_normal = ((1.0 - n) / (1.0 + n)).powi(2);

    let absorption =
        (-parameters.extinction_coefficient * parameters.thickness / cos_refraction).exp();
    let normal_absorption = (-parameters.extinction_coefficient * parameters.thickness).exp();
    let transmitted = ((1.0 - rho_s) + (1.0 - rho_p)) * absorption / 2.0;
    let normal_transmitted = (1.0 - rho_normal) * normal_absorption;

    transmitted / normal_transmitted
}

#[cfg(test)]
mod tests {
    use super::{OpticalLossParameters, aoi, physical};

    #[test]
    fn geometric_aoi_matches_pvlib() {
        let actual = aoi(50.0, 200.0, 30.0, 180.0, None);
        assert!((actual - 23.566_943_771_139).abs() < 1e-12);
    }

    #[test]
    fn physical_iam_matches_pvlib_without_ar_coating() {
        let parameters = OpticalLossParameters::default();
        let actual = aoi(50.0, 200.0, 30.0, 180.0, Some(parameters));
        assert!((actual - 0.999_133_283_938).abs() < 1e-12);
        assert!((physical(60.0, parameters) - 0.946_002_914_223).abs() < 1e-12);
    }
}
