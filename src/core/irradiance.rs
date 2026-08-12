/// Parameters for the uncoated-glass physical incidence angle modifier.
///
/// `refractive_index` is unitless, `extinction_coefficient` is in 1/m, and
/// `thickness` is in meters.
#[derive(Clone, Copy, Debug)]
pub(crate) struct OpticalLossParameters {
    pub(crate) refractive_index: f64,
    pub(crate) extinction_coefficient: f64,
    pub(crate) thickness: f64,
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
pub(crate) fn aoi(
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
        physical_optical_loss(geometric_aoi, parameters)
    })
}

pub(crate) fn etraterrestrial_radiation(day_of_year: i64) -> f64 {
    let x = (2.0 * std::f64::consts::PI * (day_of_year - 1) as f64) / 365.0;

    1366.1
        * (1.00011 + 0.034221 * x.cos() + 0.00128 * x.sin() - 0.000719 * (2.0 * x).cos()
            + 0.000077 * (2.0 * x).sin())
}

/// Clear-sky irradiance components calculated with the Ineichen/Perez model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ClearSkyIrradiance {
    pub(crate) ghi: f64,
    pub(crate) dni: f64,
    pub(crate) dhi: f64,
}

/// Calculate relative airmass using the Kasten-Young 1989 model.
///
/// `zenith` is in degrees. Values at or below the horizon return infinity.
pub(crate) fn relative_airmass_kasten(zenith: f64) -> f64 {
    if zenith >= 90.0 {
        return f64::INFINITY;
    }

    let zenith_radians = zenith.to_radians();
    1.0 / (zenith_radians.cos() + 0.50572 * (96.07995 - zenith).powf(-1.6364))
}

/// Calculate scalar clear-sky GHI, DNI, and DHI with the Ineichen/Perez model.
///
/// Angles are degrees, elevation is meters, and irradiance values are W/m2.
/// `pressure`, when supplied, is in millibars and corrects relative airmass
/// to absolute airmass. The Perez enhancement factor is not applied.
pub(crate) fn ineichen_clearsky(
    zenith: f64,
    linke_turbidity: f64,
    elevation: f64,
    pressure: Option<f64>,
    dni_extra: f64,
) -> ClearSkyIrradiance {
    let cos_zenith = zenith.to_radians().cos().max(0.0);
    if cos_zenith == 0.0 {
        return ClearSkyIrradiance {
            ghi: 0.0,
            dni: 0.0,
            dhi: 0.0,
        };
    }

    let relative_airmass = relative_airmass_kasten(zenith);
    let airmass = pressure.map_or(relative_airmass, |value| relative_airmass * value / 1013.25);
    let fh1 = (-elevation / 8000.0).exp();
    let fh2 = (-elevation / 1250.0).exp();
    let cg1 = 5.09e-5 * elevation + 0.868;
    let cg2 = 3.92e-5 * elevation + 0.0387;

    let ghi = cg1
        * dni_extra
        * cos_zenith
        * (-cg2 * airmass * (fh1 + fh2 * (linke_turbidity - 1.0)))
            .exp()
            .max(0.0);
    let b = 0.664 + 0.163 / fh1;
    let bnci = dni_extra * (b * (-0.09 * airmass * (linke_turbidity - 1.0)).exp()).max(0.0);
    let bnci_2 = ghi
        * ((1.0 - (0.1 - 0.2 * (-linke_turbidity).exp()) / (0.1 + 0.882 / fh1)) / cos_zenith)
            .clamp(0.0, 1e20);
    let dni = bnci.min(bnci_2);
    let dhi = ghi - dni * cos_zenith;

    ClearSkyIrradiance { ghi, dni, dhi }
}

/// Plane-of-array irradiance components calculated with the Hay-Davies model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct PoaIrradiance {
    pub(crate) global: f64,
    pub(crate) direct: f64,
    pub(crate) diffuse: f64,
    pub(crate) sky_diffuse: f64,
    pub(crate) ground_diffuse: f64,
}

/// Calculate scalar plane-of-array irradiance using the Hay-Davies model.
///
/// Angles are degrees and irradiance values are W/m2. `dni_extra` is the
/// extraterrestrial direct normal irradiance and `albedo` is the ground
/// reflectance factor.
pub(crate) fn poa_haydavies(
    panel_tilt: f64,
    panel_azimuth: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
    dni: f64,
    ghi: f64,
    dhi: f64,
    dni_extra: f64,
    albedo: f64,
) -> PoaIrradiance {
    let zenith_rad = solar_zenith.to_radians();
    let tilt_rad = panel_tilt.to_radians();

    let aoi_projection = zenith_rad.cos() * tilt_rad.cos()
        + zenith_rad.sin() * tilt_rad.sin() * (solar_azimuth - panel_azimuth).to_radians().cos();
    let aoi_projection = aoi_projection.max(0.0);
    let zenith_projection = zenith_rad.cos().max(0.01745);
    let projection_ratio = aoi_projection / zenith_projection;
    let anisotropy_index = dni / dni_extra;
    let sky_view_factor = (1.0 + tilt_rad.cos()) / 2.0;

    let isotropic = (dhi * (1.0 - anisotropy_index) * sky_view_factor).max(0.0);
    let circumsolar = (dhi * anisotropy_index * projection_ratio).max(0.0);
    let sky_diffuse = isotropic + circumsolar;
    let ground_diffuse = ghi * albedo * (1.0 - tilt_rad.cos()) / 2.0;
    let direct = dni * aoi_projection;
    let diffuse = sky_diffuse + ground_diffuse;

    PoaIrradiance {
        global: direct + diffuse,
        direct,
        diffuse,
        sky_diffuse,
        ground_diffuse,
    }
}

/// Calculate uncoated-glass physical incidence angle modifier.
pub(crate) fn physical_optical_loss(aoi: f64, parameters: OpticalLossParameters) -> f64 {
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
    use super::{
        OpticalLossParameters, aoi, ineichen_clearsky, physical_optical_loss, poa_haydavies,
        relative_airmass_kasten,
    };

    #[test]
    fn geometric_aoi() {
        let actual = aoi(50.0, 200.0, 30.0, 180.0, None);
        assert!((actual - 23.566_943_771_139).abs() < 1e-12);
    }

    #[test]
    fn physical_optical_loss_ar_coating() {
        let parameters = OpticalLossParameters::default();
        let actual = aoi(50.0, 200.0, 30.0, 180.0, Some(parameters));
        assert!((actual - 0.999_133_283_938).abs() < 1e-12);
        assert!((physical_optical_loss(60.0, parameters) - 0.946_002_914_223).abs() < 1e-12);
    }

    #[test]
    fn haydavies_poa_daytime_conditions() {
        let actual = poa_haydavies(30.0, 180.0, 40.0, 190.0, 800.0, 600.0, 120.0, 1367.0, 0.2);

        assert!((actual.global - 928.251_756_634_908).abs() < 1e-10);
        assert!((actual.direct - 783.940_047_158_946).abs() < 1e-10);
        assert!((actual.diffuse - 144.311_709_475_962).abs() < 1e-10);
        assert!((actual.sky_diffuse - 136.273_233_703_029).abs() < 1e-10);
        assert!((actual.ground_diffuse - 8.038_475_772_934).abs() < 1e-10);
    }

    #[test]
    fn haydavies_poa_near_the_horizon() {
        let actual = poa_haydavies(45.0, 135.0, 84.0, 110.0, 200.0, 90.0, 60.0, 1414.0, 0.25);

        assert!((actual.global - 247.262_591_285_123).abs() < 1e-10);
        assert!((actual.direct - 142.251_697_788_909).abs() < 1e-10);
        assert!((actual.diffuse - 105.010_893_496_214).abs() < 1e-10);
        assert!((actual.sky_diffuse - 101.715_844_784_563).abs() < 1e-10);
        assert!((actual.ground_diffuse - 3.295_048_711_651).abs() < 1e-10);
    }

    #[test]
    fn ineichen_clearsky_at_sea_level() {
        let actual = ineichen_clearsky(40.0, 3.0, 0.0, None, 1367.0);

        assert!((relative_airmass_kasten(40.0) - 1.304_223_540_914).abs() < 1e-12);
        assert!((actual.ghi - 781.234_083_806_219).abs() < 1e-10);
        assert!((actual.dni - 893.961_773_580_715).abs() < 1e-10);
        assert!((actual.dhi - 96.419_634_793_926).abs() < 1e-10);
    }

    #[test]
    fn ineichen_clearsky_with_millibar_pressure() {
        let actual = ineichen_clearsky(65.0, 4.2, 1830.0, Some(820.0), 1414.0);

        assert!((relative_airmass_kasten(65.0) - 2.356_019_282_209).abs() < 1e-12);
        assert!((actual.ghi - 415.671_524_712_685).abs() < 1e-10);
        assert!((actual.dni - 709.473_968_460_860).abs() < 1e-10);
        assert!((actual.dhi - 115.834_869_411_481).abs() < 1e-10);
    }
}
