use crate::core::utils::limit_deg_to_360;

/// result is dtau
/// `r = calculate_heliocentric_coeff(..., CoordType::Radius)`
pub(crate) fn calculate_aberration_correction(r: f64) -> f64 {
    -20.4898 / (3600.0 * r)
}

/// result is lambda
/// `theta (Geocentric longitude) = calculate_geocentric_coeff(c: l, CoordType::Longitude)` (L = Heliocentric longitude)
pub(crate) fn calculate_apparent_sun_longitude(theta: f64, dpsi: f64, dtau: f64) -> f64 {
    theta + dpsi + dtau
}

pub(crate) fn sidereal_time_greenwich(jd: f64, jce: f64, dpsi: f64, epsilon: f64) -> f64 {
    let v_0 = 80.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * jce.powi(2)
        - jce.powi(3) / 38710000.0;

    // limit v_0 to [0, 360] deg
    limit_deg_to_360(v_0) + dpsi * epsilon.cos()
}
