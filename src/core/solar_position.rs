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
