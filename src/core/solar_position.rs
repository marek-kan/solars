use crate::core::utils::limit_deg_to_360;

/// result is dtau (degrees)
/// `r = calculate_heliocentric_coeff(..., CoordType::Radius)`
pub(crate) fn calculate_aberration_correction(r: f64) -> f64 {
    -20.4898 / (3600.0 * r)
}

/// result is lambda (degrees)
/// `theta (Geocentric longitude) = calculate_geocentric_coeff(c: l, CoordType::Longitude)` (L = Heliocentric longitude)
pub(crate) fn calculate_apparent_sun_longitude(theta: f64, dpsi: f64, dtau: f64) -> f64 {
    theta + dpsi + dtau
}

/// v (degrees)
pub(crate) fn sidereal_time_greenwich(jd: f64, jce: f64, dpsi: f64, epsilon: f64) -> f64 {
    let v_0 = limit_deg_to_360(
        280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * jce.powi(2)
            - jce.powi(3) / 38710000.0,
    );

    limit_deg_to_360(v_0 + dpsi * epsilon.to_radians().cos())
}

/// alpha (degrees)
pub(crate) fn sun_right_ascension(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    let l = lambda.to_radians();
    let e = epsilon.to_radians();
    let b = beta.to_radians();

    (l.sin() * e.cos() - b.tan() * e.sin())
        .atan2(l.cos())
        .to_degrees()
}

/// delta (degrees)
pub(crate) fn geocentric_sun_declination(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    let b = beta.to_radians();
    let e = epsilon.to_radians();
    let l = lambda.to_radians();

    (b.sin() * e.cos() + b.cos() * e.sin() * l.sin())
        .asin()
        // .to_degrees()
}

/// H (degrees)
pub(crate) fn obs_local_hour_angle(lon: f64, v: f64, alpha: f64) -> f64 {
    limit_deg_to_360(v + lon - alpha)
}

/// delta' (in degrees)
pub(crate) fn calculate_topocentric_sun_right_ascension_declination(lat: f64, elevation: f64, r: f64, hour_angle: f64, geoc_sun_declination: f64, sun_right_ascension: f64) -> f64 {
    let lat_rad = lat.to_radians();

    let e = 8.794 / (3600.0 * r);
    let u = 0.99664719 * lat_rad.tan();
    let x = u.cos() + elevation / 6378140.0 * lat_rad.cos();
    let y = 0.99664719 * u.sin() + elevation / 6378140.0 * lat_rad.sin();

    let d_alpha = (-x * e.sin() * hour_angle.sin()).atan2(geoc_sun_declination.cos() * e.sin() * hour_angle.cos());
    
    // alpha_ = sun_right_ascension + d_alpha.to_degrees();

    ((geoc_sun_declination.sin() - y * e.sin()) * d_alpha.cos()).atan2(geoc_sun_declination.cos() - x * e.sin() * hour_angle.cos()) // 39
}
