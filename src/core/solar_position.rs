use crate::core::utils::limit_deg_to_360;

/// dtau (radians)
/// `r = calculate_heliocentric_coeff(..., CoordType::Radius)`
pub(crate) fn aberration_correction(r: f64) -> f64 {
    (-20.4898 / (3600.0 * r)).to_radians()
}

/// lambda (radians)
/// `theta (Geocentric longitude) = calculate_geocentric_coeff(c: l, CoordType::Longitude)` (L = Heliocentric longitude)
pub(crate) fn apparent_sun_longitude(theta: f64, dpsi: f64, dtau: f64) -> f64 {
    theta + dpsi + dtau
}

/// v (radians)
pub(crate) fn sidereal_time_greenwich(jd: f64, jce: f64, dpsi: f64, epsilon: f64) -> f64 {
    let v_0 = (280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * jce.powi(2)
        - jce.powi(3) / 38710000.0)
        .to_radians()
        .rem_euclid(std::f64::consts::TAU);

    (v_0 + dpsi * epsilon.cos()).rem_euclid(std::f64::consts::TAU)
}

/// alpha (radians)
pub(crate) fn geocentric_sun_right_ascension(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    ((lambda.sin() * epsilon.cos() - beta.tan() * epsilon.sin()).atan2(lambda.cos()))
        .rem_euclid(std::f64::consts::TAU)
}

/// delta (radians)
pub(crate) fn geocentric_sun_declination(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    (beta.sin() * epsilon.cos() + beta.cos() * epsilon.sin() * lambda.sin()).asin()
}

/// H (radians)
pub(crate) fn obs_local_hour_angle(lon: f64, v: f64, alpha: f64) -> f64 {
    (v + lon - alpha).rem_euclid(std::f64::consts::TAU)
}

/// delta' and topocentric hour angle (radians)
pub(crate) fn topocentric_sun_right_ascension_declination_and_local_hour_angle(
    lat: f64,
    elevation: f64,
    r: f64,
    hour_angle: f64,
    geoc_sun_declination: f64,
    // gec_sun_right_ascension: f64,
) -> (f64, f64) {
    let e_rad = (8.794 / (3600.0 * r)).to_radians();

    let u = (0.99664719 * lat.tan()).atan();
    let x = u.cos() + (elevation / 6378140.0) * lat.cos();
    let y = 0.99664719 * u.sin() + (elevation / 6378140.0) * lat.sin();

    let d_alpha_num = -x * e_rad.sin() * hour_angle.sin();
    let d_alpha_den = geoc_sun_declination.cos() - x * e_rad.sin() * hour_angle.cos();

    let d_alpha_rad = d_alpha_num.atan2(d_alpha_den);

    // delta'
    let delta_ =
        ((geoc_sun_declination.sin() - y * e_rad.sin()) * d_alpha_rad.cos()).atan2(d_alpha_den);
    // H'
    let h_ = hour_angle - d_alpha_rad;

    (delta_, h_)
}

/// temp in deg. Celsius, pres in milibars; returns radians
pub(crate) fn topocentric_elevation_angle(
    lat: f64,
    topocentric_sun_declination: f64,
    topocentric_hour_angle: f64,
    pressure: Option<f64>,
    temperature: Option<f64>,
) -> f64 {
    let e0_uncorr = topocentric_elevation_angle_wo_correction(
        lat,
        topocentric_sun_declination,
        topocentric_hour_angle,
    );

    if let (Some(p), Some(t)) = (pressure, temperature) {
        let delta_e = atmospheric_refraction_angle(e0_uncorr, p, t);

        return e0_uncorr + delta_e;
    };

    println!(
        "Calculationg uncorrected topocentric elevation angle. For correction provide pressure and temperature."
    );

    e0_uncorr
}

/// e0 (radians)
pub(crate) fn topocentric_elevation_angle_wo_correction(
    lat: f64,
    topocentric_sun_declination: f64,
    topocentric_hour_angle: f64,
) -> f64 {
    (lat.sin() * topocentric_sun_declination.sin()
        + lat.cos() * topocentric_sun_declination.cos() * topocentric_hour_angle.cos())
    .asin()
}

/// delta_e0 (radians)
/// temp in deg. Celsius, pres in milibars
pub(crate) fn atmospheric_refraction_angle(
    topocentric_elevation_angle: f64,
    pressure: f64,
    temperature: f64,
) -> f64 {
    let p = pressure / 1010.0;
    let t = 283.0 / (273.0 + temperature);
    let elevation_degrees = topocentric_elevation_angle.to_degrees();

    let e_degrees = 1.02
        / (60.0
            * (elevation_degrees + (10.3 / (elevation_degrees + 5.11)))
                .to_radians()
                .tan());

    let delta_e = p * t * e_degrees;

    delta_e.to_radians()
}

/// phi (degrees)
/// topocentric_elevation angle in degrees
pub(crate) fn topocentric_zenith_angle(topocentric_elevation_angle: f64) -> f64 {
    90.0 - topocentric_elevation_angle
}

/// capital gamma (radians westward from south)
pub(crate) fn topocentric_azimuth_angle_w_from_s(
    lat: f64,
    topocentric_hour_angle: f64,
    topocentric_sun_declination: f64,
) -> f64 {
    topocentric_hour_angle.sin().atan2(
        topocentric_hour_angle.cos() * lat.sin() - topocentric_sun_declination.tan() * lat.cos(),
    )
}

/// capital phi (degrees eastward from north)
/// for navigators and solar radiotion
/// topocentric_azimuth_angle_w_from_south: degrees
pub(crate) fn topocentric_azimuth_angle_e_from_n(
    topocentric_azimuth_angle_w_from_south: f64,
) -> f64 {
    limit_deg_to_360(topocentric_azimuth_angle_w_from_south + 180.0)
}

/// all inputs should be degrees
pub fn aoi(
    zenith_angle: f64,
    topocentric_azimuth_from_north: f64,
    panel_tilt: f64,
    panel_azimuth: f64,
) -> f64 {
    let zenith_rad = zenith_angle.to_radians();
    let tilt_rad = panel_tilt.to_radians();
    let azimuth_diff_rad = (topocentric_azimuth_from_north - panel_azimuth).to_radians();

    let cos_aoi = zenith_rad.cos() * tilt_rad.cos()
        + zenith_rad.sin() * tilt_rad.sin() * azimuth_diff_rad.cos();

    cos_aoi.clamp(-1.0, 1.0).acos().to_degrees()
}
