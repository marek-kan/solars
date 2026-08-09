use crate::core::utils::limit_deg_to_360;

/// dtau (degrees)
/// `r = calculate_heliocentric_coeff(..., CoordType::Radius)`
pub fn aberration_correction(r: f64) -> f64 {
    -20.4898 / (3600.0 * r)
}

/// lambda (degrees)
/// `theta (Geocentric longitude) = calculate_geocentric_coeff(c: l, CoordType::Longitude)` (L = Heliocentric longitude)
pub fn apparent_sun_longitude(theta: f64, dpsi: f64, dtau: f64) -> f64 {
    theta + dpsi + dtau
}

/// v (degrees)
pub fn sidereal_time_greenwich(jd: f64, jce: f64, dpsi: f64, epsilon: f64) -> f64 {
    let v_0 = limit_deg_to_360(
        280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * jce.powi(2)
            - jce.powi(3) / 38710000.0,
    );

    limit_deg_to_360(v_0 + dpsi * epsilon.to_radians().cos())
}

/// alpha (degrees)
pub fn geocentric_sun_right_ascension(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    let l = lambda.to_radians();
    let e = epsilon.to_radians();
    let b = beta.to_radians();

    limit_deg_to_360(((l.sin() * e.cos() - b.tan() * e.sin()).atan2(l.cos())).to_degrees())
}

/// delta (radians)
pub fn geocentric_sun_declination(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    let b = beta.to_radians();
    let e = epsilon.to_radians();
    let l = lambda.to_radians();

    (b.sin() * e.cos() + b.cos() * e.sin() * l.sin()).asin()
}

/// H (degrees)
pub fn obs_local_hour_angle(lon: f64, v: f64, alpha: f64) -> f64 {
    limit_deg_to_360(v + lon - alpha)
}

/// delta' (in rad)
/// topocentric_hour_angle (degrees)
pub fn topocentric_sun_right_ascension_declination_and_local_hour_angle(
    lat: f64,
    elevation: f64,
    r: f64,
    hour_angle: f64,
    geoc_sun_declination: f64,
    // gec_sun_right_ascension: f64,
) -> (f64, f64) {
    println!("Geoc sun decli: {geoc_sun_declination}");
    let lat_rad = lat.to_radians();
    let hour_rad = hour_angle.to_radians();

    let e = 8.794 / (3600.0 * r);
    let e_rad = e.to_radians();

    let u = (0.99664719 * lat_rad.tan()).atan();
    let x = u.cos() + (elevation / 6378140.0) * lat_rad.cos();
    let y = 0.99664719 * u.sin() + (elevation / 6378140.0) * lat_rad.sin();

    let d_alpha_num = -x * e_rad.sin() * hour_rad.sin();
    let d_alpha_den = geoc_sun_declination.cos() - x * e_rad.sin() * hour_rad.cos();

    let d_alpha_rad = d_alpha_num.atan2(d_alpha_den);

    // delta'
    let delta_ = ((geoc_sun_declination.sin() - y * e_rad.sin()) * d_alpha_rad.cos())
        .atan2(d_alpha_den)
        .to_degrees();
    // H'
    let h_ = hour_angle - d_alpha_rad.to_degrees();
    // alpha'
    // let alpha_ = gec_sun_right_ascension + d_alpha;
    // println!("alpha': {alpha_}");

    (delta_, h_)
}

/// temp in deg. Celsius, pres in milibars
pub fn topocentric_elevation_angle(
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

/// e0 (degrees)
pub fn topocentric_elevation_angle_wo_correction(
    lat: f64,
    topocentric_sun_declination: f64,
    topocentric_hour_angle: f64,
) -> f64 {
    let lat_rad = lat.to_radians();
    (lat_rad.sin() * topocentric_sun_declination.sin()
        + lat_rad.cos() * topocentric_sun_declination.cos() * topocentric_hour_angle.cos())
    .asin()
    .to_degrees()
}

/// delta_e0 (degrees)
/// temp in deg. Celsius, pres in milibars
pub fn atmospheric_refraction_angle(
    topocentric_elevation_angle: f64,
    pressure: f64,
    temperature: f64,
) -> f64 {
    let p = pressure / 1010.0;
    let t = 283.0 / (273.0 + temperature);

    let e = 1.02
        / (60.0
            * (topocentric_elevation_angle + (10.3 / (topocentric_elevation_angle + 5.11)))
                .to_radians()
                .tan());

    let delta_e = p * t * e;

    delta_e
}

/// phi (degrees)
/// topocentric_elevation angle in degrees
pub fn topocentric_zenith_angle(topocentric_elevation_angle: f64) -> f64 {
    90.0 - topocentric_elevation_angle
}

/// capital gamma (degrees westward from south)
pub fn topocentric_azimuth_angle_w_from_s(
    topocentric_hour_angle: f64,
    lat: f64,
    topocentric_sun_declination: f64,
) -> f64 {
    let lat_rad = lat.to_radians();

    topocentric_hour_angle
        .sin()
        .atan2(
            topocentric_hour_angle.cos() * lat_rad.sin()
                - topocentric_sun_declination.tan() * lat_rad.cos(),
        )
        .to_degrees()
}

/// capital phi (degrees eastward from north)
/// for navigators and solar radiotion
pub fn topocentric_azimuth_angle_e_from_n(topocentric_azimuth_angle_w_from_south: f64) -> f64 {
    limit_deg_to_360(topocentric_azimuth_angle_w_from_south + 180.0)
}

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
