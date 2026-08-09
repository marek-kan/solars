use crate::core::solar_position::*;
use crate::core::time::{calc_julian_day, julian_century, julian_millennium};
use crate::core::utils::round_to_decimals;
use crate::{periodic_tables::earth::*, periodic_tables::nutation::*};
use chrono::{DateTime, TimeZone, Utc};

fn get_julian_date_values() -> (f64, f64, f64, f64) {
    let test_date: DateTime<Utc> = Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap();

    let jd_ut = calc_julian_day(&test_date);

    assert_eq!(
        2452930.313,
        round_to_decimals(jd_ut, 3),
        "Julian day failure"
    );

    let jd_tt = jd_ut + 67.0 / 86400.0; // In the example they mention dT

    let jc = julian_century(jd_tt);
    let jm = julian_millennium(jc);

    (jd_ut, jd_tt, jc, jm)
}

#[test]
fn end_to_end_test() {
    let (jd_ut, _, jc, jm) = get_julian_date_values();
    let (lat, lon): (f64, f64) = (39.742476, -105.1786);
    let (lat_rad, lon_rad) = (lat.to_radians(), lon.to_radians());

    let l0 = sum_table(&L0_TABLE, &jm);
    let l1 = sum_table(&L1_TABLE, &jm);
    let l2 = sum_table(&L2_TABLE, &jm);
    let l3 = sum_table(&L3_TABLE, &jm);
    let l4 = sum_table(&L4_TABLE, &jm);
    let l5 = sum_table(&L5_TABLE, &jm);

    let b0 = sum_table(&B0_TABLE, &jm);
    let b1 = sum_table(&B1_TABLE, &jm);

    let r0 = sum_table(&R0_TABLE, &jm);
    let r1 = sum_table(&R1_TABLE, &jm);
    let r2 = sum_table(&R2_TABLE, &jm);
    let r3 = sum_table(&R3_TABLE, &jm);
    let r4 = sum_table(&R4_TABLE, &jm);

    assert_eq!(
        172067561.527,
        round_to_decimals(l0, 3),
        "L0 failure, got {l0}"
    );
    assert_eq!(
        628332010650.051,
        round_to_decimals(l1, 3),
        "L1 failure, got {l1}"
    );
    assert_eq!(61368.682, round_to_decimals(l2, 3), "L2 failure, got {l2}");
    assert_eq!(-26.903, round_to_decimals(l3, 3), "L3 failure, got {l3}");
    assert_eq!(-121.280, round_to_decimals(l4, 3), "L4 failure, got {l4}");
    assert_eq!(-1.000, round_to_decimals(l5, 3), "L5 failure, got {l5}");

    assert_eq!(-176.503, round_to_decimals(b0, 3), "B0 failure, got {b0}");
    assert_eq!(3.068, round_to_decimals(b1, 3), "B1 failure, got {b1}");

    assert_eq!(
        99653849.038,
        round_to_decimals(r0, 3),
        "R0 failure, got {r0}"
    );
    assert_eq!(100378.567, round_to_decimals(r1, 3), "R1 failure, got {r1}");
    assert_eq!(-1140.954, round_to_decimals(r2, 3), "R2 failure, got {r2}");
    assert_eq!(-141.115, round_to_decimals(r3, 3), "R3 failure, got {r3}");
    assert_eq!(1.232, round_to_decimals(r4, 3), "R4 failure, got {r4}");

    let l = calculate_heliocentric_coeff(jm, l0, l1, l2, l3, l4, l5, CoordType::Longitude);
    let b = calculate_heliocentric_coeff(jm, b0, b1, 0.0, 0.0, 0.0, 0.0, CoordType::Latitude);
    let r = calculate_heliocentric_coeff(jm, r0, r1, r2, r3, r4, 0.0, CoordType::Radius);

    assert_eq!(
        24.018,
        round_to_decimals(l.to_degrees(), 3),
        "Heliocentric longitude failure, got {l}"
    );
    assert_eq!(
        -0.000101,
        round_to_decimals(b.to_degrees(), 6),
        "Heliocentric latitude failure, got {b}"
    );
    assert_eq!(
        0.9965423,
        round_to_decimals(r, 7),
        "Earth radius vector failure, got {r}"
    );

    let theta = calculate_geocentric_coeff(l, CoordType::Longitude);
    let beta = calculate_geocentric_coeff(b, CoordType::Latitude);

    assert_eq!(
        204.0182617,
        round_to_decimals(theta.to_degrees(), 7),
        "Geocentric longitude failure, got {theta}"
    );
    assert_eq!(
        0.0001011219,
        round_to_decimals(beta.to_degrees(), 10),
        "Geocentric latitude failure, got {beta}"
    );

    let (dpsi, deps) = calculate_dpsi_depsilon(jc);
    let eps = calculate_epsilon(jm, deps);

    assert_eq!(
        -0.003998,
        round_to_decimals(dpsi.to_degrees(), 6),
        "delta psi failure, got {dpsi}"
    );
    assert_eq!(
        0.001667,
        round_to_decimals(deps.to_degrees(), 6),
        "delta epsilon failure, got {deps}"
    );
    assert_eq!(
        23.440465,
        round_to_decimals(eps.to_degrees(), 6),
        "epsilon failure, got {eps}"
    );

    let dtau = aberration_correction(r);
    let lambda = apparent_sun_longitude(theta, dpsi, dtau);

    assert_eq!(
        204.008552,
        round_to_decimals(lambda.to_degrees(), 6),
        "Apparent sun longitude failure, got {lambda}",
    );

    let alpha = geocentric_sun_right_ascension(lambda, eps, beta);
    assert_eq!(
        202.22741,
        round_to_decimals(alpha.to_degrees(), 5),
        "Geocentric sun right ascension, got {alpha}",
    );

    let delta = geocentric_sun_declination(lambda, eps, beta);
    assert_eq!(
        -9.31434,
        round_to_decimals(delta.to_degrees(), 5),
        "geocentric sun declination failure, got {delta}"
    );

    let v = sidereal_time_greenwich(jd_ut, jc, dpsi, eps);
    let hour_angle = obs_local_hour_angle(lon_rad, v, alpha);
    assert_eq!(
        11.1059,
        round_to_decimals(hour_angle.to_degrees(), 4),
        "Observer local hour angle failure, got {hour_angle}"
    );

    let (delta_, hour_angle_) = topocentric_sun_right_ascension_declination_and_local_hour_angle(
        lat_rad, 1830.14, r, hour_angle, delta,
    );
    assert_eq!(
        11.1063,
        round_to_decimals(hour_angle_.to_degrees(), 4),
        "Topocentric local hour angle failure, got {hour_angle_}"
    );
    assert_eq!(
        -9.316179,
        round_to_decimals(delta_.to_degrees(), 6),
        "Topocentric sun right ascension declination failure, got {delta_}"
    );

    let temperature = Some(11.0);
    let pressure = Some(820.0);

    let e0 = topocentric_elevation_angle(lat_rad, delta_, hour_angle_, pressure, temperature);
    let zenith = topocentric_zenith_angle(e0.to_degrees());
    let azimuth_w_s = topocentric_azimuth_angle_w_from_s(lat_rad, hour_angle_, delta_);
    let azimuth = topocentric_azimuth_angle_e_from_n(azimuth_w_s.to_degrees());
    let azimuth2 = topocentric_azimuth(lat, hour_angle_.to_degrees(), delta_.to_degrees());

    assert_eq!(
        50.11162,
        round_to_decimals(zenith, 5),
        "Topocentric zenith angle failure, got {zenith}"
    );
    assert_eq!(
        194.34024,
        round_to_decimals(azimuth, 5),
        "Topocentric azimuth failure, got {azimuth}"
    );
    assert_eq!(
        azimuth, azimuth2,
        "`topocentric_azimuth` failure, got {azimuth2}"
    );
}
