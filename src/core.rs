use std::collections::HashMap;

use chrono::{DateTime, Utc};

pub mod grid;
pub mod irradiance;
pub(crate) mod solar_position;
pub(crate) mod time;
pub(crate) mod utils;

pub use grid::{
    AoiInput, AoiResult, AtmosphericInput, ClearSkyInput, ClearSkyResult, PoaInput, PoaResult,
    SolarError, SolarPositionInput, SolarPositionResult, SpatialInput, calculate_aoi,
    calculate_clearsky, calculate_poa, calculate_solar_position,
};

use crate::core::irradiance::aoi;
use crate::core::solar_position::{
    aberration_correction, apparent_sun_longitude, geocentric_sun_declination,
    geocentric_sun_right_ascension, obs_local_hour_angle, sidereal_time_greenwich,
    topocentric_azimuth_angle_e_from_n, topocentric_azimuth_angle_w_from_s,
    topocentric_elevation_angle, topocentric_sun_right_ascension_declination_and_local_hour_angle,
    topocentric_zenith_angle,
};
use crate::core::time::{calc_julian_day, delta_t_seconds, julian_century, julian_millennium};
use crate::periodic_tables::earth::{
    B0_TABLE, B1_TABLE, CoordType, L0_TABLE, L1_TABLE, L2_TABLE, L3_TABLE, L4_TABLE, L5_TABLE,
    R0_TABLE, R1_TABLE, R2_TABLE, R3_TABLE, R4_TABLE, calculate_geocentric_coeff,
    calculate_heliocentric_coeff, sum_table,
};
use crate::periodic_tables::nutation::{calculate_dpsi_depsilon, calculate_epsilon};

/// Calculates solar position for one location and UTC instant.
///
/// Latitude, longitude, panel tilt, and panel azimuth are degrees.
/// Elevation is meters,
/// pressure is millibars
/// temperature is degrees Celsius.
/// panel tilt is measured from horizontal
/// panel azimuth is measured eastward from north; 180 means south facing
///
/// Result:
/// The returned `zenith`, `azimuth`, and optional `aoi` are in degrees.
/// Missing elevation defaults to sea level. Atmospheric refraction is
/// applied only when both pressure and temperature are provided. `aoi` is only
/// included when both panel values are provided.
pub fn calculate_scalar_solar_position(
    latitude: f64,
    longitude: f64,
    time: DateTime<Utc>,
    elevation: Option<f64>,
    temperature: Option<f64>,
    pressure: Option<f64>,
    panel_tilt: Option<f64>,
    panel_azimuth: Option<f64>,
) -> HashMap<String, f64> {
    let julian_day = calc_julian_day(&time);
    let julian_century_ut = julian_century(julian_day);
    let julian_day_tt = julian_day + delta_t_seconds(&time) / 86400.0;
    let julian_millennium_tt = julian_millennium(julian_century(julian_day_tt));

    let latitude_radians = latitude.to_radians();
    let longitude_radians = longitude.to_radians();
    let elevation = elevation.unwrap_or(0.0);

    let l0 = sum_table(&L0_TABLE, &julian_millennium_tt);
    let l1 = sum_table(&L1_TABLE, &julian_millennium_tt);
    let l2 = sum_table(&L2_TABLE, &julian_millennium_tt);
    let l3 = sum_table(&L3_TABLE, &julian_millennium_tt);
    let l4 = sum_table(&L4_TABLE, &julian_millennium_tt);
    let l5 = sum_table(&L5_TABLE, &julian_millennium_tt);
    let b0 = sum_table(&B0_TABLE, &julian_millennium_tt);
    let b1 = sum_table(&B1_TABLE, &julian_millennium_tt);
    let r0 = sum_table(&R0_TABLE, &julian_millennium_tt);
    let r1 = sum_table(&R1_TABLE, &julian_millennium_tt);
    let r2 = sum_table(&R2_TABLE, &julian_millennium_tt);
    let r3 = sum_table(&R3_TABLE, &julian_millennium_tt);
    let r4 = sum_table(&R4_TABLE, &julian_millennium_tt);

    let heliocentric_longitude = calculate_heliocentric_coeff(
        julian_millennium_tt,
        l0,
        l1,
        l2,
        l3,
        l4,
        l5,
        CoordType::Longitude,
    );
    let heliocentric_latitude = calculate_heliocentric_coeff(
        julian_millennium_tt,
        b0,
        b1,
        0.0,
        0.0,
        0.0,
        0.0,
        CoordType::Latitude,
    );
    let radius = calculate_heliocentric_coeff(
        julian_millennium_tt,
        r0,
        r1,
        r2,
        r3,
        r4,
        0.0,
        CoordType::Radius,
    );
    let geocentric_longitude =
        calculate_geocentric_coeff(heliocentric_longitude, CoordType::Longitude);
    let geocentric_latitude =
        calculate_geocentric_coeff(heliocentric_latitude, CoordType::Latitude);
    let (delta_psi, delta_epsilon) = calculate_dpsi_depsilon(julian_century(julian_day_tt));
    let epsilon = calculate_epsilon(julian_millennium_tt, delta_epsilon);
    let apparent_longitude = apparent_sun_longitude(
        geocentric_longitude,
        delta_psi,
        aberration_correction(radius),
    );
    let right_ascension =
        geocentric_sun_right_ascension(apparent_longitude, epsilon, geocentric_latitude);
    let declination = geocentric_sun_declination(apparent_longitude, epsilon, geocentric_latitude);
    let sidereal_time = sidereal_time_greenwich(julian_day, julian_century_ut, delta_psi, epsilon);

    let hour_angle = obs_local_hour_angle(longitude_radians, sidereal_time, right_ascension);
    let (topocentric_declination, topocentric_hour_angle) =
        topocentric_sun_right_ascension_declination_and_local_hour_angle(
            latitude_radians,
            elevation,
            radius,
            hour_angle,
            declination,
        );
    let elevation_angle = topocentric_elevation_angle(
        latitude_radians,
        topocentric_declination,
        topocentric_hour_angle,
        pressure,
        temperature,
    );
    let zenith = topocentric_zenith_angle(elevation_angle.to_degrees());
    let azimuth_west_from_south = topocentric_azimuth_angle_w_from_s(
        latitude_radians,
        topocentric_hour_angle,
        topocentric_declination,
    );
    let azimuth = topocentric_azimuth_angle_e_from_n(azimuth_west_from_south.to_degrees());

    let mut result = HashMap::from([
        ("latitude".to_owned(), latitude),
        ("longitude".to_owned(), longitude),
        ("zenith".to_owned(), zenith),
        ("azimuth".to_owned(), azimuth),
    ]);

    if let (Some(tilt), Some(panel_direction)) = (panel_tilt, panel_azimuth) {
        result.insert(
            "aoi".to_owned(),
            aoi(zenith, azimuth, tilt, panel_direction, None),
        );
    }

    result
}

#[cfg(test)]
mod tests {
    use super::calculate_scalar_solar_position;
    use crate::core::utils::round_to_decimals;
    use chrono::{TimeZone, Utc};

    #[test]
    fn solar_position_omits_aoi_without_complete_panel_geometry() {
        let time = Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap();
        let result = calculate_scalar_solar_position(
            39.742476,
            -105.1786,
            time,
            None,
            None,
            None,
            Some(30.0),
            None,
        );

        assert_eq!(result.len(), 4);
        assert_eq!(result["latitude"], 39.742476);
        assert_eq!(result["longitude"], -105.1786);
        assert!((0.0..=180.0).contains(&result["zenith"]));
        assert!((0.0..360.0).contains(&result["azimuth"]));
        assert!(!result.contains_key("aoi"));
    }

    #[test]
    fn solar_position_includes_aoi_with_complete_panel_geometry() {
        let time = Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap();
        let result = calculate_scalar_solar_position(
            39.742476,
            -105.1786,
            time,
            Some(1830.14),
            Some(11.0),
            Some(820.0),
            Some(30.0),
            Some(180.0),
        );

        let zenith = result.get("zenith").unwrap().to_owned();
        let azimuth = result.get("azimuth").unwrap().to_owned();

        assert!(
            (round_to_decimals(zenith, 3) - 50.112).abs() < 1e-6,
            "Zenith failure, got: {zenith}"
        );
        assert!(
            (round_to_decimals(azimuth, 3) - 194.340).abs() < 1e-6,
            "Azimuth failure, got: {azimuth}"
        );
        assert_eq!(result.len(), 5);
        assert!((0.0..=180.0).contains(&result["aoi"]));
    }
}
