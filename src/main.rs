use chrono::{DateTime, TimeZone, Utc};
use ndarray::Array1;
use ndarray::prelude::*;
use solars::core::time::calc_julian_day;

/// Outdated, do not take into account
fn main() {

    // ALL THIS IS OUTDATED, DO NOT TAKE INTO ACCOUNT!

    // let d = Utc.with_ymd_and_hms(2024, 10, 30, 12, 0, 0).unwrap();
    // let i = d.timestamp_micros();
    // let date = DateTime::from_timestamp_micros(i).unwrap();
    // let lat: Array1<f64> = arr1(&[50.1, 50.05]);
    // let lon = arr1(&[14.5, 14.55]);

    // println!("Original date: {}, reconstructed date: {}", d, date);

    // let julian_day = calc_julian_day(&date);
    // let julian_century = (&julian_day - 2451545.0) / 36525.0;
    // let julian_millennium = &julian_century / 10.0;

    // println!("Julian Day: {julian_day}");
    // println!("Julian century: {julian_century}");
    // println!("Julian millennium: {julian_millennium}");

    // let t = julian_century.clone();

    // // THIS IS JUST SOME SLOP IT GENERATED INSTEAD OF USING THE LOOKUP TABLE IN THE PDF FILE

    // // --- 1. Earth / Solar Orbit Geometry Constants (via Meeus / NREL SPA) ---
    // // Mean elongation of the Moon, Sun's mean anomaly, etc.
    // // Truncated variations give clean, precise coordinates:
    // let L0 = (280.46646 + 36000.76983 * t + 0.0003032 * t * t) % 360.0; // Mean Longitude
    // let m_anomaly = 357.52911 + 35999.05029 * t - 0.0001537 * t * t; // Mean Anomaly

    // // Sun Equation of the Center (C)
    // let m_rad = m_anomaly.to_radians();
    // let center_corr = (1.914602 - 0.004817 * t - 0.000014 * t * t) * m_rad.sin()
    //     + (0.019993 - 0.000101 * t) * (2.0 * m_rad).sin()
    //     + 0.000289 * (3.0 * m_rad).sin();

    // let true_long = L0 + center_corr; // Sun True Longitude

    // // Obliquity of the Ecliptic (e)
    // let epsilon0 = 23.0
    //     + 26.0 / 60.0
    //     + (21.448 - 46.8150 * t - 0.00059 * t * t + 0.001813 * t * t * t) / 3600.0;
    // let epsilon = epsilon0.to_radians();

    // // --- 2. Calculate Apparent Solar Equatorial Coordinates ---
    // let true_long_rad = true_long.to_radians();
    // let declination = (epsilon.sin() * true_long_rad.sin()).asin(); // Sun Declination (δ)

    // // Right Ascension (α)
    // let ra = (epsilon.cos() * true_long_rad.sin()).atan2(true_long_rad.cos());

    // // Greenwich Mean Sidereal Time (GMST) in degrees
    // let mut gmst = 280.46061837 + 360.98564736629 * (julian_day - 2451545.0) + 0.000387933 * t * t
    //     - (t * t * t) / 38710000.0;
    // gmst = (gmst % 360.0 + 360.0) % 360.0;

    // // --- 3. Compute Vectorized Coordinate Spaces Over Lat/Lon ---
    // // Zip and map over matching structural arrays via `ndarray`
    // let results = ndarray::Zip::from(&lat)
    //     .and(&lon)
    //     .map_collect(|&lat_val, &lon_val| {
    //         let lat_rad = lat_val.to_radians();

    //         // Local Hour Angle (H) = GMST + Longitude - Right Ascension
    //         let mut h_deg = gmst + lon_val - ra.to_degrees();
    //         h_deg = (h_deg % 360.0 + 360.0) % 360.0;
    //         let h_rad = h_deg.to_radians();

    //         // Solar Elevation Angle (e) math
    //         let sin_el =
    //             lat_rad.sin() * declination.sin() + lat_rad.cos() * declination.cos() * h_rad.cos();
    //         let el_rad = sin_el.clamp(-1.0, 1.0).asin();

    //         // Solar Zenith Angle (θ) = 90° - Elevation
    //         let zenith_deg = 90.0 - el_rad.to_degrees();

    //         // Solar Azimuth Angle (A) measured clockwise from North
    //         let y_az = h_rad.sin();
    //         let x_az = h_rad.cos() * lat_rad.sin() - declination.tan() * lat_rad.cos();
    //         // Convert to clockwise angle tracking North = 0°
    //         let mut azimuth_deg = y_az.atan2(x_az).to_degrees() + 180.0;
    //         azimuth_deg = (azimuth_deg % 360.0 + 360.0) % 360.0;

    //         (zenith_deg, azimuth_deg)
    //     });

    // // Unpack the matrix of tuples into individual matrices
    // let zenith = results.mapv(|(z, _)| z);
    // let azimuth = results.mapv(|(_, a)| a);

    // println!("\n--- Solar Position Matrices ---");
    // println!("Zenith Matrix (degrees):\n{:#?}", zenith);
    // println!("Azimuth Matrix (degrees from North):\n{:#?}", azimuth);
}
