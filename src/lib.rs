mod periodic_tables;

use chrono::{DateTime, Datelike, TimeZone, Timelike, Utc};
use pyo3::prelude::*;

pub(crate) fn round_to_decimals(num: f64, decimals: i32) -> f64 {
    let factor = 10.0_f64.powi(decimals);
    (factor * num).round() / factor
}

pub fn calc_julian_day(date: &DateTime<Utc>) -> f64 {
    let gregorian_start: DateTime<Utc> = Utc.with_ymd_and_hms(1582, 10, 15, 0, 0, 0).unwrap();

    let mut y = date.year() as f64;
    let mut m = date.month() as f64;

    if m <= 2.0 {
        y -= 1.0;
        m += 12.0;
    }

    let d = date.day() as f64
        + (date.hour() as f64 / 24.0)
        + (date.minute() as f64 / 1440.0)
        + (date.second() as f64 / 86400.0);

    let b = if date >= &gregorian_start {
        let a = (y / 100.0).floor();
        2.0 - a + (a / 4.0).floor()
    } else {
        0.0
    };

    let year_part = (365.25 * (y + 4716.0)).floor();
    let month_part = (30.6001 * (m + 1.0)).floor();

    year_part + month_part + d + b - 1524.5
}

pub fn julian_century(julian_day: f64) -> f64 {
    (julian_day - 2451545.0) / 36525.0
}

pub fn julian_millennium(julian_century: f64) -> f64 {
    julian_century / 10.0
}

/// X_0 (degrees)
pub fn calculate_elongation_moon_from_sun(jce: f64) -> f64 {
    297.85036 + 445267.111480 * jce - 0.0019142 * jce.powi(2) + jce.powi(3) / 189474.0
}

/// X_1 (degrees)
pub fn calculate_sun_anomaly(jce: f64) -> f64 {
    357.52772 + 35999.050340 * jce - 0.0001603 * jce.powi(2) - jce.powi(3) / 30000.0
}

/// X_2 (degrees)
pub fn calculate_moon_anomaly(jce: f64) -> f64 {
    134.96298 + 477198.867398 * jce - 0.0086972 * jce.powi(2) - jce.powi(3) / 56250.0
}

/// X_3 (degrees)
pub fn calculate_moon_lat(jce: f64) -> f64 {
    93.27191 + 483202.017538 * jce - 0.0036825 * jce.powi(2) - jce.powi(3) / 327270.0
}

/// X_4 (degrees)
/// the longitude of the ascending node of the moon’s mean orbit on the ecliptic, measured from the mean equinox of the date
pub fn calculate_moon_lon(jce: f64) -> f64 {
    125.04452 - 1934.136261 * jce - 0.0020708 * jce.powi(2) - jce.powi(3) / 450000.0
}

/// A Python module implemented in Rust.
#[pymodule]
mod solars {
    use pyo3::prelude::*;

    /// Formats the sum of two numbers as string.
    #[pyfunction]
    pub fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
        Ok((a + b).to_string())
    }
}

#[cfg(test)]
mod test {
    use chrono::{DateTime, TimeZone, Utc};

    use crate::{calc_julian_day, round_to_decimals};

    #[test]
    fn test_julian_day() {
        let test_cases: [DateTime<Utc>; 5] = [
            Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap(),
            Utc.with_ymd_and_hms(1987, 1, 27, 0, 0, 0).unwrap(),
            Utc.with_ymd_and_hms(1987, 6, 19, 12, 0, 0).unwrap(),
            Utc.with_ymd_and_hms(1988, 6, 19, 12, 0, 0).unwrap(),
            Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap(),
        ];

        let test_results: [f64; 5] = [2451545.0, 2446822.5, 2446966.0, 2447332.0, 2452930.313];

        let iter = test_cases.iter().zip(test_results);

        for (i, (date, result)) in iter.enumerate() {
            let calculated_result = round_to_decimals(calc_julian_day(date), 3);

            assert_eq!(
                result, calculated_result,
                "Calculated julian day not same for {i}"
            );
        }
    }
}
