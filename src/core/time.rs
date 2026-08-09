use chrono::{DateTime, Datelike, TimeZone, Timelike, Utc};

pub(crate) fn calc_julian_day(date: &DateTime<Utc>) -> f64 {
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

pub(crate) fn julian_century(julian_day: f64) -> f64 {
    (julian_day - 2451545.0) / 36525.0
}

pub(crate) fn julian_millennium(julian_century: f64) -> f64 {
    julian_century / 10.0
}

/// Estimates Delta T in seconds using the Espenak-Meeus NASA polynomials.
/// Ref: https://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html
pub(crate) fn delta_t_seconds(date: &DateTime<Utc>) -> f64 {
    let y = date.year() as f64 + (date.month() as f64 - 0.5) / 12.0;
    let delta_t = if y < -500.0 {
        let u = (y - 1820.0) / 100.0;
        -20.0 + 32.0 * u.powi(2)
    } else if y < 500.0 {
        let u = y / 100.0;
        10583.6 - 1014.41 * u + 33.78311 * u.powi(2) - 5.952053 * u.powi(3) - 0.1798452 * u.powi(4)
            + 0.022174192 * u.powi(5)
            + 0.0090316521 * u.powi(6)
    } else if y < 1600.0 {
        let u = (y - 1000.0) / 100.0;
        1574.2 - 556.01 * u + 71.23472 * u.powi(2) + 0.319781 * u.powi(3)
            - 0.8503463 * u.powi(4)
            - 0.005050998 * u.powi(5)
            + 0.0083572073 * u.powi(6)
    } else if y < 1700.0 {
        let t = y - 1600.0;
        120.0 - 0.9808 * t - 0.01532 * t.powi(2) + t.powi(3) / 7129.0
    } else if y < 1800.0 {
        let t = y - 1700.0;
        8.83 + 0.1603 * t - 0.0059285 * t.powi(2) + 0.00013336 * t.powi(3) - t.powi(4) / 1174000.0
    } else if y < 1860.0 {
        let t = y - 1800.0;
        13.72 - 0.332447 * t + 0.0068612 * t.powi(2) + 0.0041116 * t.powi(3)
            - 0.00037436 * t.powi(4)
            + 0.0000121272 * t.powi(5)
            - 0.0000001699 * t.powi(6)
            + 0.000000000875 * t.powi(7)
    } else if y < 1900.0 {
        let t = y - 1860.0;
        7.62 + 0.5737 * t - 0.251754 * t.powi(2) + 0.01680668 * t.powi(3) - 0.0004473624 * t.powi(4)
            + t.powi(5) / 233174.0
    } else if y < 1920.0 {
        let t = y - 1900.0;
        -2.79 + 1.494119 * t - 0.0598939 * t.powi(2) + 0.0061966 * t.powi(3) - 0.000197 * t.powi(4)
    } else if y < 1941.0 {
        let t = y - 1920.0;
        21.20 + 0.84493 * t - 0.076100 * t.powi(2) + 0.0020936 * t.powi(3)
    } else if y < 1961.0 {
        let t = y - 1950.0;
        29.07 + 0.407 * t - t.powi(2) / 233.0 + t.powi(3) / 2547.0
    } else if y < 1986.0 {
        let t = y - 1975.0;
        45.45 + 1.067 * t - t.powi(2) / 260.0 - t.powi(3) / 718.0
    } else if y < 2005.0 {
        let t = y - 2000.0;
        63.86 + 0.3345 * t - 0.060374 * t.powi(2)
            + 0.0017275 * t.powi(3)
            + 0.000651814 * t.powi(4)
            + 0.00002373599 * t.powi(5)
    } else if y < 2050.0 {
        let t = y - 2000.0;
        62.92 + 0.32217 * t + 0.005589 * t.powi(2)
    } else if y < 2150.0 {
        let u = (y - 1820.0) / 100.0;
        -20.0 + 32.0 * u.powi(2) - 0.5628 * (2150.0 - y)
    } else {
        let u = (y - 1820.0) / 100.0;
        -20.0 + 32.0 * u.powi(2)
    };

    if (1955.0..=2005.0).contains(&y) {
        delta_t
    } else {
        delta_t - 0.000012932 * (y - 1955.0).powi(2)
    }
}

#[cfg(test)]
mod test {
    use chrono::{DateTime, TimeZone, Utc};

    use super::*;
    use crate::core::utils::round_to_decimals;

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

    #[test]
    fn delta_t_uses_decimal_year_and_modern_nasa_polynomial() {
        let date1: DateTime<Utc> = Utc.with_ymd_and_hms(2000, 1, 1, 0, 0, 0).unwrap();
        let date2: DateTime<Utc> = Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap();

        let dt2 = round_to_decimals(delta_t_seconds(&date2), 2);

        assert!((delta_t_seconds(&date1) - 63.8738328).abs() < 1e-6);
        assert!((dt2 - 64.51).abs() < 1e-6);
    }
}
