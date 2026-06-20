use chrono::{DateTime, Datelike, TimeZone, Timelike, Utc};

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
}
