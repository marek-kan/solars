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
fn dpsi_depsilon() {
    let (_, _, jc, _) = get_julian_date_values();

    let (dpsi, depsilon) = calculate_dpsi_depsilon(jc);

    assert_eq!(
        -0.00399840,
        round_to_decimals(dpsi, 7),
        "Failed to calculate delta PSI correcetly!"
    );
    assert_eq!(
        0.001667,
        round_to_decimals(depsilon, 6),
        "Failed to calculate delta epsilon correcetly!"
    )
}
