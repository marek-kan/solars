use crate::core::time::{calc_julian_day, julian_century, julian_millennium};
use crate::core::utils::round_to_decimals;
use crate::{periodic_tables::earth::*, periodic_tables::nutation::*};
use chrono::{DateTime, TimeZone, Utc};

fn get_julian_date_values() -> (f64, f64, f64) {
    let test_date: DateTime<Utc> = Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap();

    let mut jd = calc_julian_day(&test_date);
    assert_eq!(2452930.313, round_to_decimals(jd, 3), "Julian day failure");

    jd += 67.0 / 86400.0; // In the example they mention dT

    let jc = julian_century(jd);
    let jm = julian_millennium(jc);

    (jd, jc, jm)
}

#[test]
fn sum_tables() {
    let (_, _, jm) = get_julian_date_values();

    let l0_res = round_to_decimals(sum_table(&L0_TABLE, &jm), 3);
    let l1_res = round_to_decimals(sum_table(&L1_TABLE, &jm), 3);
    let l2_res = round_to_decimals(sum_table(&L2_TABLE, &jm), 3);
    let l3_res = round_to_decimals(sum_table(&L3_TABLE, &jm), 3);
    let l4_res = round_to_decimals(sum_table(&L4_TABLE, &jm), 3);
    let l5_res = round_to_decimals(sum_table(&L5_TABLE, &jm), 3);
    let b0_res = round_to_decimals(sum_table(&B0_TABLE, &jm), 3);
    let b1_res = round_to_decimals(sum_table(&B1_TABLE, &jm), 3);
    let r0_res = round_to_decimals(sum_table(&R0_TABLE, &jm), 3);
    let r1_res = round_to_decimals(sum_table(&R1_TABLE, &jm), 3);
    let r2_res = round_to_decimals(sum_table(&R2_TABLE, &jm), 3);
    let r3_res = round_to_decimals(sum_table(&R3_TABLE, &jm), 3);
    let r4_res = round_to_decimals(sum_table(&R4_TABLE, &jm), 3);

    assert_eq!(172067561.527, l0_res, "L0 failure");
    assert_eq!(628332010650.051, l1_res, "L1 failure");
    assert_eq!(61368.682, l2_res, "L2 failure");
    assert_eq!(-26.903, l3_res, "L3 failure");
    assert_eq!(-121.280, l4_res, "L4 failure");
    assert_eq!(-1.000, l5_res, "L5 failure");

    assert_eq!(-176.503, b0_res, "B0 failure");
    assert_eq!(3.068, b1_res, "B1 failure");

    assert_eq!(99653849.038, r0_res, "R0 failure");
    assert_eq!(100378.567, r1_res, "R1 failure");
    assert_eq!(-1140.954, r2_res, "R2 failure");
    assert_eq!(-141.115, r3_res, "R3 failure");
    assert_eq!(1.232, r4_res, "R4 failure");
}

#[test]
fn dpsi_depsilon() {
    let (_, jc, _) = get_julian_date_values();

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
