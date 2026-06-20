pub(crate) fn round_to_decimals(num: f64, decimals: i32) -> f64 {
    let factor = 10.0_f64.powi(decimals);
    (factor * num).round() / factor
}
