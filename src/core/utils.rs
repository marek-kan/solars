pub(crate) fn round_to_decimals(num: f64, decimals: i32) -> f64 {
    let factor = 10.0_f64.powi(decimals);
    (factor * num).round() / factor
}

pub(crate) fn limit_deg_to_360(deg: f64) -> f64 {
    let f = (deg / 360.0).fract();

    if deg.ceil() >= 0.0 {
        360.0 * f
    } else {
        360.0 - 360.0 * f
    }
}
