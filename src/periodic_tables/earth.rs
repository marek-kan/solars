pub(crate) enum CoordType {
    Longitude,
    Latitude,
    Radius,
}

pub(crate) struct EarthPeriodicTermRow {
    a: f64,
    b: f64,
    c: f64,
}

impl EarthPeriodicTermRow {
    fn calculate_term(&self, jme: &f64) -> f64 {
        self.a * (self.b + self.c * jme).cos()
    }
}

pub(crate) fn sum_table(table: &[EarthPeriodicTermRow], jme: &f64) -> f64 {
    table.iter().map(|row| row.calculate_term(jme)).sum()
}

/// Returns:\
/// `CoordType::Longitude` => degrees bounded to [0, 360]\
/// `CoordType::Latitude` => degrees bounded to [0, 360]\
/// `CoordType::Radius` => Astronomical Units
pub(crate) fn calculate_heliocentric_coeff(
    jme: f64,
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
    c5: f64,
    coord_type: CoordType,
) -> f64 {
    let coord_rad =
        (c0 + c1 * jme + c2 * jme.powi(2) + c3 * jme.powi(3) + c4 * jme.powi(4) + c5 * jme.powi(5))
            / 10.0_f64.powi(8);

    match coord_type {
        CoordType::Longitude => limit_deg_to_360(coord_rad.to_degrees()),
        CoordType::Latitude => limit_deg_to_360(coord_rad.to_degrees()),
        CoordType::Radius => coord_rad,
    }
}

/// Returns:\
/// `CoordType::Longitude` => degrees bounded to [0, 360]\
/// `CoordType::Latitude` => degrees bounded to [0, 360]\
/// `CoordType::Radius` => Undefined
pub(crate) fn calculate_geocentric_coeff(c: f64, coord_type: CoordType) -> f64 {
    match coord_type {
        CoordType::Longitude => limit_deg_to_360(c + 180.0),
        CoordType::Latitude => limit_deg_to_360(-1.0 * c),
        CoordType::Radius => {
            panic!("`CoordType::Radius` is undefined for geocentric calculation")
        }
    }
}

fn limit_deg_to_360(deg: f64) -> f64 {
    let f = (deg / 360.0).fract();

    if deg.ceil() >= 0.0 {
        360.0 * f
    } else {
        360.0 - 360.0 * f
    }
}

pub(crate) const L0_TABLE: [EarthPeriodicTermRow; 64] = [
    EarthPeriodicTermRow {
        a: 175347046.0,
        b: 0.0,
        c: 0.0,
    },
    EarthPeriodicTermRow {
        a: 3341656.0,
        b: 4.6692568,
        c: 6283.07585,
    },
    EarthPeriodicTermRow {
        a: 34894.0,
        b: 4.6261,
        c: 12566.1517,
    },
    EarthPeriodicTermRow {
        a: 3497.0,
        b: 2.7441,
        c: 5753.3849,
    },
    EarthPeriodicTermRow {
        a: 3418.0,
        b: 2.8289,
        c: 3.5231,
    },
    EarthPeriodicTermRow {
        a: 3136.0,
        b: 3.6277,
        c: 77713.7715,
    },
    EarthPeriodicTermRow {
        a: 2676.0,
        b: 4.4181,
        c: 7860.4194,
    },
    EarthPeriodicTermRow {
        a: 2343.0,
        b: 6.1352,
        c: 3930.2097,
    },
    EarthPeriodicTermRow {
        a: 1324.0,
        b: 0.7425,
        c: 11506.7698,
    },
    EarthPeriodicTermRow {
        a: 1273.0,
        b: 2.0371,
        c: 529.691,
    },
    EarthPeriodicTermRow {
        a: 1199.0,
        b: 1.1096,
        c: 1577.3435,
    },
    EarthPeriodicTermRow {
        a: 990.0,
        b: 5.233,
        c: 5884.927,
    },
    EarthPeriodicTermRow {
        a: 902.0,
        b: 2.045,
        c: 26.298,
    },
    EarthPeriodicTermRow {
        a: 857.0,
        b: 3.508,
        c: 398.149,
    },
    EarthPeriodicTermRow {
        a: 780.0,
        b: 1.179,
        c: 5223.694,
    },
    EarthPeriodicTermRow {
        a: 753.0,
        b: 2.533,
        c: 5507.553,
    },
    EarthPeriodicTermRow {
        a: 505.0,
        b: 4.583,
        c: 18849.228,
    },
    EarthPeriodicTermRow {
        a: 492.0,
        b: 4.205,
        c: 775.523,
    },
    EarthPeriodicTermRow {
        a: 357.0,
        b: 2.92,
        c: 0.067,
    },
    EarthPeriodicTermRow {
        a: 317.0,
        b: 5.849,
        c: 11790.629,
    },
    EarthPeriodicTermRow {
        a: 284.0,
        b: 1.899,
        c: 796.298,
    },
    EarthPeriodicTermRow {
        a: 271.0,
        b: 0.315,
        c: 10977.079,
    },
    EarthPeriodicTermRow {
        a: 243.0,
        b: 0.345,
        c: 5486.778,
    },
    EarthPeriodicTermRow {
        a: 206.0,
        b: 4.806,
        c: 2544.314,
    },
    EarthPeriodicTermRow {
        a: 205.0,
        b: 1.869,
        c: 5573.143,
    },
    EarthPeriodicTermRow {
        a: 202.0,
        b: 2.458,
        c: 6069.777,
    },
    EarthPeriodicTermRow {
        a: 156.0,
        b: 0.833,
        c: 213.299,
    },
    EarthPeriodicTermRow {
        a: 132.0,
        b: 3.411,
        c: 2942.463,
    },
    EarthPeriodicTermRow {
        a: 126.0,
        b: 1.083,
        c: 20.775,
    },
    EarthPeriodicTermRow {
        a: 115.0,
        b: 0.645,
        c: 0.98,
    },
    EarthPeriodicTermRow {
        a: 103.0,
        b: 0.636,
        c: 4694.003,
    },
    EarthPeriodicTermRow {
        a: 102.0,
        b: 0.976,
        c: 15720.839,
    },
    EarthPeriodicTermRow {
        a: 102.0,
        b: 4.267,
        c: 7.114,
    },
    EarthPeriodicTermRow {
        a: 99.0,
        b: 6.21,
        c: 2146.17,
    },
    EarthPeriodicTermRow {
        a: 98.0,
        b: 0.68,
        c: 155.42,
    },
    EarthPeriodicTermRow {
        a: 86.0,
        b: 5.98,
        c: 161000.69,
    },
    EarthPeriodicTermRow {
        a: 85.0,
        b: 1.3,
        c: 6275.96,
    },
    EarthPeriodicTermRow {
        a: 85.0,
        b: 3.67,
        c: 71430.7,
    },
    EarthPeriodicTermRow {
        a: 80.0,
        b: 1.81,
        c: 17260.15,
    },
    EarthPeriodicTermRow {
        a: 79.0,
        b: 3.04,
        c: 12036.46,
    },
    EarthPeriodicTermRow {
        a: 75.0,
        b: 1.76,
        c: 5088.63,
    },
    EarthPeriodicTermRow {
        a: 74.0,
        b: 3.5,
        c: 3154.69,
    },
    EarthPeriodicTermRow {
        a: 74.0,
        b: 4.68,
        c: 801.82,
    },
    EarthPeriodicTermRow {
        a: 70.0,
        b: 0.83,
        c: 9437.76,
    },
    EarthPeriodicTermRow {
        a: 62.0,
        b: 3.98,
        c: 8827.39,
    },
    EarthPeriodicTermRow {
        a: 61.0,
        b: 1.82,
        c: 7084.9,
    },
    EarthPeriodicTermRow {
        a: 57.0,
        b: 2.78,
        c: 6286.6,
    },
    EarthPeriodicTermRow {
        a: 56.0,
        b: 4.39,
        c: 14143.5,
    },
    EarthPeriodicTermRow {
        a: 56.0,
        b: 3.47,
        c: 6279.55,
    },
    EarthPeriodicTermRow {
        a: 52.0,
        b: 0.19,
        c: 12139.55,
    },
    EarthPeriodicTermRow {
        a: 52.0,
        b: 1.33,
        c: 1748.02,
    },
    EarthPeriodicTermRow {
        a: 51.0,
        b: 0.28,
        c: 5856.48,
    },
    EarthPeriodicTermRow {
        a: 49.0,
        b: 0.49,
        c: 1194.45,
    },
    EarthPeriodicTermRow {
        a: 41.0,
        b: 5.37,
        c: 8429.24,
    },
    EarthPeriodicTermRow {
        a: 41.0,
        b: 2.4,
        c: 19651.05,
    },
    EarthPeriodicTermRow {
        a: 39.0,
        b: 6.17,
        c: 10447.39,
    },
    EarthPeriodicTermRow {
        a: 37.0,
        b: 6.04,
        c: 10213.29,
    },
    EarthPeriodicTermRow {
        a: 37.0,
        b: 2.57,
        c: 1059.38,
    },
    EarthPeriodicTermRow {
        a: 36.0,
        b: 1.71,
        c: 2352.87,
    },
    EarthPeriodicTermRow {
        a: 36.0,
        b: 1.78,
        c: 6812.77,
    },
    EarthPeriodicTermRow {
        a: 33.0,
        b: 0.59,
        c: 17789.85,
    },
    EarthPeriodicTermRow {
        a: 30.0,
        b: 0.44,
        c: 83996.85,
    },
    EarthPeriodicTermRow {
        a: 30.0,
        b: 2.74,
        c: 1349.87,
    },
    EarthPeriodicTermRow {
        a: 25.0,
        b: 3.16,
        c: 4690.48,
    },
];

pub(crate) const L1_TABLE: [EarthPeriodicTermRow; 34] = [
    EarthPeriodicTermRow {
        a: 628331966747.0,
        b: 0.0,
        c: 0.0,
    },
    EarthPeriodicTermRow {
        a: 206059.0,
        b: 2.678235,
        c: 6283.07585,
    },
    EarthPeriodicTermRow {
        a: 4303.0,
        b: 2.6351,
        c: 12566.1517,
    },
    EarthPeriodicTermRow {
        a: 425.0,
        b: 1.59,
        c: 3.523,
    },
    EarthPeriodicTermRow {
        a: 119.0,
        b: 5.796,
        c: 26.298,
    },
    EarthPeriodicTermRow {
        a: 109.0,
        b: 2.966,
        c: 1577.344,
    },
    EarthPeriodicTermRow {
        a: 93.0,
        b: 2.59,
        c: 18849.23,
    },
    EarthPeriodicTermRow {
        a: 72.0,
        b: 1.14,
        c: 529.69,
    },
    EarthPeriodicTermRow {
        a: 68.0,
        b: 1.87,
        c: 398.15,
    },
    EarthPeriodicTermRow {
        a: 67.0,
        b: 4.41,
        c: 5507.55,
    },
    EarthPeriodicTermRow {
        a: 59.0,
        b: 2.89,
        c: 5223.69,
    },
    EarthPeriodicTermRow {
        a: 56.0,
        b: 2.17,
        c: 155.42,
    },
    EarthPeriodicTermRow {
        a: 45.0,
        b: 0.4,
        c: 796.3,
    },
    EarthPeriodicTermRow {
        a: 36.0,
        b: 0.47,
        c: 775.52,
    },
    EarthPeriodicTermRow {
        a: 29.0,
        b: 2.65,
        c: 7.11,
    },
    EarthPeriodicTermRow {
        a: 21.0,
        b: 5.34,
        c: 0.98,
    },
    EarthPeriodicTermRow {
        a: 19.0,
        b: 1.85,
        c: 5486.78,
    },
    EarthPeriodicTermRow {
        a: 19.0,
        b: 4.97,
        c: 213.3,
    },
    EarthPeriodicTermRow {
        a: 17.0,
        b: 2.99,
        c: 6275.96,
    },
    EarthPeriodicTermRow {
        a: 16.0,
        b: 0.03,
        c: 2544.31,
    },
    EarthPeriodicTermRow {
        a: 16.0,
        b: 1.43,
        c: 2146.17,
    },
    EarthPeriodicTermRow {
        a: 15.0,
        b: 1.21,
        c: 10977.08,
    },
    EarthPeriodicTermRow {
        a: 12.0,
        b: 2.83,
        c: 1748.02,
    },
    EarthPeriodicTermRow {
        a: 12.0,
        b: 3.26,
        c: 5088.63,
    },
    EarthPeriodicTermRow {
        a: 12.0,
        b: 5.27,
        c: 1194.45,
    },
    EarthPeriodicTermRow {
        a: 12.0,
        b: 2.08,
        c: 4694.0,
    },
    EarthPeriodicTermRow {
        a: 11.0,
        b: 0.77,
        c: 553.57,
    },
    EarthPeriodicTermRow {
        a: 10.0,
        b: 1.3,
        c: 6286.6,
    },
    EarthPeriodicTermRow {
        a: 10.0,
        b: 4.24,
        c: 1349.87,
    },
    EarthPeriodicTermRow {
        a: 9.0,
        b: 2.7,
        c: 242.73,
    },
    EarthPeriodicTermRow {
        a: 9.0,
        b: 5.64,
        c: 951.72,
    },
    EarthPeriodicTermRow {
        a: 8.0,
        b: 5.3,
        c: 2352.87,
    },
    EarthPeriodicTermRow {
        a: 6.0,
        b: 2.65,
        c: 9437.76,
    },
    EarthPeriodicTermRow {
        a: 6.0,
        b: 4.67,
        c: 4690.48,
    },
];

pub(crate) const L2_TABLE: [EarthPeriodicTermRow; 20] = [
    EarthPeriodicTermRow {
        a: 52919.0,
        b: 0.0,
        c: 0.0,
    },
    EarthPeriodicTermRow {
        a: 8720.0,
        b: 1.0721,
        c: 6283.0758,
    },
    EarthPeriodicTermRow {
        a: 309.0,
        b: 0.867,
        c: 12566.152,
    },
    EarthPeriodicTermRow {
        a: 27.0,
        b: 0.05,
        c: 3.52,
    },
    EarthPeriodicTermRow {
        a: 16.0,
        b: 5.19,
        c: 26.3,
    },
    EarthPeriodicTermRow {
        a: 16.0,
        b: 3.68,
        c: 155.42,
    },
    EarthPeriodicTermRow {
        a: 10.0,
        b: 0.76,
        c: 18849.23,
    },
    EarthPeriodicTermRow {
        a: 9.0,
        b: 2.06,
        c: 77713.77,
    },
    EarthPeriodicTermRow {
        a: 7.0,
        b: 0.83,
        c: 775.52,
    },
    EarthPeriodicTermRow {
        a: 5.0,
        b: 4.66,
        c: 1577.34,
    },
    EarthPeriodicTermRow {
        a: 4.0,
        b: 1.03,
        c: 7.11,
    },
    EarthPeriodicTermRow {
        a: 4.0,
        b: 3.44,
        c: 5573.14,
    },
    EarthPeriodicTermRow {
        a: 3.0,
        b: 5.14,
        c: 796.3,
    },
    EarthPeriodicTermRow {
        a: 3.0,
        b: 6.05,
        c: 5507.55,
    },
    EarthPeriodicTermRow {
        a: 3.0,
        b: 1.19,
        c: 242.73,
    },
    EarthPeriodicTermRow {
        a: 3.0,
        b: 6.12,
        c: 529.69,
    },
    EarthPeriodicTermRow {
        a: 3.0,
        b: 0.31,
        c: 398.15,
    },
    EarthPeriodicTermRow {
        a: 3.0,
        b: 2.28,
        c: 553.57,
    },
    EarthPeriodicTermRow {
        a: 2.0,
        b: 4.38,
        c: 5223.69,
    },
    EarthPeriodicTermRow {
        a: 2.0,
        b: 3.75,
        c: 0.98,
    },
];

pub(crate) const L3_TABLE: [EarthPeriodicTermRow; 7] = [
    EarthPeriodicTermRow {
        a: 289.0,
        b: 5.844,
        c: 6283.076,
    },
    EarthPeriodicTermRow {
        a: 35.0,
        b: 0.0,
        c: 0.0,
    },
    EarthPeriodicTermRow {
        a: 17.0,
        b: 5.49,
        c: 12566.15,
    },
    EarthPeriodicTermRow {
        a: 3.0,
        b: 5.2,
        c: 155.42,
    },
    EarthPeriodicTermRow {
        a: 1.0,
        b: 4.72,
        c: 3.52,
    },
    EarthPeriodicTermRow {
        a: 1.0,
        b: 5.3,
        c: 18849.23,
    },
    EarthPeriodicTermRow {
        a: 1.0,
        b: 5.97,
        c: 242.73,
    },
];

pub(crate) const L4_TABLE: [EarthPeriodicTermRow; 3] = [
    EarthPeriodicTermRow {
        a: 114.0,
        b: 3.142,
        c: 0.0,
    },
    EarthPeriodicTermRow {
        a: 8.0,
        b: 4.13,
        c: 6283.08,
    },
    EarthPeriodicTermRow {
        a: 1.0,
        b: 3.84,
        c: 12566.15,
    },
];

pub(crate) const L5_TABLE: [EarthPeriodicTermRow; 1] = [EarthPeriodicTermRow {
    a: 1.0,
    b: 3.14,
    c: 0.0,
}];

pub(crate) const B0_TABLE: [EarthPeriodicTermRow; 5] = [
    EarthPeriodicTermRow {
        a: 280.0,
        b: 3.199,
        c: 84334.662,
    },
    EarthPeriodicTermRow {
        a: 102.0,
        b: 5.422,
        c: 5507.553,
    },
    EarthPeriodicTermRow {
        a: 80.0,
        b: 3.88,
        c: 5223.69,
    },
    EarthPeriodicTermRow {
        a: 44.0,
        b: 3.7,
        c: 2352.87,
    },
    EarthPeriodicTermRow {
        a: 32.0,
        b: 4.0,
        c: 1577.34,
    },
];

pub(crate) const B1_TABLE: [EarthPeriodicTermRow; 2] = [
    EarthPeriodicTermRow {
        a: 9.0,
        b: 3.9,
        c: 5507.55,
    },
    EarthPeriodicTermRow {
        a: 6.0,
        b: 1.73,
        c: 5223.69,
    },
];

pub(crate) const R0_TABLE: [EarthPeriodicTermRow; 40] = [
    EarthPeriodicTermRow {
        a: 100013989.0,
        b: 0.0,
        c: 0.0,
    },
    EarthPeriodicTermRow {
        a: 1670700.0,
        b: 3.0984635,
        c: 6283.07585,
    },
    EarthPeriodicTermRow {
        a: 13956.0,
        b: 3.05525,
        c: 12566.1517,
    },
    EarthPeriodicTermRow {
        a: 3084.0,
        b: 5.1985,
        c: 77713.7715,
    },
    EarthPeriodicTermRow {
        a: 1628.0,
        b: 1.1739,
        c: 5753.3849,
    },
    EarthPeriodicTermRow {
        a: 1576.0,
        b: 2.8469,
        c: 7860.4194,
    },
    EarthPeriodicTermRow {
        a: 925.0,
        b: 5.453,
        c: 11506.77,
    },
    EarthPeriodicTermRow {
        a: 542.0,
        b: 4.564,
        c: 3930.21,
    },
    EarthPeriodicTermRow {
        a: 472.0,
        b: 3.661,
        c: 5884.927,
    },
    EarthPeriodicTermRow {
        a: 346.0,
        b: 0.964,
        c: 5507.553,
    },
    EarthPeriodicTermRow {
        a: 329.0,
        b: 5.9,
        c: 5223.694,
    },
    EarthPeriodicTermRow {
        a: 307.0,
        b: 0.299,
        c: 5573.143,
    },
    EarthPeriodicTermRow {
        a: 243.0,
        b: 4.273,
        c: 11790.629,
    },
    EarthPeriodicTermRow {
        a: 212.0,
        b: 5.847,
        c: 1577.344,
    },
    EarthPeriodicTermRow {
        a: 186.0,
        b: 5.022,
        c: 10977.079,
    },
    EarthPeriodicTermRow {
        a: 175.0,
        b: 3.012,
        c: 18849.228,
    },
    EarthPeriodicTermRow {
        a: 110.0,
        b: 5.055,
        c: 5486.778,
    },
    EarthPeriodicTermRow {
        a: 98.0,
        b: 0.89,
        c: 6069.78,
    },
    EarthPeriodicTermRow {
        a: 86.0,
        b: 5.69,
        c: 15720.84,
    },
    EarthPeriodicTermRow {
        a: 86.0,
        b: 1.27,
        c: 161000.69,
    },
    EarthPeriodicTermRow {
        a: 65.0,
        b: 0.27,
        c: 17260.15,
    },
    EarthPeriodicTermRow {
        a: 63.0,
        b: 0.92,
        c: 529.69,
    },
    EarthPeriodicTermRow {
        a: 57.0,
        b: 2.01,
        c: 83996.85,
    },
    EarthPeriodicTermRow {
        a: 56.0,
        b: 5.24,
        c: 71430.7,
    },
    EarthPeriodicTermRow {
        a: 49.0,
        b: 3.25,
        c: 2544.31,
    },
    EarthPeriodicTermRow {
        a: 47.0,
        b: 2.58,
        c: 775.52,
    },
    EarthPeriodicTermRow {
        a: 45.0,
        b: 5.54,
        c: 9437.76,
    },
    EarthPeriodicTermRow {
        a: 43.0,
        b: 6.01,
        c: 6275.96,
    },
    EarthPeriodicTermRow {
        a: 39.0,
        b: 5.36,
        c: 4694.0,
    },
    EarthPeriodicTermRow {
        a: 38.0,
        b: 2.39,
        c: 8827.39,
    },
    EarthPeriodicTermRow {
        a: 37.0,
        b: 0.83,
        c: 19651.05,
    },
    EarthPeriodicTermRow {
        a: 37.0,
        b: 4.9,
        c: 12139.55,
    },
    EarthPeriodicTermRow {
        a: 36.0,
        b: 1.67,
        c: 12036.46,
    },
    EarthPeriodicTermRow {
        a: 35.0,
        b: 1.84,
        c: 2942.46,
    },
    EarthPeriodicTermRow {
        a: 33.0,
        b: 0.24,
        c: 7084.9,
    },
    EarthPeriodicTermRow {
        a: 32.0,
        b: 0.18,
        c: 5088.63,
    },
    EarthPeriodicTermRow {
        a: 32.0,
        b: 1.78,
        c: 398.15,
    },
    EarthPeriodicTermRow {
        a: 28.0,
        b: 1.21,
        c: 6286.6,
    },
    EarthPeriodicTermRow {
        a: 28.0,
        b: 1.9,
        c: 6279.55,
    },
    EarthPeriodicTermRow {
        a: 26.0,
        b: 4.59,
        c: 10447.39,
    },
];

pub(crate) const R1_TABLE: [EarthPeriodicTermRow; 10] = [
    EarthPeriodicTermRow {
        a: 103019.0,
        b: 1.10749,
        c: 6283.07585,
    },
    EarthPeriodicTermRow {
        a: 1721.0,
        b: 1.0644,
        c: 12566.1517,
    },
    EarthPeriodicTermRow {
        a: 702.0,
        b: 3.142,
        c: 0.0,
    },
    EarthPeriodicTermRow {
        a: 32.0,
        b: 1.02,
        c: 18849.23,
    },
    EarthPeriodicTermRow {
        a: 31.0,
        b: 2.84,
        c: 5507.55,
    },
    EarthPeriodicTermRow {
        a: 25.0,
        b: 1.32,
        c: 5223.69,
    },
    EarthPeriodicTermRow {
        a: 18.0,
        b: 1.42,
        c: 1577.34,
    },
    EarthPeriodicTermRow {
        a: 10.0,
        b: 5.91,
        c: 10977.08,
    },
    EarthPeriodicTermRow {
        a: 9.0,
        b: 1.42,
        c: 6275.96,
    },
    EarthPeriodicTermRow {
        a: 9.0,
        b: 0.27,
        c: 5486.78,
    },
];

pub(crate) const R2_TABLE: [EarthPeriodicTermRow; 6] = [
    EarthPeriodicTermRow {
        a: 4359.0,
        b: 5.7846,
        c: 6283.0758,
    },
    EarthPeriodicTermRow {
        a: 124.0,
        b: 5.579,
        c: 12566.152,
    },
    EarthPeriodicTermRow {
        a: 12.0,
        b: 3.14,
        c: 0.0,
    },
    EarthPeriodicTermRow {
        a: 9.0,
        b: 3.63,
        c: 77713.77,
    },
    EarthPeriodicTermRow {
        a: 6.0,
        b: 1.87,
        c: 5573.14,
    },
    EarthPeriodicTermRow {
        a: 3.0,
        b: 5.47,
        c: 18849.23,
    },
];

pub(crate) const R3_TABLE: [EarthPeriodicTermRow; 2] = [
    EarthPeriodicTermRow {
        a: 145.0,
        b: 4.273,
        c: 6283.076,
    },
    EarthPeriodicTermRow {
        a: 7.0,
        b: 3.92,
        c: 12566.15,
    },
];

pub(crate) const R4_TABLE: [EarthPeriodicTermRow; 1] = [EarthPeriodicTermRow {
    a: 4.0,
    b: 2.56,
    c: 6283.08,
}];
