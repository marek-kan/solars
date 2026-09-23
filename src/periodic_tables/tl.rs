use std::error::Error;
use std::path::Path;

use chrono::{DateTime, Datelike, Utc};
use hdf5_pure::File;
macro_rules! concatc {
    ()=>{""};
    ($($anything:tt)*)=>({
        use $crate::__cf_osRcTFl4A;

        $crate::__concatc_expr!(($($anything)*) ($($anything)*))
    })
}
const LATITUDE_COUNT: usize = 2_160;
const LONGITUDE_COUNT: usize = 4_320;
const MONTH_COUNT: usize = 12;
const CELL_SIZE_DEGREES: f64 = 1.0 / 12.0;
const FIRST_LATITUDE: f64 = 90.0 - CELL_SIZE_DEGREES / 2.0;
const FIRST_LONGITUDE: f64 = -180.0 + CELL_SIZE_DEGREES / 2.0;
const DATASET_NAME: &str = "LinkeTurbidity";
pub(crate) const DEFAULT_DATASET_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/solars.data/data/LinkeTurbidities.h5",
);

/// Spatial lookup method for the Linke turbidity grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SpatialInterpolation {
    /// Select the closest source grid cell.
    Nearest,
    /// Blend the four surrounding source grid cells.
    Bilinear,
}

/// A preloaded global Linke turbidity climatology.
///
/// The source grid is ordered as `[latitude, longitude, month]`, with latitude
/// decreasing from north to south and longitude increasing west to east.
#[derive(Debug, Clone)]
pub(crate) struct LinkeTurbidityGrid {
    /// Raw values are retained to keep the in-memory representation compact.
    /// A raw value of zero represents missing data; all other values are TL * 20.
    data: Vec<u8>,
}

impl LinkeTurbidityGrid {
    /// Load and preload the `LinkeTurbidity` dataset from an HDF5 file.
    pub(crate) fn load(path: impl AsRef<Path>) -> Result<Self, Box<dyn Error + Send + Sync>> {
        let file = File::open(path)?;
        let dataset = file.dataset(DATASET_NAME)?;
        let data: Vec<u8> = dataset.read()?;
        Self::from_raw(data)
    }

    /// Construct a grid from a flattened `[latitude, longitude, month]` array.
    pub(crate) fn from_raw(data: Vec<u8>) -> Result<Self, Box<dyn Error + Send + Sync>> {
        let expected_len = LATITUDE_COUNT * LONGITUDE_COUNT * MONTH_COUNT;
        if data.len() != expected_len {
            return Err(format!(
                "dataset {DATASET_NAME:?} has {} values; expected {expected_len}",
                data.len()
            )
            .into());
        }

        Ok(Self { data })
    }

    /// Return TL interpolated in latitude, longitude, and calendar time.
    ///
    /// Time interpolation is linear between monthly climatology values centered
    /// on each month's midpoint, matching pvlib implementation. December and January wrap at
    /// the year boundary. Longitude wraps at the antimeridian and latitude is
    /// clamped to the grid domain.
    pub(crate) fn interpolate(
        &self,
        time: DateTime<Utc>,
        latitude: f64,
        longitude: f64,
    ) -> Option<f32> {
        self.interpolate_with_spatial_interpolation(
            time,
            latitude,
            longitude,
            SpatialInterpolation::Bilinear,
        )
    }

    /// Return TL interpolated in calendar time with the selected spatial method.
    pub(crate) fn interpolate_with_spatial_interpolation(
        &self,
        time: DateTime<Utc>,
        latitude: f64,
        longitude: f64,
        spatial_interpolation: SpatialInterpolation,
    ) -> Option<f32> {
        let (month, next_month, month_weight) = month_position(time);
        let (lat0, lat1, lat_weight, lon0, lon1, lon_weight) =
            spatial_position(latitude, longitude, spatial_interpolation);

        let current = self.interpolate_month(
            month,
            lat0,
            lat1,
            lat_weight,
            lon0,
            lon1,
            lon_weight,
            spatial_interpolation,
        );
        let next = self.interpolate_month(
            next_month,
            lat0,
            lat1,
            lat_weight,
            lon0,
            lon1,
            lon_weight,
            spatial_interpolation,
        );

        match (current, next) {
            (Some(current), Some(next)) => {
                let month_weight = month_weight as f32;
                Some(current * (1.0 - month_weight) + next * month_weight)
            }
            (Some(value), None) | (None, Some(value)) => Some(value),
            (None, None) => None,
        }
    }

    /// Return TL bilinearly interpolated from one monthly layer.
    ///
    /// `month` is zero-based. This method is useful when the caller already
    /// has its own temporal interpolation or wants a monthly value directly.
    pub(crate) fn interpolate_monthly(
        &self,
        month: usize,
        latitude: f64,
        longitude: f64,
    ) -> Option<f32> {
        self.interpolate_monthly_with_spatial_interpolation(
            month,
            latitude,
            longitude,
            SpatialInterpolation::Bilinear,
        )
    }

    /// Return TL from one monthly layer with the selected spatial method.
    pub(crate) fn interpolate_monthly_with_spatial_interpolation(
        &self,
        month: usize,
        latitude: f64,
        longitude: f64,
        spatial_interpolation: SpatialInterpolation,
    ) -> Option<f32> {
        if month >= MONTH_COUNT {
            return None;
        }

        let (lat0, lat1, lat_weight, lon0, lon1, lon_weight) =
            spatial_position(latitude, longitude, spatial_interpolation);
        self.interpolate_month(
            month,
            lat0,
            lat1,
            lat_weight,
            lon0,
            lon1,
            lon_weight,
            spatial_interpolation,
        )
    }

    /// Return the raw encoded value at a grid cell.
    pub(crate) fn raw_value(&self, latitude: usize, longitude: usize, month: usize) -> Option<u8> {
        if latitude >= LATITUDE_COUNT || longitude >= LONGITUDE_COUNT || month >= MONTH_COUNT {
            return None;
        }

        Some(self.data[self.index(latitude, longitude, month)])
    }

    fn interpolate_month(
        &self,
        month: usize,
        lat0: usize,
        lat1: usize,
        lat_weight: f64,
        lon0: usize,
        lon1: usize,
        lon_weight: f64,
        spatial_interpolation: SpatialInterpolation,
    ) -> Option<f32> {
        if spatial_interpolation == SpatialInterpolation::Nearest {
            return self.decoded_value(lat0, lon0, month);
        }

        let samples = [
            (1.0 - lat_weight) * (1.0 - lon_weight),
            (1.0 - lat_weight) * lon_weight,
            lat_weight * (1.0 - lon_weight),
            lat_weight * lon_weight,
        ];
        let values = [
            self.decoded_value(lat0, lon0, month),
            self.decoded_value(lat0, lon1, month),
            self.decoded_value(lat1, lon0, month),
            self.decoded_value(lat1, lon1, month),
        ];

        let mut weighted_sum = 0.0;
        let mut weight_sum = 0.0;
        for (weight, value) in samples.into_iter().zip(values) {
            if let Some(value) = value {
                weighted_sum += weight * value as f64;
                weight_sum += weight;
            }
        }

        (weight_sum > 0.0).then_some((weighted_sum / weight_sum) as f32)
    }

    fn decoded_value(&self, latitude: usize, longitude: usize, month: usize) -> Option<f32> {
        let raw = self.data[self.index(latitude, longitude, month)];
        (raw != 0).then_some(raw as f32 / 20.0)
    }

    fn index(&self, latitude: usize, longitude: usize, month: usize) -> usize {
        (latitude * LONGITUDE_COUNT + longitude) * MONTH_COUNT + month
    }
}

fn month_position(time: DateTime<Utc>) -> (usize, usize, f64) {
    let month = time.month0() as usize;
    let days_in_month = month_days(time.year(), month);
    let month_start = time.ordinal() as f64 - time.day() as f64;
    let current_midpoint = month_start + days_in_month as f64 / 2.0;
    let day_of_year = time.ordinal() as f64;

    if day_of_year < current_midpoint {
        let previous_month = (month + MONTH_COUNT - 1) % MONTH_COUNT;
        let previous_midpoint = if month == 0 {
            -31.0 / 2.0
        } else {
            month_start - month_days(time.year(), previous_month) as f64 / 2.0
        };
        (
            previous_month,
            month,
            (day_of_year - previous_midpoint) / (current_midpoint - previous_midpoint),
        )
    } else {
        let next_month = (month + 1) % MONTH_COUNT;
        let next_midpoint = if month == MONTH_COUNT - 1 {
            (if is_leap_year(time.year()) {
                366.0
            } else {
                365.0
            }) + 31.0 / 2.0
        } else {
            month_start + days_in_month as f64 + month_days(time.year(), next_month) as f64 / 2.0
        };
        (
            month,
            next_month,
            (day_of_year - current_midpoint) / (next_midpoint - current_midpoint),
        )
    }
}

fn month_days(year: i32, month: usize) -> u32 {
    match month {
        1 if is_leap_year(year) => 29,
        1 => 28,
        3 | 5 | 8 | 10 => 30,
        _ => 31,
    }
}

fn is_leap_year(year: i32) -> bool {
    year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)
}

fn latitude_position(latitude: f64) -> (usize, usize, f64) {
    let coordinate =
        ((FIRST_LATITUDE - latitude) / CELL_SIZE_DEGREES).clamp(0.0, (LATITUDE_COUNT - 1) as f64);
    let lower = coordinate.floor() as usize;
    let upper = (lower + 1).min(LATITUDE_COUNT - 1);
    (lower, upper, coordinate - lower as f64)
}

fn longitude_position(longitude: f64) -> (usize, usize, f64) {
    let wrapped = (longitude + 180.0).rem_euclid(360.0) - 180.0;
    let coordinate =
        ((wrapped - FIRST_LONGITUDE) / CELL_SIZE_DEGREES).rem_euclid(LONGITUDE_COUNT as f64);
    let lower = coordinate.floor() as usize;
    let upper = (lower + 1) % LONGITUDE_COUNT;
    (lower, upper, coordinate - lower as f64)
}

fn spatial_position(
    latitude: f64,
    longitude: f64,
    spatial_interpolation: SpatialInterpolation,
) -> (usize, usize, f64, usize, usize, f64) {
    let (lat0, lat1, lat_weight) = latitude_position(latitude);
    let (lon0, lon1, lon_weight) = longitude_position(longitude);

    match spatial_interpolation {
        SpatialInterpolation::Nearest => {
            let latitude = (lat0 as f64 + lat_weight).round_ties_even() as usize;
            let longitude =
                ((lon0 as f64 + lon_weight).round_ties_even() as usize) % LONGITUDE_COUNT;
            (latitude, latitude, 0.0, longitude, longitude, 0.0)
        }
        SpatialInterpolation::Bilinear => (lat0, lat1, lat_weight, lon0, lon1, lon_weight),
    }
}

#[cfg(test)]
mod tests {
    use std::sync::OnceLock;

    use chrono::{TimeZone, Utc};

    use crate::periodic_tables::tl::DEFAULT_DATASET_PATH;

    use super::{LinkeTurbidityGrid, SpatialInterpolation};

    fn grid() -> &'static LinkeTurbidityGrid {
        println!("Path: {DEFAULT_DATASET_PATH}");
        static GRID: OnceLock<LinkeTurbidityGrid> = OnceLock::new();
        GRID.get_or_init(|| {
            LinkeTurbidityGrid::load(DEFAULT_DATASET_PATH)
                .expect("the bundled Linke turbidity dataset should load from")
        })
    }

    fn assert_matches_pvlib(
        time: chrono::DateTime<Utc>,
        latitude: f64,
        longitude: f64,
        pvlib_value: f32,
    ) {
        let actual = grid()
            .interpolate_with_spatial_interpolation(
                time,
                latitude,
                longitude,
                SpatialInterpolation::Nearest,
            )
            .expect("pvlib reference locations should have valid turbidity values");

        assert!(
            (actual - pvlib_value).abs() <= 0.15,
            "actual {actual} differs too much from pvlib reference {pvlib_value}"
        );
    }

    #[test]
    fn nearest_lookup_is_close_to_pvlib_for_berlin_in_winter() {
        assert_matches_pvlib(
            Utc.with_ymd_and_hms(2024, 1, 15, 12, 0, 0).unwrap(),
            52.5,
            13.416_666_666_666_657,
            2.699_193_5,
        );
    }

    #[test]
    fn nearest_lookup_is_close_to_pvlib_for_phoenix_in_summer() {
        assert_matches_pvlib(
            Utc.with_ymd_and_hms(2024, 7, 15, 12, 0, 0).unwrap(),
            33.458_333_333_333_33,
            -112.041_666_666_666_67,
            3.448_360_7,
        );
    }
}
