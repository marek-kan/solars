use chrono::DateTime;
use ndarray::{Array3, ArrayView1, ArrayView2, ArrayView3};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::fmt::{Display, Formatter};

use crate::core::irradiance::OpticalLossParameters;
use crate::core::solar_position::{
    ObserverLatitudeGeometry, aberration_correction, apparent_sun_longitude,
    geocentric_sun_declination, geocentric_sun_right_ascension, obs_local_hour_angle,
    sidereal_time_greenwich, topocentric_azimuth_angle_e_from_n,
    topocentric_azimuth_angle_w_from_s_with_geometry, topocentric_elevation_angle_with_geometry,
    topocentric_sun_right_ascension_declination_and_local_hour_angle_with_geometry,
    topocentric_zenith_angle,
};
use crate::core::time::{calc_julian_day, delta_t_seconds, julian_century, julian_millennium};
use crate::periodic_tables::earth::{
    B0_TABLE, B1_TABLE, CoordType, L0_TABLE, L1_TABLE, L2_TABLE, L3_TABLE, L4_TABLE, L5_TABLE,
    R0_TABLE, R1_TABLE, R2_TABLE, R3_TABLE, R4_TABLE, calculate_geocentric_coeff,
    calculate_heliocentric_coeff, sum_table,
};
use crate::periodic_tables::nutation::{calculate_dpsi_depsilon, calculate_epsilon};

/// A value that is either constant over the spatial grid or specified per location.
pub enum SpatialInput<'a> {
    Scalar(f64),
    Grid(ArrayView2<'a, f64>),
}

/// An atmospheric value that changes per timestamp or per grid cell and timestamp.
pub enum AtmosphericInput<'a> {
    Time(ArrayView1<'a, f64>),
    Grid(ArrayView3<'a, f64>),
}

/// Borrowed inputs for array-first solar-position calculations.
///
/// The output shape is always `(time, lat, lon)`. Latitude and longitude are
/// independent one-dimensional axes. Elevation is scalar or `(lat, lon)`;
/// pressure and temperature are `(time)` or `(time, lat, lon)`. Time values
/// are Unix timestamps in nanoseconds, matching NumPy `datetime64[ns]` storage.
pub struct SolarPositionInput<'a> {
    pub latitude: ArrayView1<'a, f64>,
    pub longitude: ArrayView1<'a, f64>,
    pub time: ArrayView1<'a, i64>,
    pub elevation: SpatialInput<'a>,
    pub pressure: AtmosphericInput<'a>,
    pub temperature: AtmosphericInput<'a>,
}

/// Zenith and azimuth angles, in degrees, shaped `(time, lat, lon)`.
pub struct SolarPositionResult {
    pub zenith: Array3<f64>,
    pub azimuth: Array3<f64>,
}

/// Borrowed inputs for calculating the angle of incidence from solar-position output.
///
/// Panel geometry is scalar or `(lat, lon)` and is constant over time. Zenith
/// and azimuth must share the `(time, lat, lon)` output shape.
pub struct AoiInput<'a> {
    pub zenith: ArrayView3<'a, f64>,
    pub azimuth: ArrayView3<'a, f64>,
    pub panel_tilt: SpatialInput<'a>,
    pub panel_azimuth: SpatialInput<'a>,
    pub optical_loss_params: Option<OpticalLossParameters>,
}

/// Angle-of-incidence values, in degrees, shaped `(time, lat, lon)`.
pub struct AoiResult {
    pub aoi: Array3<f64>,
}

#[derive(Debug, PartialEq, Eq)]
pub enum SolarError {
    InvalidShape(String),
    InvalidTimestamp(i64),
    InvalidThreadCount,
    ThreadPool(String),
}

impl Display for SolarError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidShape(message) => write!(formatter, "invalid array shape: {message}"),
            Self::InvalidTimestamp(value) => {
                write!(formatter, "invalid Unix timestamp in nanoseconds: {value}")
            }
            Self::InvalidThreadCount => write!(formatter, "thread count must be greater than zero"),
            Self::ThreadPool(message) => {
                write!(formatter, "failed to build thread pool: {message}")
            }
        }
    }
}

impl std::error::Error for SolarError {}

impl<'a> SolarPositionInput<'a> {
    pub(crate) fn output_shape(&self) -> Result<(usize, usize, usize), SolarError> {
        let shape = (self.time.len(), self.latitude.len(), self.longitude.len());

        validate_spatial_input("elevation", &self.elevation, (shape.1, shape.2))?;
        validate_atmospheric_input("pressure", &self.pressure, shape)?;
        validate_atmospheric_input("temperature", &self.temperature, shape)?;

        Ok(shape)
    }
}

impl<'a> AoiInput<'a> {
    pub(crate) fn output_shape(&self) -> Result<(usize, usize, usize), SolarError> {
        let shape = self.zenith.dim();
        if self.azimuth.dim() != shape {
            return Err(SolarError::InvalidShape(format!(
                "azimuth must have shape (time, lat, lon) = {shape:?}; got {:?}",
                self.azimuth.dim()
            )));
        }

        validate_spatial_input("panel_tilt", &self.panel_tilt, (shape.1, shape.2))?;
        validate_spatial_input("panel_azimuth", &self.panel_azimuth, (shape.1, shape.2))?;

        Ok(shape)
    }
}

struct TimeGeometry {
    radius: f64,
    right_ascension: f64,
    declination: f64,
    sidereal_time: f64,
}

/// Calculates zenith and azimuth for all requested times and grid locations.
///
/// Time-dependent solar geometry is evaluated once per timestamp and reused for
/// every latitude/longitude pair in that time slice.
pub fn calculate_solar_position(
    input: SolarPositionInput<'_>,
    num_threads: usize,
) -> Result<SolarPositionResult, SolarError> {
    if num_threads == 0 {
        return Err(SolarError::InvalidThreadCount);
    }

    let shape = input.output_shape()?;
    if shape.1 == 0 || shape.2 == 0 {
        return Err(SolarError::InvalidShape(
            "latitude and longitude must each contain at least one value".to_owned(),
        ));
    }

    let time_geometry = input
        .time
        .iter()
        .map(calculate_time_geometry)
        .collect::<Result<Vec<_>, _>>()?;
    let latitude_geometry = input
        .latitude
        .iter()
        .copied()
        .map(ObserverLatitudeGeometry::from_degrees)
        .collect::<Vec<_>>();
    let longitude_radians = input
        .longitude
        .iter()
        .copied()
        .map(f64::to_radians)
        .collect::<Vec<_>>();

    let cells_per_time = shape.1 * shape.2;
    let mut zenith_values = vec![0.0; shape.0 * cells_per_time];
    let mut azimuth_values = vec![0.0; shape.0 * cells_per_time];

    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|error| SolarError::ThreadPool(error.to_string()))?;

    thread_pool.install(|| {
        zenith_values
            .par_chunks_mut(2048)
            .zip(azimuth_values.par_chunks_mut(2048))
            .enumerate()
            .for_each(|(chunk_index, (zenith_chunk, azimuth_chunk))| {
                let output_start = chunk_index * 2048;

                for (offset, (zenith, azimuth)) in zenith_chunk
                    .iter_mut()
                    .zip(azimuth_chunk.iter_mut())
                    .enumerate()
                {
                    let output_index = output_start + offset;
                    let time_index = output_index / cells_per_time;
                    let cell_index = output_index % cells_per_time;
                    let latitude_index = cell_index / shape.2;
                    let longitude_index = cell_index % shape.2;
                    let elevation =
                        spatial_value(&input.elevation, latitude_index, longitude_index);
                    let pressure = atmospheric_value(
                        &input.pressure,
                        time_index,
                        latitude_index,
                        longitude_index,
                    );
                    let temperature = atmospheric_value(
                        &input.temperature,
                        time_index,
                        latitude_index,
                        longitude_index,
                    );
                    (*zenith, *azimuth) = calculate_cell(
                        &latitude_geometry[latitude_index],
                        longitude_radians[longitude_index],
                        elevation,
                        pressure,
                        temperature,
                        &time_geometry[time_index],
                    );
                }
            });
    });

    Ok(SolarPositionResult {
        zenith: Array3::from_shape_vec(shape, zenith_values)
            .expect("output buffer length must match the requested shape"),
        azimuth: Array3::from_shape_vec(shape, azimuth_values)
            .expect("output buffer length must match the requested shape"),
    })
}

/// Calculates angle of incidence without recomputing solar position.
pub fn calculate_aoi(input: AoiInput<'_>, num_threads: usize) -> Result<AoiResult, SolarError> {
    if num_threads == 0 {
        return Err(SolarError::InvalidThreadCount);
    }

    let shape = input.output_shape()?;
    if shape.1 == 0 || shape.2 == 0 {
        return Err(SolarError::InvalidShape(
            "latitude and longitude dimensions must each be greater than zero".to_owned(),
        ));
    }

    let cells_per_time = shape.1 * shape.2;
    let mut aoi_values = vec![0.0; shape.0 * cells_per_time];
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|error| SolarError::ThreadPool(error.to_string()))?;

    thread_pool.install(|| {
        aoi_values
            .par_chunks_mut(2048)
            .enumerate()
            .for_each(|(chunk_index, aoi_chunk)| {
                let output_start = chunk_index * 2048;

                for (offset, aoi) in aoi_chunk.iter_mut().enumerate() {
                    let output_index = output_start + offset;
                    let time_index = output_index / cells_per_time;
                    let cell_index = output_index % cells_per_time;
                    let latitude_index = cell_index / shape.2;
                    let longitude_index = cell_index % shape.2;
                    *aoi = crate::core::irradiance::aoi(
                        input.zenith[[time_index, latitude_index, longitude_index]],
                        input.azimuth[[time_index, latitude_index, longitude_index]],
                        spatial_value(&input.panel_tilt, latitude_index, longitude_index),
                        spatial_value(&input.panel_azimuth, latitude_index, longitude_index),
                        input.optical_loss_params,
                    );
                }
            });
    });

    Ok(AoiResult {
        aoi: Array3::from_shape_vec(shape, aoi_values)
            .expect("output buffer length must match the requested shape"),
    })
}

fn spatial_value(input: &SpatialInput<'_>, latitude_index: usize, longitude_index: usize) -> f64 {
    match input {
        SpatialInput::Scalar(value) => *value,
        SpatialInput::Grid(values) => values[[latitude_index, longitude_index]],
    }
}

fn atmospheric_value(
    input: &AtmosphericInput<'_>,
    time_index: usize,
    latitude_index: usize,
    longitude_index: usize,
) -> f64 {
    match input {
        AtmosphericInput::Time(values) => values[time_index],
        AtmosphericInput::Grid(values) => values[[time_index, latitude_index, longitude_index]],
    }
}

fn calculate_time_geometry(timestamp_nanoseconds: &i64) -> Result<TimeGeometry, SolarError> {
    let seconds = timestamp_nanoseconds.div_euclid(1_000_000_000);
    let nanoseconds = timestamp_nanoseconds.rem_euclid(1_000_000_000) as u32;
    let time = DateTime::from_timestamp(seconds, nanoseconds)
        .ok_or(SolarError::InvalidTimestamp(*timestamp_nanoseconds))?;
    let julian_day = calc_julian_day(&time);
    let julian_century_ut = julian_century(julian_day);
    let julian_day_tt = julian_day + delta_t_seconds(&time) / 86400.0;
    let julian_century_tt = julian_century(julian_day_tt);
    let julian_millennium_tt = julian_millennium(julian_century_tt);

    let l0 = sum_table(&L0_TABLE, &julian_millennium_tt);
    let l1 = sum_table(&L1_TABLE, &julian_millennium_tt);
    let l2 = sum_table(&L2_TABLE, &julian_millennium_tt);
    let l3 = sum_table(&L3_TABLE, &julian_millennium_tt);
    let l4 = sum_table(&L4_TABLE, &julian_millennium_tt);
    let l5 = sum_table(&L5_TABLE, &julian_millennium_tt);
    let b0 = sum_table(&B0_TABLE, &julian_millennium_tt);
    let b1 = sum_table(&B1_TABLE, &julian_millennium_tt);
    let r0 = sum_table(&R0_TABLE, &julian_millennium_tt);
    let r1 = sum_table(&R1_TABLE, &julian_millennium_tt);
    let r2 = sum_table(&R2_TABLE, &julian_millennium_tt);
    let r3 = sum_table(&R3_TABLE, &julian_millennium_tt);
    let r4 = sum_table(&R4_TABLE, &julian_millennium_tt);

    let heliocentric_longitude = calculate_heliocentric_coeff(
        julian_millennium_tt,
        l0,
        l1,
        l2,
        l3,
        l4,
        l5,
        CoordType::Longitude,
    );
    let heliocentric_latitude = calculate_heliocentric_coeff(
        julian_millennium_tt,
        b0,
        b1,
        0.0,
        0.0,
        0.0,
        0.0,
        CoordType::Latitude,
    );
    let radius = calculate_heliocentric_coeff(
        julian_millennium_tt,
        r0,
        r1,
        r2,
        r3,
        r4,
        0.0,
        CoordType::Radius,
    );
    let geocentric_longitude =
        calculate_geocentric_coeff(heliocentric_longitude, CoordType::Longitude);
    let geocentric_latitude =
        calculate_geocentric_coeff(heliocentric_latitude, CoordType::Latitude);
    let (delta_psi, delta_epsilon) = calculate_dpsi_depsilon(julian_century_tt);
    let epsilon = calculate_epsilon(julian_millennium_tt, delta_epsilon);
    let apparent_longitude = apparent_sun_longitude(
        geocentric_longitude,
        delta_psi,
        aberration_correction(radius),
    );
    let right_ascension =
        geocentric_sun_right_ascension(apparent_longitude, epsilon, geocentric_latitude);

    Ok(TimeGeometry {
        radius,
        right_ascension,
        declination: geocentric_sun_declination(apparent_longitude, epsilon, geocentric_latitude),
        sidereal_time: sidereal_time_greenwich(julian_day, julian_century_ut, delta_psi, epsilon),
    })
}

fn calculate_cell(
    latitude: &ObserverLatitudeGeometry,
    longitude_radians: f64,
    elevation: f64,
    pressure: f64,
    temperature: f64,
    geometry: &TimeGeometry,
) -> (f64, f64) {
    let hour_angle = obs_local_hour_angle(
        longitude_radians,
        geometry.sidereal_time,
        geometry.right_ascension,
    );
    let (topocentric_declination, topocentric_hour_angle) =
        topocentric_sun_right_ascension_declination_and_local_hour_angle_with_geometry(
            latitude,
            elevation,
            geometry.radius,
            hour_angle,
            geometry.declination,
        );
    let elevation_angle = topocentric_elevation_angle_with_geometry(
        latitude,
        topocentric_declination,
        topocentric_hour_angle,
        Some(pressure),
        Some(temperature),
    );
    let zenith = topocentric_zenith_angle(elevation_angle.to_degrees());
    let azimuth_west_from_south = topocentric_azimuth_angle_w_from_s_with_geometry(
        latitude,
        topocentric_hour_angle,
        topocentric_declination,
    );
    let azimuth = topocentric_azimuth_angle_e_from_n(azimuth_west_from_south.to_degrees());

    (zenith, azimuth)
}

#[cfg(test)]
mod tests {
    use super::{
        AoiInput, AtmosphericInput, SolarPositionInput, SpatialInput, calculate_aoi,
        calculate_solar_position,
    };
    use crate::core::calculate_scalar_solar_position;
    use chrono::{TimeZone, Utc};
    use ndarray::{Array3, arr1, arr2};

    #[test]
    fn array_calculation_matches_scalar_calculation_for_one_cell() {
        let latitude = arr1(&[39.742476]);
        let longitude = arr1(&[-105.1786]);
        let time = Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap();
        let times = arr1(&[time.timestamp_nanos_opt().unwrap()]);
        let result = calculate_solar_position(
            SolarPositionInput {
                latitude: latitude.view(),
                longitude: longitude.view(),
                time: times.view(),
                elevation: SpatialInput::Scalar(1830.14),
                pressure: AtmosphericInput::Time(arr1(&[820.0]).view()),
                temperature: AtmosphericInput::Time(arr1(&[11.0]).view()),
            },
            1,
        )
        .unwrap();
        let scalar_result = calculate_scalar_solar_position(
            39.742476,
            -105.1786,
            time,
            Some(1830.14),
            Some(11.0),
            Some(820.0),
            None,
            None,
        );

        assert_eq!(result.zenith.dim(), (1, 1, 1));
        assert_eq!(result.azimuth.dim(), (1, 1, 1));
        assert!((result.zenith[[0, 0, 0]] - scalar_result["zenith"]).abs() < 1e-12);
        assert!((result.azimuth[[0, 0, 0]] - scalar_result["azimuth"]).abs() < 1e-12);
    }

    #[test]
    fn aoi_uses_precomputed_array_solar_position() {
        let latitude = arr1(&[39.742476]);
        let longitude = arr1(&[-105.1786]);
        let time = Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap();
        let times = arr1(&[time.timestamp_nanos_opt().unwrap()]);
        let solar_position = calculate_solar_position(
            SolarPositionInput {
                latitude: latitude.view(),
                longitude: longitude.view(),
                time: times.view(),
                elevation: SpatialInput::Scalar(1830.14),
                pressure: AtmosphericInput::Time(arr1(&[820.0]).view()),
                temperature: AtmosphericInput::Time(arr1(&[11.0]).view()),
            },
            1,
        )
        .unwrap();
        let aoi = calculate_aoi(
            AoiInput {
                zenith: solar_position.zenith.view(),
                azimuth: solar_position.azimuth.view(),
                panel_tilt: SpatialInput::Scalar(30.0),
                panel_azimuth: SpatialInput::Scalar(180.0),
                optical_loss_params: None,
            },
            1,
        )
        .unwrap();
        let scalar_result = calculate_scalar_solar_position(
            39.742476,
            -105.1786,
            time,
            Some(1830.14),
            Some(11.0),
            Some(820.0),
            Some(30.0),
            Some(180.0),
        );

        assert_eq!(aoi.aoi.dim(), (1, 1, 1));
        assert!((aoi.aoi[[0, 0, 0]] - scalar_result["aoi"]).abs() < 1e-12);
    }

    #[test]
    fn solar_position_accepts_spatial_elevation_and_three_dimensional_atmosphere() {
        let time = Utc.with_ymd_and_hms(2003, 10, 17, 19, 30, 30).unwrap();
        let times = arr1(&[
            time.timestamp_nanos_opt().unwrap(),
            (time + chrono::Duration::hours(1))
                .timestamp_nanos_opt()
                .unwrap(),
        ]);
        let latitude = arr1(&[39.742476, 40.0]);
        let longitude = arr1(&[-105.1786, -105.0]);
        let elevation = arr2(&[[1830.14, 1800.0], [1750.0, 1700.0]]);
        let pressure = Array3::from_shape_vec(
            (2, 2, 2),
            vec![820.0, 815.0, 810.0, 805.0, 800.0, 795.0, 790.0, 785.0],
        )
        .unwrap();
        let temperature =
            Array3::from_shape_vec((2, 2, 2), vec![11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0])
                .unwrap();

        let result = calculate_solar_position(
            SolarPositionInput {
                latitude: latitude.view(),
                longitude: longitude.view(),
                time: times.view(),
                elevation: SpatialInput::Grid(elevation.view()),
                pressure: AtmosphericInput::Grid(pressure.view()),
                temperature: AtmosphericInput::Grid(temperature.view()),
            },
            2,
        )
        .unwrap();

        assert_eq!(result.zenith.dim(), (2, 2, 2));
        assert_eq!(result.azimuth.dim(), (2, 2, 2));
        assert!(result.zenith.iter().all(|value| value.is_finite()));
        assert!(
            result
                .azimuth
                .iter()
                .all(|value| (0.0..360.0).contains(value))
        );
    }

    #[test]
    fn aoi_accepts_spatial_panel_geometry() {
        let zenith = Array3::zeros((1, 2, 2));
        let azimuth = Array3::zeros((1, 2, 2));
        let panel_tilt = arr2(&[[0.0, 15.0], [30.0, 45.0]]);
        let panel_azimuth = arr2(&[[0.0, 90.0], [180.0, 270.0]]);

        let result = calculate_aoi(
            AoiInput {
                zenith: zenith.view(),
                azimuth: azimuth.view(),
                panel_tilt: SpatialInput::Grid(panel_tilt.view()),
                panel_azimuth: SpatialInput::Grid(panel_azimuth.view()),
                optical_loss_params: None,
            },
            2,
        )
        .unwrap();

        assert_eq!(result.aoi.dim(), (1, 2, 2));
        assert_eq!(result.aoi[[0, 0, 0]], 0.0);
        assert!((result.aoi[[0, 1, 1]] - 45.0).abs() < 1e-12);
    }
}

pub(crate) fn validate_spatial_input(
    name: &str,
    input: &SpatialInput<'_>,
    expected_shape: (usize, usize),
) -> Result<(), SolarError> {
    if let SpatialInput::Grid(values) = input
        && values.dim() != expected_shape
    {
        return Err(SolarError::InvalidShape(format!(
            "{name} must be scalar or have shape (lat, lon) = {expected_shape:?}; got {:?}",
            values.dim()
        )));
    }

    Ok(())
}

pub(crate) fn validate_atmospheric_input(
    name: &str,
    input: &AtmosphericInput<'_>,
    expected_shape: (usize, usize, usize),
) -> Result<(), SolarError> {
    match input {
        AtmosphericInput::Time(values) if values.len() != expected_shape.0 => {
            Err(SolarError::InvalidShape(format!(
                "{name} must have shape (time) = ({},); got ({},)",
                expected_shape.0,
                values.len()
            )))
        }
        AtmosphericInput::Grid(values) if values.dim() != expected_shape => {
            Err(SolarError::InvalidShape(format!(
                "{name} must have shape (time, lat, lon) = {expected_shape:?}; got {:?}",
                values.dim()
            )))
        }
        _ => Ok(()),
    }
}
