use chrono::{DateTime, Datelike};
use ndarray::{Array3, ArrayView1, ArrayView2, ArrayView3};
use rayon::prelude::*;
use std::fmt::{Display, Formatter};

use crate::core::irradiance::etraterrestrial_radiation;
use crate::core::irradiance::poa_haydavies;
use crate::core::irradiance::{
    ClearSkyIrradiance, OpticalLossParameters, PoaIrradiance, ineichen_clearsky,
};
use crate::core::solar_position::{
    ObserverLatitudeGeometry, aberration_correction, apparent_sun_longitude,
    geocentric_sun_declination, geocentric_sun_right_ascension, obs_local_hour_angle,
    sidereal_time_greenwich, topocentric_azimuth_angle_e_from_n,
    topocentric_azimuth_angle_w_from_s_with_geometry, topocentric_elevation_angle_with_geometry,
    topocentric_sun_right_ascension_declination_and_local_hour_angle_with_geometry,
    topocentric_zenith_angle,
};
use crate::core::time::{calc_julian_day, delta_t_seconds, julian_century, julian_millennium};
use crate::core::with_thread_pool;
use crate::periodic_tables::earth::{
    B0_TABLE, B1_TABLE, CoordType, L0_TABLE, L1_TABLE, L2_TABLE, L3_TABLE, L4_TABLE, L5_TABLE,
    R0_TABLE, R1_TABLE, R2_TABLE, R3_TABLE, R4_TABLE, calculate_geocentric_coeff,
    calculate_heliocentric_coeff, sum_table,
};
use crate::periodic_tables::nutation::{calculate_dpsi_depsilon, calculate_epsilon};
use crate::periodic_tables::tl::{LinkeTurbidityGrid, SpatialInterpolation};

/// A value that is either constant over the spatial grid or specified per location.
pub(crate) enum SpatialInput<'a> {
    Scalar(f64),
    Grid(ArrayView2<'a, f64>),
}

/// An atmospheric value that changes per timestamp or per grid cell and timestamp.
pub(crate) enum AtmosphericInput<'a> {
    Time(ArrayView1<'a, f64>),
    Grid(ArrayView3<'a, f64>),
}

/// Borrowed inputs for array-first solar-position calculations.
///
/// The output shape is always `(time, lat, lon)`. Latitude and longitude are
/// independent one-dimensional axes. Elevation is scalar or `(lat, lon)`;
/// pressure and temperature are optional `(time)` or `(time, lat, lon)` arrays.
/// Time values are Unix timestamps in nanoseconds, matching NumPy
/// `datetime64[ns]` storage.
pub(crate) struct SolarPositionInput<'a> {
    pub(crate) latitude: ArrayView1<'a, f64>,
    pub(crate) longitude: ArrayView1<'a, f64>,
    pub(crate) time: ArrayView1<'a, i64>,
    pub(crate) elevation: SpatialInput<'a>,
    pub(crate) pressure: Option<AtmosphericInput<'a>>,
    pub(crate) temperature: Option<AtmosphericInput<'a>>,
}

/// Zenith and azimuth angles, in degrees, shaped `(time, lat, lon)`.
pub(crate) struct SolarPositionResult {
    pub(crate) zenith: Array3<f64>,
    pub(crate) azimuth: Array3<f64>,
}

/// Borrowed inputs for calculating the angle of incidence from solar-position output.
///
/// Panel geometry is scalar or `(lat, lon)` and is constant over time. Zenith
/// and azimuth must share the `(time, lat, lon)` output shape.
pub(crate) struct AoiInput<'a> {
    pub(crate) zenith: ArrayView3<'a, f64>,
    pub(crate) azimuth: ArrayView3<'a, f64>,
    pub(crate) panel_tilt: SpatialInput<'a>,
    pub(crate) panel_azimuth: SpatialInput<'a>,
    pub(crate) optical_loss_params: Option<OpticalLossParameters>,
}

/// Angle-of-incidence values, in degrees, shaped `(time, lat, lon)`.
pub(crate) struct AoiResult {
    pub(crate) aoi: Array3<f64>,
}

/// Borrowed inputs for Hay-Davies plane-of-array irradiance calculations.
///
/// Zenith, geometric AOI, DNI, GHI, and DHI use `(time, lat, lon)`. Panel
/// tilt and albedo are scalar or `(lat, lon)`. Extraterrestrial DNI is either
/// `(time)` or `(time, lat, lon)`.
pub(crate) struct PoaInput<'a> {
    pub(crate) zenith: ArrayView3<'a, f64>,
    pub(crate) aoi: ArrayView3<'a, f64>,
    pub(crate) panel_tilt: SpatialInput<'a>,
    pub(crate) dni: ArrayView3<'a, f64>,
    pub(crate) ghi: ArrayView3<'a, f64>,
    pub(crate) dhi: ArrayView3<'a, f64>,
    pub(crate) dni_extra: AtmosphericInput<'a>,
    pub(crate) albedo: AtmosphericInput<'a>,
}

/// Hay-Davies plane-of-array irradiance components, in W/m2.
pub(crate) struct PoaResult {
    pub(crate) global: Array3<f64>,
    pub(crate) direct: Array3<f64>,
    pub(crate) diffuse: Array3<f64>,
    pub(crate) sky_diffuse: Array3<f64>,
    pub(crate) ground_diffuse: Array3<f64>,
}

/// Borrowed inputs for Ineichen/Perez clear-sky irradiance calculations.
///
/// Zenith uses `(time, lat, lon)`. Time is Unix nanoseconds, latitude and
/// longitude are degree axes, elevation is scalar or `(lat, lon)`, and
/// pressure is optional in millibars.
pub(crate) struct ClearSkyInput<'a> {
    pub(crate) time: ArrayView1<'a, i64>,
    pub(crate) latitude: ArrayView1<'a, f64>,
    pub(crate) longitude: ArrayView1<'a, f64>,
    pub(crate) zenith: ArrayView3<'a, f64>,
    pub(crate) elevation: SpatialInput<'a>,
    pub(crate) pressure: Option<AtmosphericInput<'a>>,
}

/// Ineichen/Perez clear-sky irradiance components, in W/m2.
pub(crate) struct ClearSkyResult {
    pub(crate) ghi: Array3<f64>,
    pub(crate) dni: Array3<f64>,
    pub(crate) dhi: Array3<f64>,
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum SolarError {
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
        if let Some(pressure) = &self.pressure {
            validate_atmospheric_input("pressure", pressure, shape)?;
        }
        if let Some(temperature) = &self.temperature {
            validate_atmospheric_input("temperature", temperature, shape)?;
        }

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

impl<'a> PoaInput<'a> {
    pub(crate) fn output_shape(&self) -> Result<(usize, usize, usize), SolarError> {
        let shape = self.zenith.dim();
        for (name, values) in [
            ("aoi", self.aoi),
            ("dni", self.dni),
            ("ghi", self.ghi),
            ("dhi", self.dhi),
        ] {
            if values.dim() != shape {
                return Err(SolarError::InvalidShape(format!(
                    "{name} must have shape (time, lat, lon) = {shape:?}; got {:?}",
                    values.dim()
                )));
            }
        }

        validate_spatial_input("panel_tilt", &self.panel_tilt, (shape.1, shape.2))?;
        validate_atmospheric_input("albedo", &self.albedo, shape)?;
        validate_atmospheric_input("dni_extra", &self.dni_extra, shape)?;

        Ok(shape)
    }
}

impl<'a> ClearSkyInput<'a> {
    pub(crate) fn output_shape(&self) -> Result<(usize, usize, usize), SolarError> {
        let shape = self.zenith.dim();
        if self.time.len() != shape.0 {
            return Err(SolarError::InvalidShape(format!(
                "time must have shape (time) = ({},); got ({},)",
                shape.0,
                self.time.len()
            )));
        }
        if self.latitude.len() != shape.1 || self.longitude.len() != shape.2 {
            return Err(SolarError::InvalidShape(format!(
                "latitude and longitude must have lengths ({}, {}); got ({}, {})",
                shape.1,
                shape.2,
                self.latitude.len(),
                self.longitude.len()
            )));
        }

        validate_spatial_input("elevation", &self.elevation, (shape.1, shape.2))?;
        if let Some(pressure) = &self.pressure {
            validate_atmospheric_input("pressure", pressure, shape)?;
        }

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
pub(crate) fn calculate_solar_position(
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

    with_thread_pool(num_threads, || {
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
                    let pressure = map_atmospheric_input_to_value(
                        &input.pressure,
                        time_index,
                        latitude_index,
                        longitude_index,
                    );
                    let temperature = map_atmospheric_input_to_value(
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
    })?;

    Ok(SolarPositionResult {
        zenith: Array3::from_shape_vec(shape, zenith_values)
            .expect("output buffer length must match the requested shape"),
        azimuth: Array3::from_shape_vec(shape, azimuth_values)
            .expect("output buffer length must match the requested shape"),
    })
}

fn map_atmospheric_input_to_value(
    input: &Option<AtmosphericInput>,
    time_index: usize,
    latitude_index: usize,
    longitude_index: usize,
) -> Option<f64> {
    input
        .as_ref()
        .map(|values| atmospheric_value(values, time_index, latitude_index, longitude_index))
}

/// Calculates angle of incidence without recomputing solar position.
pub(crate) fn calculate_aoi(
    input: AoiInput<'_>,
    num_threads: usize,
) -> Result<AoiResult, SolarError> {
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

    with_thread_pool(num_threads, || {
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
    })?;

    Ok(AoiResult {
        aoi: Array3::from_shape_vec(shape, aoi_values)
            .expect("output buffer length must match the requested shape"),
    })
}

/// Calculates Ineichen/Perez clear-sky GHI, DNI, and DHI from precomputed zenith.
pub(crate) fn calculate_clearsky(
    input: ClearSkyInput<'_>,
    num_threads: usize,
    linke_turbidity: &LinkeTurbidityGrid,
) -> Result<ClearSkyResult, SolarError> {
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
    let mut ghi_values = vec![0.0; shape.0 * cells_per_time];
    let mut dni_values = vec![0.0; shape.0 * cells_per_time];
    let mut dhi_values = vec![0.0; shape.0 * cells_per_time];

    with_thread_pool(num_threads, || {
        ghi_values
            .par_iter_mut()
            .zip(dni_values.par_iter_mut())
            .zip(dhi_values.par_iter_mut())
            .enumerate()
            .for_each(|(output_index, ((ghi, dni), dhi))| {
                let time_index = output_index / cells_per_time;
                let cell_index = output_index % cells_per_time;
                let latitude_index = cell_index / shape.2;
                let longitude_index = cell_index % shape.2;
                let pressure = input.pressure.as_ref().map(|values| {
                    atmospheric_value(values, time_index, latitude_index, longitude_index)
                });
                let time = DateTime::from_timestamp(
                    input.time[time_index].div_euclid(1_000_000_000),
                    input.time[time_index].rem_euclid(1_000_000_000) as u32,
                )
                .expect("ClearSkyInput timestamps must be valid");
                let linke_turbidity = linke_turbidity
                    .interpolate_with_spatial_interpolation(
                        time,
                        input.latitude[latitude_index],
                        input.longitude[longitude_index],
                        SpatialInterpolation::Nearest,
                    )
                    .unwrap_or(0.0);
                let result: ClearSkyIrradiance = ineichen_clearsky(
                    input.zenith[[time_index, latitude_index, longitude_index]],
                    linke_turbidity as f64,
                    spatial_value(&input.elevation, latitude_index, longitude_index),
                    pressure,
                    etraterrestrial_radiation(time.ordinal() as i64),
                );
                *ghi = result.ghi;
                *dni = result.dni;
                *dhi = result.dhi;
            });
    })?;

    Ok(ClearSkyResult {
        ghi: Array3::from_shape_vec(shape, ghi_values)
            .expect("output buffer length must match the requested shape"),
        dni: Array3::from_shape_vec(shape, dni_values)
            .expect("output buffer length must match the requested shape"),
        dhi: Array3::from_shape_vec(shape, dhi_values)
            .expect("output buffer length must match the requested shape"),
    })
}

/// Calculates Hay-Davies plane-of-array irradiance from precomputed AOI and zenith.
pub(crate) fn calculate_poa(
    input: PoaInput<'_>,
    num_threads: usize,
) -> Result<PoaResult, SolarError> {
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
    let mut global_values = vec![0.0; shape.0 * cells_per_time];
    let mut direct_values = vec![0.0; shape.0 * cells_per_time];
    let mut diffuse_values = vec![0.0; shape.0 * cells_per_time];
    let mut sky_diffuse_values = vec![0.0; shape.0 * cells_per_time];
    let mut ground_diffuse_values = vec![0.0; shape.0 * cells_per_time];

    with_thread_pool(num_threads, || {
        global_values
            .par_iter_mut()
            .zip(direct_values.par_iter_mut())
            .zip(diffuse_values.par_iter_mut())
            .zip(sky_diffuse_values.par_iter_mut())
            .zip(ground_diffuse_values.par_iter_mut())
            .enumerate()
            .for_each(
                |(output_index, ((((global, direct), diffuse), sky_diffuse), ground_diffuse))| {
                    let time_index = output_index / cells_per_time;
                    let cell_index = output_index % cells_per_time;
                    let latitude_index = cell_index / shape.2;
                    let longitude_index = cell_index % shape.2;
                    let result = poa_haydavies(
                        spatial_value(&input.panel_tilt, latitude_index, longitude_index),
                        input.zenith[[time_index, latitude_index, longitude_index]],
                        input.ghi[[time_index, latitude_index, longitude_index]],
                        input.dni[[time_index, latitude_index, longitude_index]],
                        input.dhi[[time_index, latitude_index, longitude_index]],
                        atmospheric_value(
                            &input.dni_extra,
                            time_index,
                            latitude_index,
                            longitude_index,
                        ),
                        input.aoi[[time_index, latitude_index, longitude_index]],
                        atmospheric_value(
                            &input.albedo,
                            time_index,
                            latitude_index,
                            longitude_index,
                        ),
                    );

                    *global = result.global;
                    *direct = result.direct;
                    *diffuse = result.diffuse;
                    *sky_diffuse = result.sky_diffuse;
                    *ground_diffuse = result.ground_diffuse;
                },
            );
    })?;

    Ok(PoaResult {
        global: Array3::from_shape_vec(shape, global_values)
            .expect("output buffer length must match the requested shape"),
        direct: Array3::from_shape_vec(shape, direct_values)
            .expect("output buffer length must match the requested shape"),
        diffuse: Array3::from_shape_vec(shape, diffuse_values)
            .expect("output buffer length must match the requested shape"),
        sky_diffuse: Array3::from_shape_vec(shape, sky_diffuse_values)
            .expect("output buffer length must match the requested shape"),
        ground_diffuse: Array3::from_shape_vec(shape, ground_diffuse_values)
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
    pressure: Option<f64>,
    temperature: Option<f64>,
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
        pressure,
        temperature,
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
    use std::sync::OnceLock;

    use crate::periodic_tables::tl::DEFAULT_DATASET_PATH;

use super::{
        AoiInput, AtmosphericInput, ClearSkyInput, LinkeTurbidityGrid, PoaInput,
        SolarPositionInput, SpatialInput, calculate_aoi, calculate_clearsky, calculate_poa,
        calculate_solar_position,
    };
    use chrono::{TimeZone, Utc};
    use ndarray::{Array3, arr1, arr2};

    fn linke_turbidity_grid() -> &'static LinkeTurbidityGrid {
        static GRID: OnceLock<LinkeTurbidityGrid> = OnceLock::new();
        GRID.get_or_init(|| {
            LinkeTurbidityGrid::load(DEFAULT_DATASET_PATH)
            .expect("the bundled Linke turbidity dataset should load")
        })
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
                pressure: Some(AtmosphericInput::Grid(pressure.view())),
                temperature: Some(AtmosphericInput::Grid(temperature.view())),
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
    fn solar_position_accepts_missing_atmosphere() {
        let time = Utc.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();
        let latitude = arr1(&[40.0]);
        let longitude = arr1(&[-105.0]);
        let times = arr1(&[time.timestamp_nanos_opt().unwrap()]);

        let result = calculate_solar_position(
            SolarPositionInput {
                latitude: latitude.view(),
                longitude: longitude.view(),
                time: times.view(),
                elevation: SpatialInput::Scalar(1600.0),
                pressure: None,
                temperature: None,
            },
            1,
        )
        .unwrap();

        assert_eq!(result.zenith.dim(), (1, 1, 1));
        assert!(result.zenith[[0, 0, 0]].is_finite());
        assert!(result.azimuth[[0, 0, 0]].is_finite());
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

    #[test]
    fn grid_haydavies_poa_daytime_conditions() {
        let zenith = Array3::from_elem((1, 1, 1), 40.0);
        let azimuth = Array3::from_elem((1, 1, 1), 190.0);
        let geometric_aoi = calculate_aoi(
            AoiInput {
                zenith: zenith.view(),
                azimuth: azimuth.view(),
                panel_tilt: SpatialInput::Scalar(30.0),
                panel_azimuth: SpatialInput::Scalar(180.0),
                optical_loss_params: None,
            },
            1,
        )
        .unwrap();
        let dni = Array3::from_elem((1, 1, 1), 800.0);
        let ghi = Array3::from_elem((1, 1, 1), 600.0);
        let dhi = Array3::from_elem((1, 1, 1), 120.0);
        let dni_extra = arr1(&[1367.0]);
        let albedo = arr1(&[0.2]);

        let result = calculate_poa(
            PoaInput {
                zenith: zenith.view(),
                aoi: geometric_aoi.aoi.view(),
                panel_tilt: SpatialInput::Scalar(30.0),
                dni: dni.view(),
                ghi: ghi.view(),
                dhi: dhi.view(),
                dni_extra: AtmosphericInput::Time(dni_extra.view()),
                albedo: AtmosphericInput::Time(albedo.view()),
            },
            1,
        )
        .unwrap();

        assert!((result.global[[0, 0, 0]] - 928.251_756_634_908).abs() < 1e-10);
        assert!((result.direct[[0, 0, 0]] - 783.940_047_158_946).abs() < 1e-10);
        assert!((result.diffuse[[0, 0, 0]] - 144.311_709_475_962).abs() < 1e-10);
        assert!((result.sky_diffuse[[0, 0, 0]] - 136.273_233_703_029).abs() < 1e-10);
        assert!((result.ground_diffuse[[0, 0, 0]] - 8.038_475_772_934).abs() < 1e-10);
    }

    #[test]
    fn grid_ineichen_clearsky_with_millibar_pressure() {
        let time = Utc.with_ymd_and_hms(2024, 1, 15, 12, 0, 0).unwrap();
        let times = arr1(&[time.timestamp_nanos_opt().unwrap()]);
        let latitude = arr1(&[52.5]);
        let longitude = arr1(&[13.416_666_666_666_657]);
        let zenith = Array3::from_elem((1, 1, 1), 50.0);
        let pressure = arr1(&[1013.25]);

        let result = calculate_clearsky(
            ClearSkyInput {
                time: times.view(),
                latitude: latitude.view(),
                longitude: longitude.view(),
                zenith: zenith.view(),
                elevation: SpatialInput::Scalar(34.0),
                pressure: Some(AtmosphericInput::Time(pressure.view())),
            },
            1,
            linke_turbidity_grid(),
        )
        .unwrap();

        assert!((result.ghi[[0, 0, 0]] - 669.116_526_792_958).abs() < 1e-5);
        assert!((result.dni[[0, 0, 0]] - 921.089_422_814_756).abs() < 1e-5);
        assert!((result.dhi[[0, 0, 0]] - 77.051_658_394_307).abs() < 1e-5);
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
