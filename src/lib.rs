pub(crate) mod core;
pub(crate) mod periodic_tables;

use chrono::{DateTime, Datelike};
use ndarray::{Array1, Ix1, Ix2, Ix3};
use numpy::{
    IntoPyArray, PyArray3, PyReadonlyArray1, PyReadonlyArray3, PyReadonlyArrayDyn,
    datetime::{Datetime, units},
};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyModule};

use crate::core::grid::{
    AoiInput, AtmosphericInput, ClearSkyInput, PoaInput, SolarError, SolarPositionInput,
    SpatialInput, calculate_aoi, calculate_clearsky, calculate_poa, calculate_solar_position,
};
use crate::core::irradiance::{OpticalLossParameters, etraterrestrial_radiation};

#[pyclass(name = "SolarPositionResult")]
struct PySolarPositionResult {
    zenith: Py<PyArray3<f64>>,
    azimuth: Py<PyArray3<f64>>,
}

#[pymethods]
impl PySolarPositionResult {
    #[getter]
    fn zenith(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.zenith.clone_ref(py)
    }

    #[getter]
    fn azimuth(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.azimuth.clone_ref(py)
    }
}

#[pyclass(name = "AoiResult")]
struct PyAoiResult {
    aoi: Py<PyArray3<f64>>,
}

#[pyclass(name = "ClearSkyResult")]
struct PyClearSkyResult {
    ghi: Py<PyArray3<f64>>,
    dni: Py<PyArray3<f64>>,
    dhi: Py<PyArray3<f64>>,
}

#[pymethods]
impl PyClearSkyResult {
    #[getter]
    fn ghi(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.ghi.clone_ref(py)
    }

    #[getter]
    fn dni(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.dni.clone_ref(py)
    }

    #[getter]
    fn dhi(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.dhi.clone_ref(py)
    }
}

#[pyclass(name = "PoaResult")]
struct PyPoaResult {
    global: Py<PyArray3<f64>>,
    direct: Py<PyArray3<f64>>,
    diffuse: Py<PyArray3<f64>>,
    sky_diffuse: Py<PyArray3<f64>>,
    ground_diffuse: Py<PyArray3<f64>>,
}

#[pymethods]
impl PyPoaResult {
    #[getter]
    fn global(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.global.clone_ref(py)
    }

    #[getter]
    fn direct(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.direct.clone_ref(py)
    }

    #[getter]
    fn diffuse(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.diffuse.clone_ref(py)
    }

    #[getter]
    fn sky_diffuse(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.sky_diffuse.clone_ref(py)
    }

    #[getter]
    fn ground_diffuse(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.ground_diffuse.clone_ref(py)
    }
}

#[pymethods]
impl PyAoiResult {
    #[getter]
    fn aoi(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.aoi.clone_ref(py)
    }
}

/// Compute solar zenith and azimuth arrays shaped `(time, lat, lon)`.
///
/// `latitude` and `longitude` are one-dimensional degree axes; `time` is a
/// one-dimensional `numpy.datetime64[ns]` array. Elevation is scalar or
/// `(lat, lon)`. Pressure and temperature are `(time)` or `(time, lat, lon)`.
#[pyfunction(
    name = "calculate_solar_position",
    signature = (
        latitude,
        longitude,
        time,
        elevation,
        pressure,
        temperature,
        num_threads = 1
    )
)]
fn calculate_solar_position_numpy<'py>(
    py: Python<'py>,
    latitude: PyReadonlyArray1<'py, f64>,
    longitude: PyReadonlyArray1<'py, f64>,
    time: PyReadonlyArray1<'py, Datetime<units::Nanoseconds>>,
    elevation: &Bound<'py, PyAny>,
    pressure: PyReadonlyArrayDyn<'py, f64>,
    temperature: PyReadonlyArrayDyn<'py, f64>,
    num_threads: usize,
) -> PyResult<PySolarPositionResult> {
    let elevation_scalar = elevation.extract::<f64>().ok();
    let elevation_array = if elevation_scalar.is_none() {
        Some(elevation.extract::<PyReadonlyArrayDyn<f64>>()?)
    } else {
        None
    };
    let elevation_input = spatial_input("elevation", elevation_scalar, &elevation_array)?;
    let pressure_input = atmospheric_input("pressure", &pressure)?;
    let temperature_input = atmospheric_input("temperature", &temperature)?;
    let time_values = time.as_array().mapv(i64::from);

    let result = calculate_solar_position(
        SolarPositionInput {
            latitude: latitude.as_array(),
            longitude: longitude.as_array(),
            time: time_values.view(),
            elevation: elevation_input,
            pressure: pressure_input,
            temperature: temperature_input,
        },
        num_threads,
    )
    .map_err(to_python_error)?;

    Ok(PySolarPositionResult {
        zenith: result.zenith.into_pyarray(py).unbind(),
        azimuth: result.azimuth.into_pyarray(py).unbind(),
    })
}

/// Compute angle of incidence from precomputed solar-position arrays.
///
/// Zenith and azimuth use `(time, lat, lon)`. Panel tilt and panel azimuth are
/// scalar or `(lat, lon)` degree arrays and are constant over time. Set
/// `apply_optical_loss` to return physical IAM instead of geometric AOI.
#[pyfunction(
    name = "calculate_aoi",
    signature = (
        zenith,
        azimuth,
        panel_tilt,
        panel_azimuth,
        num_threads = 1,
        apply_optical_loss = false,
        refractive_index = 1.526,
        extinction_coefficient = 4.0,
        thickness = 0.002,
    )
)]
fn calculate_aoi_numpy<'py>(
    py: Python<'py>,
    zenith: PyReadonlyArray3<'py, f64>,
    azimuth: PyReadonlyArray3<'py, f64>,
    panel_tilt: &Bound<'py, PyAny>,
    panel_azimuth: &Bound<'py, PyAny>,
    num_threads: usize,
    apply_optical_loss: bool,
    refractive_index: f64,
    extinction_coefficient: f64,
    thickness: f64,
) -> PyResult<PyAoiResult> {
    let panel_tilt_scalar = panel_tilt.extract::<f64>().ok();
    let panel_tilt_array = if panel_tilt_scalar.is_none() {
        Some(panel_tilt.extract::<PyReadonlyArrayDyn<f64>>()?)
    } else {
        None
    };
    let panel_azimuth_scalar = panel_azimuth.extract::<f64>().ok();
    let panel_azimuth_array = if panel_azimuth_scalar.is_none() {
        Some(panel_azimuth.extract::<PyReadonlyArrayDyn<f64>>()?)
    } else {
        None
    };

    let result = calculate_aoi(
        AoiInput {
            zenith: zenith.as_array(),
            azimuth: azimuth.as_array(),
            panel_tilt: spatial_input("panel_tilt", panel_tilt_scalar, &panel_tilt_array)?,
            panel_azimuth: spatial_input(
                "panel_azimuth",
                panel_azimuth_scalar,
                &panel_azimuth_array,
            )?,
            optical_loss_params: apply_optical_loss.then_some(OpticalLossParameters {
                refractive_index,
                extinction_coefficient,
                thickness,
            }),
        },
        num_threads,
    )
    .map_err(to_python_error)?;

    Ok(PyAoiResult {
        aoi: result.aoi.into_pyarray(py).unbind(),
    })
}

/// Compute Ineichen/Perez clear-sky irradiance arrays shaped `(time, lat, lon)`.
///
/// `time` is a one-dimensional `numpy.datetime64[ns]` array. Zenith uses
/// `(time, lat, lon)`. Elevation is scalar or `(lat, lon)`. Pressure is
/// optional, in millibars, and may be `(time)` or `(time, lat, lon)`.
#[pyfunction(
    name = "calculate_clearsky",
    signature = (latitude, longitude, time, zenith, elevation, pressure = None, num_threads = 1)
)]
fn calculate_clearsky_numpy<'py>(
    py: Python<'py>,
    latitude: PyReadonlyArray1<'py, f64>,
    longitude: PyReadonlyArray1<'py, f64>,
    time: PyReadonlyArray1<'py, Datetime<units::Nanoseconds>>,
    zenith: PyReadonlyArray3<'py, f64>,
    elevation: &Bound<'py, PyAny>,
    pressure: Option<PyReadonlyArrayDyn<'py, f64>>,
    num_threads: usize,
) -> PyResult<PyClearSkyResult> {
    let elevation_scalar = elevation.extract::<f64>().ok();
    let elevation_array = if elevation_scalar.is_none() {
        Some(elevation.extract::<PyReadonlyArrayDyn<f64>>()?)
    } else {
        None
    };
    let time_values = time.as_array().mapv(i64::from);
    let pressure_input = pressure
        .as_ref()
        .map(|values| atmospheric_input("pressure", values))
        .transpose()?;
    let result = calculate_clearsky(
        ClearSkyInput {
            time: time_values.view(),
            latitude: latitude.as_array(),
            longitude: longitude.as_array(),
            zenith: zenith.as_array(),
            elevation: spatial_input("elevation", elevation_scalar, &elevation_array)?,
            pressure: pressure_input,
        },
        num_threads,
    )
    .map_err(to_python_error)?;

    Ok(PyClearSkyResult {
        ghi: result.ghi.into_pyarray(py).unbind(),
        dni: result.dni.into_pyarray(py).unbind(),
        dhi: result.dhi.into_pyarray(py).unbind(),
    })
}

/// Compute Hay-Davies plane-of-array irradiance arrays shaped `(time, lat, lon)`.
///
/// Zenith, AOI, DNI, GHI, and DHI use `(time, lat, lon)`. Panel tilt is
/// scalar or `(lat, lon)` and albedo is scalar. `time` derives extraterrestrial DNI.
#[pyfunction(
    name = "calculate_poa",
    signature = (time, zenith, aoi, panel_tilt, dni, ghi, dhi, albedo, num_threads = 1)
)]
fn calculate_poa_numpy<'py>(
    py: Python<'py>,
    time: PyReadonlyArray1<'py, Datetime<units::Nanoseconds>>,
    zenith: PyReadonlyArray3<'py, f64>,
    aoi: PyReadonlyArray3<'py, f64>,
    panel_tilt: &Bound<'py, PyAny>,
    dni: PyReadonlyArray3<'py, f64>,
    ghi: PyReadonlyArray3<'py, f64>,
    dhi: PyReadonlyArray3<'py, f64>,
    albedo: PyReadonlyArrayDyn<'py, f64>,
    num_threads: usize,
) -> PyResult<PyPoaResult> {
    let panel_tilt_scalar = panel_tilt.extract::<f64>().ok();
    let panel_tilt_array = if panel_tilt_scalar.is_none() {
        Some(panel_tilt.extract::<PyReadonlyArrayDyn<f64>>()?)
    } else {
        None
    };

    let albedo_input = atmospheric_input("albedo", &albedo)?;

    let time_values = time.as_array().mapv(i64::from);
    let dni_extra = time_values
        .iter()
        .map(|timestamp| {
            let seconds = timestamp.div_euclid(1_000_000_000);
            let nanoseconds = timestamp.rem_euclid(1_000_000_000) as u32;
            let time = DateTime::from_timestamp(seconds, nanoseconds).ok_or_else(|| {
                PyValueError::new_err(format!("invalid Unix timestamp: {timestamp}"))
            })?;
            Ok(etraterrestrial_radiation(time.ordinal() as i64))
        })
        .collect::<PyResult<Vec<_>>>()?;

    let dni_extra = Array1::from(dni_extra);

    let result = calculate_poa(
        PoaInput {
            zenith: zenith.as_array(),
            aoi: aoi.as_array(),
            panel_tilt: spatial_input("panel_tilt", panel_tilt_scalar, &panel_tilt_array)?,
            dni: dni.as_array(),
            ghi: ghi.as_array(),
            dhi: dhi.as_array(),
            dni_extra: AtmosphericInput::Time(dni_extra.view()),
            albedo: albedo_input,
        },
        num_threads,
    )
    .map_err(to_python_error)?;

    Ok(PyPoaResult {
        global: result.global.into_pyarray(py).unbind(),
        direct: result.direct.into_pyarray(py).unbind(),
        diffuse: result.diffuse.into_pyarray(py).unbind(),
        sky_diffuse: result.sky_diffuse.into_pyarray(py).unbind(),
        ground_diffuse: result.ground_diffuse.into_pyarray(py).unbind(),
    })
}

fn spatial_input<'a, 'py>(
    name: &str,
    scalar: Option<f64>,
    array: &'a Option<PyReadonlyArrayDyn<'py, f64>>,
) -> PyResult<SpatialInput<'a>> {
    if let Some(value) = scalar {
        return Ok(SpatialInput::Scalar(value));
    }

    let values = array
        .as_ref()
        .expect("an array is required when no scalar value is provided")
        .as_array();
    let values = values.into_dimensionality::<Ix2>().map_err(|_| {
        PyValueError::new_err(format!("{name} must be a scalar or a 2D (lat, lon) array"))
    })?;

    Ok(SpatialInput::Grid(values))
}

fn atmospheric_input<'a, 'py>(
    name: &str,
    values: &'a PyReadonlyArrayDyn<'py, f64>,
) -> PyResult<AtmosphericInput<'a>> {
    let values = values.as_array();

    match values.ndim() {
        1 => Ok(AtmosphericInput::Time(
            values
                .into_dimensionality::<Ix1>()
                .expect("array dimensionality was checked"),
        )),
        3 => Ok(AtmosphericInput::Grid(
            values
                .into_dimensionality::<Ix3>()
                .expect("array dimensionality was checked"),
        )),
        _ => Err(PyValueError::new_err(format!(
            "{name} must be a 1D (time) or 3D (time, lat, lon) array"
        ))),
    }
}

fn to_python_error(error: SolarError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[pymodule]
fn solars(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySolarPositionResult>()?;
    module.add_class::<PyAoiResult>()?;
    module.add_class::<PyClearSkyResult>()?;
    module.add_class::<PyPoaResult>()?;
    module.add_function(wrap_pyfunction!(calculate_solar_position_numpy, module)?)?;
    module.add_function(wrap_pyfunction!(calculate_aoi_numpy, module)?)?;
    module.add_function(wrap_pyfunction!(calculate_clearsky_numpy, module)?)?;
    module.add_function(wrap_pyfunction!(calculate_poa_numpy, module)?)?;
    module.add(
        "__all__",
        (
            "SolarPositionResult",
            "AoiResult",
            "ClearSkyResult",
            "PoaResult",
            "calculate_solar_position",
            "calculate_aoi",
            "calculate_clearsky",
            "calculate_poa",
        ),
    )?;
    Ok(())
}

#[cfg(test)]
mod test;
