pub mod core; // for main.rs
mod periodic_tables;

#[cfg(test)]
mod test;
use ndarray::{Ix1, Ix2, Ix3};
use numpy::{
    IntoPyArray, PyArray3, PyReadonlyArray1, PyReadonlyArray3, PyReadonlyArrayDyn,
    datetime::{Datetime, units},
};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyModule};

use crate::core::irradiance::OpticalLossParameters;
use crate::core::{
    AoiInput, AtmosphericInput, SolarError, SolarPositionInput, SpatialInput, calculate_aoi,
    calculate_solar_position,
};

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
#[pyfunction(name = "calculate_solar_position")]
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
        num_threads,
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
    module.add_function(wrap_pyfunction!(calculate_solar_position_numpy, module)?)?;
    module.add_function(wrap_pyfunction!(calculate_aoi_numpy, module)?)?;
    Ok(())
}
