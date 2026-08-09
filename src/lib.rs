pub mod core; // for main.rs
mod periodic_tables;

#[cfg(test)]
mod test;
use pyo3::prelude::*;

#[pymodule]
mod solars {
    use pyo3::prelude::*;

    /// Formats the sum of two numbers as string.
    #[pyfunction]
    pub fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
        Ok((a + b).to_string())
    }
}
