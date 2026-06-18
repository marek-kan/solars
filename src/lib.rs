use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
mod solars {
    use pyo3::prelude::*;

    /// Formats the sum of two numbers as string.
    #[pyfunction]
    pub fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
        Ok((a + b).to_string())
    }
}

#[cfg(test)]
mod test {
    use crate::solars::*;

    #[test]
    fn as_string() {
        let res = sum_as_string(10, 20).unwrap();

        assert_eq!("30", res);
    }
}