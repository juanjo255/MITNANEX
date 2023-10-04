use pyo3::prelude::*;
use crate::kmer_cnt::get_kmer_profiles;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok(get_kmer_profiles())
}

/// A Python module implemented in Rust.
#[pymodule]
fn src(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    Ok(())
}