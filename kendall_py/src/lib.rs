use kendall;
use pyo3::prelude::*;

#[pyfunction]
fn tau_raw(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    Ok(kendall::tau_raw(&x, &y).unwrap())
}

#[pyfunction]
fn tau(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    Ok(kendall::tau(&x, &y).unwrap())
}

#[pyfunction]
fn zb(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    Ok(kendall::zb(&x, &y).unwrap())
}

#[pyfunction]
fn tau_avltree(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    Ok(kendall::tau_avltree(&x, &y).unwrap())
}

#[pyfunction]
fn zb_avltree(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    Ok(kendall::zb_avltree(&x, &y).unwrap())
}

#[pyfunction]
fn tau_binary(x: Vec<bool>, y: Vec<bool>) -> PyResult<f64> {
    Ok(kendall::tau_binary(&x, &y).unwrap())
}

#[pyfunction]
fn zb_binary(x: Vec<bool>, y: Vec<bool>) -> PyResult<f64> {
    Ok(kendall::zb_binary(&x, &y).unwrap())
}

#[pyfunction]
fn tau_zero(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    Ok(kendall::tau_zero(&x, &y).unwrap())
}

#[pyfunction]
fn zb_zero(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    Ok(kendall::zb_zero(&x, &y).unwrap())
}

#[pymodule]
fn kendall_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(tau_raw, m)?)?;
    m.add_function(wrap_pyfunction!(tau, m)?)?;
    m.add_function(wrap_pyfunction!(zb, m)?)?;
    m.add_function(wrap_pyfunction!(tau_avltree, m)?)?;
    m.add_function(wrap_pyfunction!(zb_avltree, m)?)?;
    m.add_function(wrap_pyfunction!(tau_binary, m)?)?;
    m.add_function(wrap_pyfunction!(zb_binary, m)?)?;
    m.add_function(wrap_pyfunction!(tau_zero, m)?)?;
    m.add_function(wrap_pyfunction!(zb_zero, m)?)?;
    Ok(())
}
