use crate::{Error, Result};

/// `value` when it is finite and `>= 0`, else `InvalidInput` naming `name`.
pub fn check(name: &str, value: f64) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(Error::invalid(format!(
            "{name} must be finite and >= 0, got {value}"
        )))
    }
}

/// `x` with `-0.0` mapped to `+0.0`, so equal values have equal bits.
pub fn canonical(x: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_accepts_finite_non_negative() {
        assert_eq!(check("t", 0.0).unwrap(), 0.0);
        assert_eq!(check("t", 2.5).unwrap(), 2.5);
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let e = check("weld", bad).unwrap_err().to_string();
            assert!(e.starts_with("weld must be"), "{e}");
        }
    }

    #[test]
    fn canonical_merges_signed_zero() {
        assert_eq!(canonical(-0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(canonical(1.5), 1.5);
    }
}
