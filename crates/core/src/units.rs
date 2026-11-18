use crate::{Error, Result};

/// Trimmed CRS label; EPSG codes in `epsg:N`, `EPSG:N` or WKT ending in an EPSG
/// authority become `EPSG:N`. Anything else is kept verbatim.
pub fn normalize_crs(crs: &str) -> String {
    let s = crs.trim();
    if let Some(code) = s
        .get(..5)
        .filter(|p| p.eq_ignore_ascii_case("epsg:"))
        .map(|_| s[5..].trim())
        && !code.is_empty()
        && code.bytes().all(|b| b.is_ascii_digit())
    {
        return format!("EPSG:{code}");
    }
    if let Some(code) = wkt_epsg(s) {
        return format!("EPSG:{code}");
    }
    s.to_string()
}

fn wkt_epsg(s: &str) -> Option<&str> {
    let body = s.strip_suffix(']')?.trim_end();
    for key in ["AUTHORITY[", "ID["] {
        let Some(at) = body.rfind(key) else { continue };
        let inner = body[at + key.len()..].trim_end_matches(']');
        let mut fields = inner.split(',').map(|f| f.trim().trim_matches('"'));
        if fields.next()?.eq_ignore_ascii_case("EPSG") {
            let code = fields.next()?;
            if !code.is_empty()
                && code.bytes().all(|b| b.is_ascii_digit())
                && fields.next().is_none()
            {
                return Some(code);
            }
        }
    }
    None
}

/// Metres per `unit`.
pub fn length_factor(unit: &str) -> Result<f64> {
    Ok(match unit.trim().to_ascii_lowercase().as_str() {
        "m" | "metre" | "meter" | "metres" | "meters" => 1.0,
        "km" => 1000.0,
        "cm" => 0.01,
        "mm" => 0.001,
        "ft" | "foot" | "feet" => 0.3048,
        "us_ft" | "us_survey_foot" => 1200.0 / 3937.0,
        "in" => 0.0254,
        "yd" => 0.9144,
        _ => return Err(Error::invalid(format!("unknown length unit {unit:?}"))),
    })
}

fn same_unit(a: &str, b: &str) -> bool {
    match (length_factor(a), length_factor(b)) {
        (Ok(x), Ok(y)) => x == y,
        _ => a.trim() == b.trim(),
    }
}

/// The merged `(crs, length_unit)` of two labelled geometries. An undeclared
/// label matches anything; two declared labels must agree after normalization.
pub fn check_same(
    what: &str,
    a: (&Option<String>, &Option<String>),
    b: (&Option<String>, &Option<String>),
) -> Result<(Option<String>, Option<String>)> {
    let crs = match (a.0, b.0) {
        (Some(x), Some(y)) if normalize_crs(x) != normalize_crs(y) => {
            return Err(Error::CrsMismatch {
                left: format!("{what}: {x}"),
                right: y.clone(),
            });
        }
        (Some(x), _) | (None, Some(x)) => Some(x.clone()),
        (None, None) => None,
    };
    let unit = match (a.1, b.1) {
        (Some(x), Some(y)) if !same_unit(x, y) => {
            return Err(Error::UnitMismatch {
                left: format!("{what}: {x}"),
                right: y.clone(),
            });
        }
        (Some(x), _) | (None, Some(x)) => Some(x.clone()),
        (None, None) => None,
    };
    Ok((crs, unit))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crs_labels_normalize_epsg() {
        assert_eq!(normalize_crs(" epsg:4326 "), "EPSG:4326");
        assert_eq!(
            normalize_crs(r#"GEOGCS["WGS 84",AUTHORITY["EPSG","4326"]]"#),
            "EPSG:4326"
        );
        assert_eq!(
            normalize_crs(r#"GEOGCRS["WGS 84",ID["EPSG",4326]]"#),
            "EPSG:4326"
        );
        assert_eq!(normalize_crs("local grid"), "local grid");
        assert_eq!(normalize_crs("epsg:abc"), "epsg:abc");
    }

    #[test]
    fn length_factors_and_aliases() {
        assert_eq!(length_factor("m").unwrap(), 1.0);
        assert_eq!(length_factor("Feet").unwrap(), 0.3048);
        assert_eq!(
            length_factor("us_survey_foot").unwrap(),
            length_factor("us_ft").unwrap()
        );
        assert!(length_factor("furlong").is_err());
    }

    #[test]
    fn check_same_merges_or_refuses() {
        let m = Some("m".to_string());
        let ft = Some("ft".to_string());
        let e = Some("EPSG:31983".to_string());
        let e2 = Some("epsg:31983".to_string());
        let none = None;
        assert_eq!(
            check_same("a", (&none, &m), (&e, &none)).unwrap(),
            (e.clone(), m.clone())
        );
        assert_eq!(
            check_same("a", (&e, &m), (&e2, &Some("metre".into()))).unwrap(),
            (e.clone(), m.clone())
        );
        assert!(matches!(
            check_same("a", (&e, &m), (&e, &ft)),
            Err(Error::UnitMismatch { .. })
        ));
        assert!(matches!(
            check_same("a", (&e, &m), (&Some("EPSG:4326".into()), &m)),
            Err(Error::CrsMismatch { .. })
        ));
    }
}
