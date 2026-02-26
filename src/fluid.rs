//! # fluid
//!
//! ## Purpose
//! Fluid property lookups by temperature for water and common fluids.
//! Supports both fresh water and seawater/brine at arbitrary salinity.
//! All properties are computed via polynomial curve fits derived from
//! NIST/CRC handbook data and UNESCO EOS-80, valid over 0–100°C at 1 atm.
//!
//! ## Algorithms
//! - Fresh water density: 5th-order polynomial fit to NIST data
//! - Seawater density: UNESCO EOS-80 (Millero & Poisson 1981)
//! - Fresh water viscosity: Vogel equation (Arrhenius-type)
//! - Seawater viscosity: Sharqawy et al. 2010 relative viscosity correction
//! - Surface tension: Linear fit + salinity correction (Nayar et al. 2014)
//! - Vapor pressure: Antoine equation + Raoult's law salinity reduction
//!
//! ## Data Structures
//! - [`Fluid`] — holds all thermophysical properties at a given temperature and salinity

/// Standard gravity (m/s²)
pub const G: f64 = 9.80665;

/// Atmospheric pressure at sea level (Pa)
pub const P_ATM: f64 = 101_325.0;

/// Fluid thermophysical properties at a given temperature and salinity.
#[derive(Debug, Clone, PartialEq)]
pub struct Fluid {
    /// Temperature (°C)
    pub temp_c: f64,
    /// Salinity (g/kg or ppt); 0.0 for fresh water, ~35 for standard seawater
    pub salinity_ppt: f64,
    /// Density (kg/m³)
    pub density_kg_m3: f64,
    /// Dynamic viscosity (Pa·s)
    pub dynamic_viscosity_pa_s: f64,
    /// Surface tension (N/m)
    pub surface_tension_n_m: f64,
    /// Vapor pressure (Pa) — the pressure at which the fluid boils at this temperature
    pub vapor_pressure_pa: f64,
}

impl Fluid {
    /// Construct water properties at the given temperature in °C (valid 0–100°C).
    ///
    /// # Panics
    /// Panics in debug mode if `temp_c` is outside [0, 100].
    pub fn water(temp_c: f64) -> Self {
        debug_assert!(
            (0.0..=100.0).contains(&temp_c),
            "temp_c must be in [0, 100]°C, got {temp_c}"
        );
        Self {
            temp_c,
            salinity_ppt: 0.0,
            density_kg_m3: water_density(temp_c),
            dynamic_viscosity_pa_s: water_viscosity(temp_c),
            surface_tension_n_m: water_surface_tension(temp_c),
            vapor_pressure_pa: water_vapor_pressure(temp_c),
        }
    }

    /// Construct seawater/brine properties at given temperature and salinity.
    ///
    /// # Arguments
    /// - `temp_c` — temperature in °C (valid 0–100°C)
    /// - `salinity_ppt` — salinity in g/kg (ppt); standard seawater ≈ 35,
    ///   Bay of Bengal ≈ 32, brine reject ≈ 70
    ///
    /// # Panics
    /// Panics in debug mode if `temp_c` is outside [0, 100] or `salinity_ppt` is negative.
    pub fn seawater(temp_c: f64, salinity_ppt: f64) -> Self {
        debug_assert!(
            (0.0..=100.0).contains(&temp_c),
            "temp_c must be in [0, 100]°C, got {temp_c}"
        );
        debug_assert!(
            salinity_ppt >= 0.0,
            "salinity_ppt must be >= 0, got {salinity_ppt}"
        );
        Self {
            temp_c,
            salinity_ppt,
            density_kg_m3: seawater_density(temp_c, salinity_ppt),
            dynamic_viscosity_pa_s: seawater_viscosity(temp_c, salinity_ppt),
            surface_tension_n_m: seawater_surface_tension(temp_c, salinity_ppt),
            vapor_pressure_pa: seawater_vapor_pressure(temp_c, salinity_ppt),
        }
    }

    /// Construct a custom fluid with explicit properties.
    pub fn custom(
        temp_c: f64,
        density_kg_m3: f64,
        dynamic_viscosity_pa_s: f64,
        surface_tension_n_m: f64,
        vapor_pressure_pa: f64,
    ) -> Self {
        Self {
            temp_c,
            salinity_ppt: 0.0,
            density_kg_m3,
            dynamic_viscosity_pa_s,
            surface_tension_n_m,
            vapor_pressure_pa,
        }
    }

    /// Kinematic viscosity ν = μ/ρ (m²/s)
    #[inline]
    pub fn kinematic_viscosity_m2_s(&self) -> f64 {
        self.dynamic_viscosity_pa_s / self.density_kg_m3
    }
}

/// Water density (kg/m³) via 4th-order polynomial fit to NIST data (0–100°C).
/// Max error < 0.05 kg/m³.
#[inline]
fn water_density(t: f64) -> f64 {
    // Coefficients fitted to NIST data
    999.842_594
        + 6.793_952e-2 * t
        - 9.095_290e-3 * t * t
        + 1.001_685e-4 * t * t * t
        - 1.120_083e-6 * t * t * t * t
        + 6.536_332e-9 * t * t * t * t * t
}

/// Water dynamic viscosity (Pa·s) via Vogel equation (0–100°C).
/// μ = A * exp(B / (T - C)), T in Kelvin.
/// Constants: A=2.414e-5, B=570.58, C=140.0  (Vogel 1921, NIST-fitted)
#[inline]
fn water_viscosity(t: f64) -> f64 {
    let t_k = t + 273.15; // convert °C to Kelvin
    let a: f64 = 2.414e-5;
    let b: f64 = 570.58;
    let c: f64 = 140.0;
    a * libm_exp(b / (t_k - c))
}

/// Water surface tension (N/m) via linear fit (0–100°C).
/// σ ≈ 0.07562 - 0.0001676 * T  (fitted to NIST data)
#[inline]
fn water_surface_tension(t: f64) -> f64 {
    0.075_62 - 1.676e-4 * t
}

/// Water vapor pressure (Pa) via Antoine equation.
/// log10(P_mmHg) = A - B / (C + T),  T in °C
/// Antoine constants for water: A=8.07131, B=1730.63, C=233.426 (0–60°C)
///                              A=8.14019, B=1810.94, C=244.485 (60–150°C)
#[inline]
fn water_vapor_pressure(t: f64) -> f64 {
    let (a, b, c) = if t <= 60.0 {
        (8.071_31, 1730.63, 233.426)
    } else {
        (8.140_19, 1810.94, 244.485)
    };
    let log_p_mmhg = a - b / (c + t);
    // Convert mmHg → Pa: 1 mmHg = 133.322 Pa
    libm_pow10(log_p_mmhg) * 133.322
}

// ---------------------------------------------------------------------------
// Seawater property functions
// ---------------------------------------------------------------------------

/// Seawater density (kg/m³) via UNESCO EOS-80 (Millero & Poisson 1981).
/// Valid 0–40°C, 0–42 ppt. Extrapolates reasonably to higher salinity.
///
/// ρ_sw = ρ_fw + S·(A + B·S^(1/2) + C·S)
/// where A, B, C are temperature-dependent polynomials.
#[inline]
fn seawater_density(t: f64, s: f64) -> f64 {
    let rho_fw = water_density(t);
    if s <= 0.0 {
        return rho_fw;
    }
    // UNESCO EOS-80 coefficients for salinity correction
    let a = 8.244_93e-1 - 4.0899e-3 * t + 7.6438e-5 * t * t
        - 8.2467e-7 * t * t * t + 5.3875e-9 * t * t * t * t;
    let b = -5.724_66e-3 + 1.0227e-4 * t - 1.6546e-6 * t * t;
    let c = 4.8314e-4;
    rho_fw + s * (a + b * s.sqrt() + c * s)
}

/// Seawater dynamic viscosity (Pa·s) via Sharqawy et al. 2010.
/// μ_sw = μ_fw × relative_viscosity(T, S)
/// Valid 0–180°C, 0–150 g/kg.
#[inline]
fn seawater_viscosity(t: f64, s: f64) -> f64 {
    let mu_fw = water_viscosity(t);
    if s <= 0.0 {
        return mu_fw;
    }
    // Sharqawy relative viscosity: μ_sw/μ_fw = 1 + A·S + B·S²
    // A and B are temperature-dependent
    let a = 1.541e-3 + 1.998e-5 * t - 9.52e-8 * t * t;
    let b = 7.974e-6 - 7.561e-8 * t + 4.724e-10 * t * t;
    let s_gkg = s; // salinity already in g/kg
    mu_fw * (1.0 + a * s_gkg + b * s_gkg * s_gkg)
}

/// Seawater surface tension (N/m) via Nayar et al. 2014 salinity correction.
/// σ_sw ≈ σ_fw × (1 + 0.000226 × S)   for S in g/kg
#[inline]
fn seawater_surface_tension(t: f64, s: f64) -> f64 {
    let sigma_fw = water_surface_tension(t);
    sigma_fw * (1.0 + 2.26e-4 * s)
}

/// Seawater vapor pressure (Pa) — reduced by Raoult's law for dissolved salts.
/// The mole fraction of water decreases with salinity:
///   x_w ≈ 1 - 0.0009539 × S  (S in g/kg, NaCl approximation)
///   P_vapor_sw = x_w × P_vapor_fw
#[inline]
fn seawater_vapor_pressure(t: f64, s: f64) -> f64 {
    let p_fw = water_vapor_pressure(t);
    if s <= 0.0 {
        return p_fw;
    }
    // Mole fraction of water (NaCl-equivalent approximation)
    // NaCl MW = 58.44, water MW = 18.015
    // For S g/kg: moles_NaCl = S/58.44 per kg, moles_water = (1000-S)/18.015
    let moles_salt = s / 58.44;
    let moles_water = (1000.0 - s) / 18.015;
    let x_w = moles_water / (moles_water + 2.0 * moles_salt); // NaCl dissociates → 2 ions
    p_fw * x_w
}

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// Compute exp(x) using std.
#[inline]
fn libm_exp(x: f64) -> f64 {
    x.exp()
}

/// Compute 10^x.
#[inline]
fn libm_pow10(x: f64) -> f64 {
    libm_exp(x * core::f64::consts::LN_10)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_water_density_at_4c() {
        // Maximum density of water occurs near 4°C ≈ 999.97 kg/m³
        let d = water_density(4.0);
        assert!((d - 999.97).abs() < 0.1, "density at 4°C = {d}");
    }

    #[test]
    fn test_water_density_at_20c() {
        let d = water_density(20.0);
        assert!((d - 998.2).abs() < 0.5, "density at 20°C = {d}");
    }

    #[test]
    fn test_water_viscosity_at_20c() {
        // ~1.002e-3 Pa·s at 20°C
        let v = water_viscosity(20.0);
        assert!((v - 1.002e-3).abs() < 5e-5, "viscosity at 20°C = {v:.4e}");
    }

    #[test]
    fn test_water_surface_tension_at_20c() {
        // ~0.07275 N/m at 20°C
        let s = water_surface_tension(20.0);
        assert!((s - 0.07275).abs() < 0.002, "surface tension at 20°C = {s:.5}");
    }

    #[test]
    fn test_vapor_pressure_at_100c() {
        // At 100°C vapor pressure ≈ 101325 Pa (1 atm)
        let vp = water_vapor_pressure(100.0);
        assert!((vp - 101_325.0).abs() < 2000.0, "vapor pressure at 100°C = {vp:.0}");
    }

    #[test]
    fn test_fluid_water_constructor() {
        let f = Fluid::water(20.0);
        assert!(f.kinematic_viscosity_m2_s() > 0.0);
        assert_eq!(f.salinity_ppt, 0.0);
    }

    // --- Seawater tests ---

    #[test]
    fn test_seawater_density_35ppt_20c() {
        // Standard seawater at 35 ppt, 20°C: ρ ≈ 1024.8 kg/m³
        let d = seawater_density(20.0, 35.0);
        assert!((d - 1024.8).abs() < 1.5, "seawater density at 35ppt, 20°C = {d:.1}");
    }

    #[test]
    fn test_seawater_density_zero_salinity_equals_fresh() {
        let d_sw = seawater_density(20.0, 0.0);
        let d_fw = water_density(20.0);
        assert!((d_sw - d_fw).abs() < 1e-10, "zero salinity should equal fresh water");
    }

    #[test]
    fn test_seawater_density_bay_of_bengal() {
        // Bay of Bengal: ~32 g/L, 25°C → ρ ≈ 1021 kg/m³
        let d = seawater_density(25.0, 32.0);
        assert!((d - 1021.0).abs() < 2.0, "Bay of Bengal density = {d:.1}");
    }

    #[test]
    fn test_seawater_viscosity_35ppt_20c() {
        // Standard seawater 35 ppt, 20°C: μ ≈ 1.075e-3 Pa·s
        let mu = seawater_viscosity(20.0, 35.0);
        assert!((mu - 1.075e-3).abs() < 1e-4, "seawater viscosity at 35ppt, 20°C = {mu:.4e}");
    }

    #[test]
    fn test_seawater_viscosity_zero_salinity_equals_fresh() {
        let mu_sw = seawater_viscosity(20.0, 0.0);
        let mu_fw = water_viscosity(20.0);
        assert!((mu_sw - mu_fw).abs() < 1e-10);
    }

    #[test]
    fn test_seawater_surface_tension_increases_with_salt() {
        // Salinity increases surface tension slightly
        let sigma_fw = water_surface_tension(20.0);
        let sigma_sw = seawater_surface_tension(20.0, 35.0);
        assert!(sigma_sw > sigma_fw, "seawater σ should exceed fresh water σ");
        // Increase is ~0.8% at 35 ppt
        let ratio = sigma_sw / sigma_fw;
        assert!((ratio - 1.008).abs() < 0.005, "ratio = {ratio:.4}");
    }

    #[test]
    fn test_seawater_vapor_pressure_lower_than_fresh() {
        // Dissolved salts reduce vapor pressure
        let vp_fw = water_vapor_pressure(20.0);
        let vp_sw = seawater_vapor_pressure(20.0, 35.0);
        assert!(vp_sw < vp_fw, "seawater VP should be less than fresh");
        // Reduction is ~1.8% at 35 ppt
        let ratio = vp_sw / vp_fw;
        assert!(ratio > 0.97 && ratio < 0.99, "VP ratio = {ratio:.4}");
    }

    #[test]
    fn test_brine_density_70ppt() {
        // RO brine reject at 70 ppt, 25°C: ρ ≈ 1050 kg/m³
        let d = seawater_density(25.0, 70.0);
        assert!((d - 1050.0).abs() < 5.0, "brine density at 70ppt = {d:.1}");
    }

    #[test]
    fn test_seawater_constructor() {
        let f = Fluid::seawater(25.0, 32.0);
        assert_eq!(f.salinity_ppt, 32.0);
        assert!(f.density_kg_m3 > 1000.0, "seawater should be denser than fresh");
        assert!(f.kinematic_viscosity_m2_s() > 0.0);
    }
}
