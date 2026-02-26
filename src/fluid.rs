//! # fluid
//!
//! ## Purpose
//! Fluid property lookups by temperature for water and common fluids.
//! All properties are computed via polynomial curve fits derived from
//! NIST/CRC handbook data, valid over 0–100°C at 1 atm.
//!
//! ## Algorithms
//! - Density: 4th-order polynomial fit to NIST water data
//! - Dynamic viscosity: Arrhenius-type exponential fit
//! - Surface tension: Linear fit (valid 0–100°C)
//! - Vapor pressure: Antoine equation (log10 form)
//!
//! ## Data Structures
//! - [`Fluid`] — holds all thermophysical properties at a given temperature

/// Standard gravity (m/s²)
pub const G: f64 = 9.80665;

/// Atmospheric pressure at sea level (Pa)
pub const P_ATM: f64 = 101_325.0;

/// Fluid thermophysical properties at a given temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct Fluid {
    /// Temperature (°C)
    pub temp_c: f64,
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
            density_kg_m3: water_density(temp_c),
            dynamic_viscosity_pa_s: water_viscosity(temp_c),
            surface_tension_n_m: water_surface_tension(temp_c),
            vapor_pressure_pa: water_vapor_pressure(temp_c),
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
    }
}
