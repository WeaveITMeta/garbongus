//! # capillary
//!
//! ## Purpose
//! Capillary rise height in circular tubes using the Jurin's Law equation.
//!
//! ## Algorithms
//! Jurin's Law: h = (2 * γ * cos(θ)) / (ρ * g * r)
//!
//! where:
//! - γ = surface tension (N/m)
//! - θ = contact angle (radians)
//! - ρ = fluid density (kg/m³)
//! - g = gravitational acceleration (m/s²)
//! - r = inner tube radius (m)
//!
//! ## Data Structures
//! - [`CapillaryRise`] — input parameters
//! - [`CapillaryResult`] — computed outputs

use crate::fluid::{Fluid, G};

/// Input parameters for capillary rise calculation.
#[derive(Debug, Clone)]
pub struct CapillaryRise {
    /// Fluid properties
    pub fluid: Fluid,
    /// Inner tube radius (m)
    pub radius_m: f64,
    /// Contact angle between fluid and tube wall (degrees)
    pub contact_angle_deg: f64,
}

/// Results of a capillary rise calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct CapillaryResult {
    /// Capillary rise height (m); negative means fluid is depressed
    pub height_m: f64,
    /// Capillary pressure driving the rise (Pa)
    pub capillary_pressure_pa: f64,
    /// Weight of fluid column balanced at equilibrium (N/m²)
    pub hydrostatic_pressure_pa: f64,
}

impl CapillaryRise {
    /// Create a new capillary rise calculation.
    ///
    /// # Arguments
    /// - `fluid` — fluid properties (use [`Fluid::water`])
    /// - `radius_m` — inner tube radius in meters
    /// - `contact_angle_deg` — contact angle in degrees (0° = fully wetting)
    pub fn new(fluid: Fluid, radius_m: f64, contact_angle_deg: f64) -> Self {
        Self {
            fluid,
            radius_m,
            contact_angle_deg,
        }
    }

    /// Compute the capillary rise height via Jurin's Law.
    pub fn calculate(&self) -> CapillaryResult {
        let cos_theta = cos_deg(self.contact_angle_deg);
        // Capillary pressure: ΔP = 2γ·cos(θ) / r
        let capillary_pressure_pa =
            2.0 * self.fluid.surface_tension_n_m * cos_theta / self.radius_m;
        // Equilibrium height: h = ΔP / (ρ·g)
        let height_m = capillary_pressure_pa / (self.fluid.density_kg_m3 * G);
        let hydrostatic_pressure_pa = self.fluid.density_kg_m3 * G * height_m;

        CapillaryResult {
            height_m,
            capillary_pressure_pa,
            hydrostatic_pressure_pa,
        }
    }
}

/// Compute cosine of angle given in degrees.
#[inline]
fn cos_deg(deg: f64) -> f64 {
    let rad = deg * core::f64::consts::PI / 180.0;
    rad.cos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid::Fluid;

    #[test]
    fn test_capillary_rise_1mm_water_20c() {
        // 1mm radius tube, water at 20°C, contact angle 0°
        // Expected: h = 2 * 0.07275 / (998.2 * 9.80665 * 0.001) ≈ 0.01486 m ≈ 14.86 mm
        let fluid = Fluid::water(20.0);
        let calc = CapillaryRise::new(fluid, 0.001, 0.0);
        let result = calc.calculate();
        assert!(
            (result.height_m - 0.01486).abs() < 0.001,
            "height = {:.5} m",
            result.height_m
        );
    }

    #[test]
    fn test_capillary_depression_mercury() {
        // Mercury has contact angle ~140° — it should be depressed
        let fluid = Fluid::custom(20.0, 13_534.0, 1.526e-3, 0.4865, 0.16);
        let calc = CapillaryRise::new(fluid, 0.001, 140.0);
        let result = calc.calculate();
        assert!(result.height_m < 0.0, "Mercury should be depressed: {}", result.height_m);
    }

    #[test]
    fn test_capillary_pressure_positive_wetting() {
        let fluid = Fluid::water(20.0);
        let calc = CapillaryRise::new(fluid, 0.001, 0.0);
        let result = calc.calculate();
        assert!(result.capillary_pressure_pa > 0.0);
    }
}
