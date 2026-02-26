//! # pipe
//!
//! ## Purpose
//! Pipe flow calculations: Reynolds number, flow regime detection,
//! Darcy-Weisbach friction factor, and pressure loss.
//!
//! ## Algorithms
//!
//! ### Reynolds Number
//! ```text
//! Re = ρ·v·D / μ = v·D / ν
//! ```
//!
//! ### Flow Regimes
//! - Re < 2300 → Laminar
//! - 2300 ≤ Re < 4000 → Transitional
//! - Re ≥ 4000 → Turbulent
//!
//! ### Friction Factor
//! - **Laminar**: f = 64 / Re  (Hagen-Poiseuille exact)
//! - **Turbulent**: Colebrook-White equation solved via Swamee-Jain approximation:
//!   ```text
//!   f = 0.25 / [log10(ε/(3.7·D) + 5.74/Re^0.9)]²
//!   ```
//!   Or iteratively via Colebrook-White for higher accuracy.
//!
//! ### Darcy-Weisbach Pressure Loss
//! ```text
//! ΔP = f · (L/D) · (ρ·v²/2)
//! ```
//!
//! ## Data Structures
//! - [`DarcyWeisbach`] — input parameters
//! - [`PipeFlowResult`] — computed outputs including Re, f, ΔP

use crate::fluid::Fluid;

/// Flow regime classification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FlowRegime {
    /// Re < 2300: viscous, ordered flow
    Laminar,
    /// 2300 ≤ Re < 4000: unstable, unpredictable
    Transitional,
    /// Re ≥ 4000: chaotic, momentum-dominated flow
    Turbulent,
}

impl FlowRegime {
    /// Classify flow regime from Reynolds number.
    pub fn from_reynolds(re: f64) -> Self {
        if re < 2300.0 {
            FlowRegime::Laminar
        } else if re < 4000.0 {
            FlowRegime::Transitional
        } else {
            FlowRegime::Turbulent
        }
    }
}

/// Input parameters for Darcy-Weisbach pipe flow calculation.
#[derive(Debug, Clone)]
pub struct DarcyWeisbach {
    /// Fluid properties
    pub fluid: Fluid,
    /// Inner pipe diameter (m)
    pub diameter_m: f64,
    /// Pipe length (m)
    pub length_m: f64,
    /// Mean flow velocity (m/s)
    pub velocity_m_s: f64,
    /// Pipe wall absolute roughness (m); 0.0 for smooth
    pub roughness_m: f64,
}

/// Results of a Darcy-Weisbach pipe flow calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct PipeFlowResult {
    /// Reynolds number (dimensionless)
    pub reynolds_number: f64,
    /// Flow regime
    pub flow_regime: FlowRegime,
    /// Darcy friction factor (dimensionless)
    pub friction_factor: f64,
    /// Pressure loss over the pipe length (Pa)
    pub pressure_loss_pa: f64,
    /// Head loss (m of fluid)
    pub head_loss_m: f64,
    /// Relative roughness ε/D (dimensionless)
    pub relative_roughness: f64,
}

impl DarcyWeisbach {
    /// Create a new Darcy-Weisbach calculation.
    ///
    /// # Arguments
    /// - `fluid` — fluid properties
    /// - `diameter_m` — inner pipe diameter (m)
    /// - `length_m` — pipe length (m)
    /// - `velocity_m_s` — mean flow velocity (m/s)
    /// - `roughness_m` — absolute pipe roughness (m); use 0.0 for smooth
    pub fn new(
        fluid: &Fluid,
        diameter_m: f64,
        length_m: f64,
        velocity_m_s: f64,
        roughness_m: f64,
    ) -> Self {
        Self {
            fluid: fluid.clone(),
            diameter_m,
            length_m,
            velocity_m_s,
            roughness_m,
        }
    }

    /// Compute Reynolds number.
    #[inline]
    pub fn reynolds_number(&self) -> f64 {
        self.fluid.density_kg_m3 * self.velocity_m_s * self.diameter_m
            / self.fluid.dynamic_viscosity_pa_s
    }

    /// Compute Darcy friction factor.
    ///
    /// - Laminar: f = 64/Re (exact)
    /// - Turbulent: Colebrook-White iterated 10× (converges in ~3-4 iterations)
    /// - Transitional: interpolated between laminar and turbulent values
    pub fn friction_factor(&self, re: f64) -> f64 {
        if re < 2300.0 {
            // Laminar: Hagen-Poiseuille exact
            64.0 / re.max(1e-10)
        } else if re < 4000.0 {
            // Transitional: linear interpolation
            let f_lam = 64.0 / 2300.0;
            let f_turb = self.colebrook_white(4000.0);
            let t = (re - 2300.0) / (4000.0 - 2300.0);
            f_lam + t * (f_turb - f_lam)
        } else {
            self.colebrook_white(re)
        }
    }

    /// Colebrook-White friction factor via fixed-point iteration.
    /// Converges to <1e-10 relative error in ≤10 iterations.
    fn colebrook_white(&self, re: f64) -> f64 {
        let rel_rough = self.roughness_m / self.diameter_m;
        if rel_rough == 0.0 && re > 0.0 {
            // Smooth pipe: use Prandtl-Kármán equation iteratively
            // 1/√f = 2·log10(Re·√f) - 0.8
            return smooth_pipe_friction(re);
        }
        // Swamee-Jain as starting estimate
        let log_arg = rel_rough / 3.7 + 5.74 / pow_f64(re, 0.9);
        let f_init = 0.25 / (log10(log_arg) * log10(log_arg));
        // Colebrook-White iteration: 1/√f = -2·log10(ε/(3.7D) + 2.51/(Re·√f))
        let mut f = f_init;
        for _ in 0..15 {
            let sqrt_f = f.sqrt();
            let rhs = -2.0 * log10(rel_rough / 3.7 + 2.51 / (re * sqrt_f));
            let f_new = 1.0 / (rhs * rhs);
            if (f_new - f).abs() < 1e-12 * f {
                return f_new;
            }
            f = f_new;
        }
        f
    }

    /// Perform the full Darcy-Weisbach calculation.
    pub fn calculate(&self) -> PipeFlowResult {
        let re = self.reynolds_number();
        let regime = FlowRegime::from_reynolds(re);
        let f = self.friction_factor(re);
        let relative_roughness = self.roughness_m / self.diameter_m;

        // ΔP = f · (L/D) · (ρ·v²/2)
        let dynamic_pressure = 0.5 * self.fluid.density_kg_m3 * self.velocity_m_s * self.velocity_m_s;
        let pressure_loss_pa = f * (self.length_m / self.diameter_m) * dynamic_pressure;

        // Head loss = ΔP / (ρ·g)
        let head_loss_m = pressure_loss_pa / (self.fluid.density_kg_m3 * crate::fluid::G);

        PipeFlowResult {
            reynolds_number: re,
            flow_regime: regime,
            friction_factor: f,
            pressure_loss_pa,
            head_loss_m,
            relative_roughness,
        }
    }
}

/// Smooth pipe friction factor via Prandtl-Kármán iteration.
fn smooth_pipe_friction(re: f64) -> f64 {
    // Blasius approximation as initial guess (valid Re 4000-100000)
    let mut f = 0.316 / re.powf(0.25);
    for _ in 0..20 {
        let sqrt_f = f.sqrt();
        let rhs = 2.0 * log10(re * sqrt_f) - 0.8;
        let f_new = 1.0 / (rhs * rhs);
        if (f_new - f).abs() < 1e-12 * f {
            return f_new;
        }
        f = f_new;
    }
    f
}

/// log10.
#[inline]
fn log10(x: f64) -> f64 {
    x.log10()
}

/// x^n for f64 exponent.
#[inline]
fn pow_f64(x: f64, n: f64) -> f64 {
    x.powf(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid::Fluid;

    #[test]
    fn test_laminar_flow_regime() {
        assert_eq!(FlowRegime::from_reynolds(1000.0), FlowRegime::Laminar);
        assert_eq!(FlowRegime::from_reynolds(2299.0), FlowRegime::Laminar);
    }

    #[test]
    fn test_turbulent_flow_regime() {
        assert_eq!(FlowRegime::from_reynolds(4000.0), FlowRegime::Turbulent);
        assert_eq!(FlowRegime::from_reynolds(100_000.0), FlowRegime::Turbulent);
    }

    #[test]
    fn test_laminar_friction_factor() {
        // f = 64/Re for laminar (Re < 2300)
        // water at 20°C: ν ≈ 1.004e-6 m²/s
        // Re = v·D/ν → for Re=500: v = 500 * 1.004e-6 / 0.05 = 0.01004 m/s
        let fluid = Fluid::water(20.0);
        let dw = DarcyWeisbach::new(&fluid, 0.05, 100.0, 0.005, 0.0);
        let re = dw.reynolds_number();
        assert!(re < 2300.0, "Re must be laminar, got {re:.1}");
        let f = dw.friction_factor(re);
        // Exact Hagen-Poiseuille: f = 64/Re
        assert!((f - 64.0 / re).abs() < 1e-10, "laminar f = {f:.6}, expected {:.6}", 64.0 / re);
    }

    #[test]
    fn test_turbulent_reynolds() {
        // Re = ρ·v·D/μ, water at 20°C: ρ=998.2, μ=1.002e-3
        // v=1 m/s, D=0.05m → Re ≈ 49820
        let fluid = Fluid::water(20.0);
        let dw = DarcyWeisbach::new(&fluid, 0.05, 100.0, 1.0, 1.5e-6);
        let re = dw.reynolds_number();
        assert!((re - 49_800.0).abs() < 500.0, "Re = {re:.0}");
    }

    #[test]
    fn test_pressure_loss_positive() {
        let fluid = Fluid::water(20.0);
        let dw = DarcyWeisbach::new(&fluid, 0.05, 100.0, 1.0, 1.5e-6);
        let result = dw.calculate();
        assert!(result.pressure_loss_pa > 0.0);
        assert!(result.head_loss_m > 0.0);
    }

    #[test]
    fn test_moody_chart_known_point() {
        // Moody chart: Re=1e5, ε/D=0.001 → f ≈ 0.022
        let fluid = Fluid::water(20.0);
        // v·D/ν = 1e5, ν = 1.004e-6 → v·D = 0.1004
        // D = 0.1m → v = 1.004 m/s, ε = 1e-4m
        let dw = DarcyWeisbach::new(&fluid, 0.1, 100.0, 1.004, 1.0e-4);
        let re = dw.reynolds_number();
        let f = dw.friction_factor(re);
        assert!((f - 0.022).abs() < 0.003, "f at Re=1e5, ε/D=0.001 = {f:.4}");
    }
}
