//! # vacuum
//!
//! ## Purpose
//! Suction/vacuum lift calculations for moving water up a pipe over any
//! vertical distance. Models the pressure balance between atmospheric
//! pressure driving the fluid up, gravity pulling it down, friction losses,
//! and the vapor pressure limit that causes cavitation.
//!
//! ## Physics
//!
//! The fundamental equation for suction lift:
//!
//! ```text
//! P_atm = P_vapor + ρ·g·h + ΔP_friction + P_vacuum_applied
//! ```
//!
//! Rearranging for maximum achievable lift height:
//!
//! ```text
//! h_max = (P_atm - P_vapor - ΔP_friction) / (ρ·g)
//! ```
//!
//! The **theoretical absolute maximum** with no friction and at 20°C:
//! ```text
//! h_max = (101325 - 2337) / (998.2 * 9.80665) ≈ 10.18 m
//! ```
//!
//! With a vacuum pump applying additional pressure differential, the lift
//! can extend **beyond** the atmospheric limit by the pump pressure:
//!
//! ```text
//! h_total = h_atm_max + P_pump / (ρ·g)
//! ```
//!
//! For any distance, the required vacuum pump pressure is:
//! ```text
//! P_required = ρ·g·h + ΔP_friction + P_vapor - P_atm
//! ```
//! (positive value means a pump is needed beyond atmospheric)
//!
//! ## Data Structures
//! - [`VacuumLift`] — input parameters
//! - [`VacuumResult`] — comprehensive output with all pressures and limits

use crate::fluid::{Fluid, G, P_ATM};
use crate::pipe::DarcyWeisbach;

/// Input parameters for a vacuum lift system.
#[derive(Debug, Clone)]
pub struct VacuumLift {
    /// Fluid being lifted
    pub fluid: Fluid,
    /// Inner pipe radius (m)
    pub pipe_radius_m: f64,
    /// Target lift height (m) — vertical distance to raise fluid
    pub target_height_m: f64,
    /// Pipe wall roughness (m); 0.0 for smooth pipe (default: 1.5e-6 for commercial steel)
    pub roughness_m: f64,
    /// Flow velocity (m/s); 0.0 means static (no flow friction)
    pub flow_velocity_m_s: f64,
    /// Ambient atmospheric pressure (Pa); defaults to sea-level 101325 Pa
    pub ambient_pressure_pa: f64,
}

/// Comprehensive results of a vacuum lift calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct VacuumResult {
    /// Target lift height (m)
    pub target_height_m: f64,
    /// Maximum lift achievable using atmosphere alone (no pump) (m)
    pub atmospheric_max_lift_m: f64,
    /// Maximum theoretical lift if perfect vacuum applied (m)
    pub theoretical_max_lift_m: f64,
    /// Hydrostatic pressure of the target column (Pa)
    pub hydrostatic_pressure_pa: f64,
    /// Friction pressure loss over the pipe length (Pa)
    pub friction_loss_pa: f64,
    /// Vapor pressure of the fluid at its temperature (Pa)
    pub vapor_pressure_pa: f64,
    /// Net available atmospheric pressure margin after vapor pressure (Pa)
    pub net_atmospheric_pa: f64,
    /// Whether the target height is achievable with atmosphere alone
    pub achievable_by_atmosphere: bool,
    /// Whether the target height is achievable at all (below cavitation limit)
    pub achievable_with_pump: bool,
    /// Vacuum pump pressure required to reach target height (Pa)
    /// 0.0 if achievable by atmosphere alone
    pub required_pump_pressure_pa: f64,
    /// Total suction pressure at pipe inlet (Pa absolute)
    pub suction_pressure_pa_abs: f64,
    /// Reynolds number of the flow
    pub reynolds_number: f64,
    /// Darcy friction factor
    pub friction_factor: f64,
}

impl VacuumLift {
    /// Create a new vacuum lift calculation.
    ///
    /// # Arguments
    /// - `fluid` — fluid properties (use [`crate::fluid::Fluid::water`])
    /// - `pipe_radius_m` — inner pipe radius (m)
    /// - `target_height_m` — vertical distance to lift fluid (m); can be any positive value
    pub fn new(fluid: Fluid, pipe_radius_m: f64, target_height_m: f64) -> Self {
        Self {
            fluid,
            pipe_radius_m,
            target_height_m,
            roughness_m: 1.5e-6,
            flow_velocity_m_s: 0.0,
            ambient_pressure_pa: P_ATM,
        }
    }

    /// Set pipe wall roughness (m). Default: 1.5e-6 m (commercial steel).
    pub fn roughness(mut self, roughness_m: f64) -> Self {
        self.roughness_m = roughness_m;
        self
    }

    /// Set flow velocity (m/s). Default: 0.0 (static, no friction loss).
    pub fn flow_velocity(mut self, velocity_m_s: f64) -> Self {
        self.flow_velocity_m_s = velocity_m_s;
        self
    }

    /// Set ambient atmospheric pressure (Pa). Default: 101325 Pa (sea level).
    /// Use lower values for high-altitude installations.
    pub fn ambient_pressure(mut self, pressure_pa: f64) -> Self {
        self.ambient_pressure_pa = pressure_pa;
        self
    }

    /// Perform the vacuum lift calculation.
    pub fn calculate(&self) -> VacuumResult {
        let rho = self.fluid.density_kg_m3;
        let p_vap = self.fluid.vapor_pressure_pa;
        let h = self.target_height_m;
        let diameter_m = 2.0 * self.pipe_radius_m;

        // Compute friction loss if flow velocity is non-zero
        let (friction_loss_pa, reynolds_number, friction_factor) =
            if self.flow_velocity_m_s > 0.0 {
                let dw = DarcyWeisbach::new(
                    &self.fluid,
                    diameter_m,
                    h, // pipe length approximated as lift height
                    self.flow_velocity_m_s,
                    self.roughness_m,
                );
                let res = dw.calculate();
                (res.pressure_loss_pa, res.reynolds_number, res.friction_factor)
            } else {
                (0.0, 0.0, 0.0)
            };

        // Net atmospheric pressure available (atmosphere minus vapor pressure)
        let net_atmospheric_pa = self.ambient_pressure_pa - p_vap;

        // Hydrostatic pressure for target column
        let hydrostatic_pressure_pa = rho * G * h;

        // Maximum lift by atmosphere alone (static, no friction)
        let atmospheric_max_lift_m = net_atmospheric_pa / (rho * G);

        // Theoretical max = full vacuum applied (P_abs = 0 at inlet)
        // Limited only by vapor pressure to prevent cavitation
        // With pump: h_max = (P_pump + P_atm - P_vapor) / (ρ·g)
        // There is no physical upper limit with a pump; we report the unconstrained value
        let theoretical_max_lift_m = f64::INFINITY;

        // Total pressure needed at bottom to push fluid to height h
        let total_required_pa = hydrostatic_pressure_pa + friction_loss_pa + p_vap;

        // What the atmosphere provides
        let atmosphere_provides_pa = self.ambient_pressure_pa;

        // Surplus = what atmosphere gives minus what's needed
        let surplus_pa = atmosphere_provides_pa - total_required_pa;

        let achievable_by_atmosphere = surplus_pa >= 0.0;

        // With a pump, it's always achievable (no physical upper bound)
        // as long as we don't exceed absolute vacuum (can't go below 0 Pa absolute)
        // The suction side pressure = P_atm - ρgh - friction - pump_pressure
        // For cavitation: suction_pressure > P_vapor
        // So pump can add unlimited pressure on discharge side — achievable_with_pump always true
        let achievable_with_pump = true;

        // Required pump pressure: how much beyond atmosphere do we need?
        let required_pump_pressure_pa = if achievable_by_atmosphere {
            0.0
        } else {
            -surplus_pa // positive value = pump must supply this
        };

        // Absolute pressure at suction inlet (pipe bottom)
        // P_suction = P_atm + P_pump - ρ·g·h - friction
        let suction_pressure_pa_abs =
            self.ambient_pressure_pa + required_pump_pressure_pa
                - hydrostatic_pressure_pa
                - friction_loss_pa;

        VacuumResult {
            target_height_m: h,
            atmospheric_max_lift_m,
            theoretical_max_lift_m,
            hydrostatic_pressure_pa,
            friction_loss_pa,
            vapor_pressure_pa: p_vap,
            net_atmospheric_pa,
            achievable_by_atmosphere,
            achievable_with_pump,
            required_pump_pressure_pa,
            suction_pressure_pa_abs,
            reynolds_number,
            friction_factor,
        }
    }
}

impl VacuumResult {
    /// Required pump pressure in bar (1 bar = 100,000 Pa).
    #[inline]
    pub fn required_pump_pressure_bar(&self) -> f64 {
        self.required_pump_pressure_pa / 100_000.0
    }

    /// Required pump pressure in PSI (1 PSI ≈ 6894.76 Pa).
    #[inline]
    pub fn required_pump_pressure_psi(&self) -> f64 {
        self.required_pump_pressure_pa / 6_894.757
    }

    /// Required pump pressure in meters of water head.
    #[inline]
    pub fn required_pump_head_m(&self, density_kg_m3: f64) -> f64 {
        self.required_pump_pressure_pa / (density_kg_m3 * crate::fluid::G)
    }

    /// Atmospheric max lift in feet.
    #[inline]
    pub fn atmospheric_max_lift_ft(&self) -> f64 {
        self.atmospheric_max_lift_m * 3.280_84
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid::Fluid;

    #[test]
    fn test_atmospheric_max_lift_20c() {
        // At 20°C, atmospheric max lift ≈ 10.18 m
        let fluid = Fluid::water(20.0);
        let lift = VacuumLift::new(fluid, 0.05, 5.0);
        let result = lift.calculate();
        assert!(
            (result.atmospheric_max_lift_m - 10.18).abs() < 0.2,
            "atm max lift = {:.3} m",
            result.atmospheric_max_lift_m
        );
    }

    #[test]
    fn test_5m_lift_achievable_by_atmosphere() {
        let fluid = Fluid::water(20.0);
        let lift = VacuumLift::new(fluid, 0.05, 5.0);
        let result = lift.calculate();
        assert!(result.achievable_by_atmosphere, "5m should be achievable by atmosphere");
        assert_eq!(result.required_pump_pressure_pa, 0.0);
    }

    #[test]
    fn test_50m_lift_requires_pump() {
        let fluid = Fluid::water(20.0);
        let lift = VacuumLift::new(fluid, 0.05, 50.0);
        let result = lift.calculate();
        assert!(!result.achievable_by_atmosphere, "50m requires a pump");
        assert!(result.achievable_with_pump, "50m is achievable with pump");
        assert!(result.required_pump_pressure_pa > 0.0);
    }

    #[test]
    fn test_1000m_lift_requires_pump() {
        // Any distance is achievable with enough pump pressure
        let fluid = Fluid::water(20.0);
        let lift = VacuumLift::new(fluid, 0.05, 1000.0);
        let result = lift.calculate();
        assert!(!result.achievable_by_atmosphere);
        assert!(result.achievable_with_pump);
        // Required pump pressure for 1000m ≈ ρ·g·h ≈ 998.2 * 9.80665 * 1000 ≈ 9.79 MPa
        assert!(
            (result.required_pump_pressure_pa - 9_79_0000.0).abs() < 100_000.0,
            "pump pressure for 1000m = {:.0} Pa",
            result.required_pump_pressure_pa
        );
    }

    #[test]
    fn test_pump_head_conversion() {
        let fluid = Fluid::water(20.0);
        let lift = VacuumLift::new(fluid.clone(), 0.05, 50.0);
        let result = lift.calculate();
        // Head should be approximately (target - atm_max) meters
        let head = result.required_pump_head_m(fluid.density_kg_m3);
        assert!((head - (50.0 - result.atmospheric_max_lift_m)).abs() < 1.0);
    }

    #[test]
    fn test_high_altitude_reduces_max_lift() {
        // At 2500m altitude, P_atm ≈ 74,682 Pa — less lift available
        let fluid_sea = Fluid::water(20.0);
        let fluid_alt = Fluid::water(20.0);
        let sea_result = VacuumLift::new(fluid_sea, 0.05, 5.0).calculate();
        let alt_result = VacuumLift::new(fluid_alt, 0.05, 5.0)
            .ambient_pressure(74_682.0)
            .calculate();
        assert!(
            alt_result.atmospheric_max_lift_m < sea_result.atmospheric_max_lift_m,
            "altitude reduces max lift"
        );
    }
}
