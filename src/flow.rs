//! # flow
//!
//! ## Purpose
//! Flow rate calculations, pipe sizing, unit conversions, and Bernoulli's equation.
//! Provides the core hydraulic relationships used throughout pipeline and tunnel design.
//!
//! ## Algorithms
//! - Continuity: Q = A·v, A = π·(D/2)²
//! - Pipe sizing: D = 2·√(Q / (v·π))
//! - Bernoulli: P₁ + ½ρv₁² + ρgh₁ = P₂ + ½ρv₂² + ρgh₂
//! - Pump power: P = ρ·g·Q·H / η
//! - Pressure→depth: d = P / (ρ·g) (well meter sensor conversion)
//!
//! ## Data Structures
//! - [`FlowRate`] — volume flow rate with unit conversions
//! - [`PumpPower`] — pump power calculation result

use core::f64::consts::PI;
use crate::fluid::G;

// ---------------------------------------------------------------------------
// Unit conversion constants
// ---------------------------------------------------------------------------

/// Cubic meters per second → million gallons per day
const M3S_TO_MGD: f64 = 22.824_5;

/// Cubic meters per second → million liters per day
const M3S_TO_MLD: f64 = 86.4;

/// Cubic meters per second → liters per minute
const M3S_TO_LPM: f64 = 60_000.0;

/// Cubic meters → acre-feet
const M3_TO_AF: f64 = 1.0 / 1233.48;

/// Seconds per year (365.25 days)
const SECONDS_PER_YEAR: f64 = 365.25 * 24.0 * 3600.0;

// ---------------------------------------------------------------------------
// Flow rate struct with unit conversions
// ---------------------------------------------------------------------------

/// Volume flow rate with built-in unit conversions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlowRate {
    /// Flow rate in m³/s (SI base unit)
    pub m3_per_s: f64,
}

impl FlowRate {
    /// Create from m³/s.
    #[inline]
    pub fn from_m3s(q: f64) -> Self {
        Self { m3_per_s: q }
    }

    /// Create from liters per minute.
    #[inline]
    pub fn from_lpm(lpm: f64) -> Self {
        Self { m3_per_s: lpm / M3S_TO_LPM }
    }

    /// Create from million gallons per day.
    #[inline]
    pub fn from_mgd(mgd: f64) -> Self {
        Self { m3_per_s: mgd / M3S_TO_MGD }
    }

    /// Convert to million gallons per day (MGD).
    #[inline]
    pub fn to_mgd(&self) -> f64 {
        self.m3_per_s * M3S_TO_MGD
    }

    /// Convert to million liters per day (MLD).
    #[inline]
    pub fn to_mld(&self) -> f64 {
        self.m3_per_s * M3S_TO_MLD
    }

    /// Convert to liters per minute (L/min).
    #[inline]
    pub fn to_lpm(&self) -> f64 {
        self.m3_per_s * M3S_TO_LPM
    }

    /// Annual volume in m³.
    #[inline]
    pub fn annual_volume_m3(&self) -> f64 {
        self.m3_per_s * SECONDS_PER_YEAR
    }

    /// Annual volume in acre-feet.
    #[inline]
    pub fn annual_volume_af(&self) -> f64 {
        self.annual_volume_m3() * M3_TO_AF
    }
}

// ---------------------------------------------------------------------------
// Pipe geometry and continuity
// ---------------------------------------------------------------------------

/// Cross-sectional area of a circular pipe (m²).
/// A = π·(D/2)²
#[inline]
pub fn pipe_area(diameter_m: f64) -> f64 {
    PI * (diameter_m / 2.0) * (diameter_m / 2.0)
}

/// Volume flow rate from pipe diameter and velocity.
/// Q = A·v (m³/s)
#[inline]
pub fn flow_rate(diameter_m: f64, velocity_m_s: f64) -> f64 {
    pipe_area(diameter_m) * velocity_m_s
}

/// Velocity from pipe diameter and flow rate.
/// v = Q / A (m/s)
#[inline]
pub fn velocity(diameter_m: f64, flow_rate_m3s: f64) -> f64 {
    flow_rate_m3s / pipe_area(diameter_m)
}

/// Required pipe diameter for a target flow rate and velocity.
/// D = 2·√(Q / (v·π))
#[inline]
pub fn required_diameter(flow_rate_m3s: f64, velocity_m_s: f64) -> f64 {
    let area = flow_rate_m3s / velocity_m_s;
    2.0 * (area / PI).sqrt()
}

/// Mass flow rate (kg/s).
/// ṁ = ρ·Q = ρ·A·v
#[inline]
pub fn mass_flow_rate(density_kg_m3: f64, flow_rate_m3s: f64) -> f64 {
    density_kg_m3 * flow_rate_m3s
}

// ---------------------------------------------------------------------------
// Bernoulli's equation
// ---------------------------------------------------------------------------

/// Bernoulli pressure at point 2 given conditions at point 1.
///
/// P₁ + ½ρv₁² + ρgh₁ = P₂ + ½ρv₂² + ρgh₂
///
/// Solves for P₂:
/// P₂ = P₁ + ½ρ(v₁² - v₂²) + ρg(h₁ - h₂)
#[inline]
pub fn bernoulli_pressure(
    p1_pa: f64,
    v1_m_s: f64,
    h1_m: f64,
    v2_m_s: f64,
    h2_m: f64,
    density_kg_m3: f64,
) -> f64 {
    p1_pa
        + 0.5 * density_kg_m3 * (v1_m_s * v1_m_s - v2_m_s * v2_m_s)
        + density_kg_m3 * G * (h1_m - h2_m)
}

// ---------------------------------------------------------------------------
// Pump power
// ---------------------------------------------------------------------------

/// Pump power calculation result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PumpPower {
    /// Hydraulic power (W) — ideal power without efficiency losses
    pub hydraulic_w: f64,
    /// Shaft power (W) — actual power accounting for pump efficiency
    pub shaft_w: f64,
    /// Pump efficiency used (dimensionless, 0–1)
    pub efficiency: f64,
}

impl PumpPower {
    /// Shaft power in kilowatts.
    #[inline]
    pub fn shaft_kw(&self) -> f64 {
        self.shaft_w / 1e3
    }

    /// Shaft power in megawatts.
    #[inline]
    pub fn shaft_mw(&self) -> f64 {
        self.shaft_w / 1e6
    }

    /// Shaft power in gigawatts.
    #[inline]
    pub fn shaft_gw(&self) -> f64 {
        self.shaft_w / 1e9
    }
}

/// Calculate pump power required to move fluid.
///
/// P = ρ·g·Q·H / η
///
/// # Arguments
/// - `density_kg_m3` — fluid density (kg/m³)
/// - `flow_rate_m3s` — volume flow rate (m³/s)
/// - `total_head_m` — total dynamic head (m) — elevation + friction
/// - `efficiency` — pump efficiency (0–1, typically 0.70–0.90)
pub fn pump_power(
    density_kg_m3: f64,
    flow_rate_m3s: f64,
    total_head_m: f64,
    efficiency: f64,
) -> PumpPower {
    let hydraulic_w = density_kg_m3 * G * flow_rate_m3s * total_head_m;
    let shaft_w = hydraulic_w / efficiency;
    PumpPower {
        hydraulic_w,
        shaft_w,
        efficiency,
    }
}

// ---------------------------------------------------------------------------
// Well meter sensor conversion
// ---------------------------------------------------------------------------

/// Convert pressure to water column depth (m).
/// d = P / (ρ·g)
///
/// Used by WELL_METER pressure transducers to convert gauge pressure to water level.
#[inline]
pub fn pressure_to_depth(pressure_pa: f64, density_kg_m3: f64) -> f64 {
    pressure_pa / (density_kg_m3 * G)
}

/// Convert water column depth to pressure (Pa).
/// P = ρ·g·d
#[inline]
pub fn depth_to_pressure(depth_m: f64, density_kg_m3: f64) -> f64 {
    density_kg_m3 * G * depth_m
}

// ---------------------------------------------------------------------------
// Water balance / required flow rate
// ---------------------------------------------------------------------------

/// Calculate required replenishment flow rate to offset annual depletion.
///
/// Q = (V_loss_per_year / recharge_efficiency) / seconds_per_year
///
/// # Arguments
/// - `annual_loss_m3` — annual volume lost (m³/year)
/// - `recharge_efficiency` — fraction of injected water that reaches aquifer (0–1)
///
/// # Returns
/// Required continuous flow rate in m³/s.
#[inline]
pub fn required_replenishment_flow(annual_loss_m3: f64, recharge_efficiency: f64) -> f64 {
    (annual_loss_m3 / recharge_efficiency) / SECONDS_PER_YEAR
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Pipe geometry ---

    #[test]
    fn test_pipe_area_1m() {
        let a = pipe_area(1.0);
        assert!((a - PI / 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_flow_rate_continuity() {
        // Q = A·v, D=0.05m, v=3 m/s
        let q = flow_rate(0.05, 3.0);
        let expected = PI * 0.025 * 0.025 * 3.0;
        assert!((q - expected).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_from_flow() {
        let q = 1.0; // m³/s
        let d = 1.0; // m
        let v = velocity(d, q);
        // v = Q/A = 1.0 / (π/4) ≈ 1.2732
        assert!((v - 4.0 / PI).abs() < 1e-10);
    }

    // --- Pipe diameter sizing (WATER README §Pipe Specs) ---

    #[test]
    fn test_required_diameter_igbwp() {
        // IGBWP: Q=595 m³/s, v=3 m/s → D ≈ 15.9 m
        let d = required_diameter(595.0, 3.0);
        assert!((d - 15.9).abs() < 0.5, "IGBWP diameter = {d:.1}");
    }

    #[test]
    fn test_required_diameter_california() {
        // README: Q=1172 m³/s, v=3 m/s → D ≈ 22.3 m
        let d = required_diameter(1172.0, 3.0);
        assert!((d - 22.3).abs() < 0.5, "California diameter = {d:.1}");
    }

    // --- Unit conversions ---

    #[test]
    fn test_flow_rate_mgd_conversion() {
        let fr = FlowRate::from_m3s(1.0);
        assert!((fr.to_mgd() - 22.8245).abs() < 0.01);
    }

    #[test]
    fn test_flow_rate_mld_conversion() {
        let fr = FlowRate::from_m3s(1.0);
        assert!((fr.to_mld() - 86.4).abs() < 0.1);
    }

    #[test]
    fn test_flow_rate_lpm_conversion() {
        let fr = FlowRate::from_m3s(1.0);
        assert!((fr.to_lpm() - 60_000.0).abs() < 0.1);
    }

    #[test]
    fn test_flow_rate_roundtrip_mgd() {
        let fr = FlowRate::from_mgd(64.0);
        assert!((fr.to_mgd() - 64.0).abs() < 0.01);
    }

    #[test]
    fn test_tunnel_flow_mgd() {
        // TUNNEL.md: 2.8 m³/s = 64 MGD
        let fr = FlowRate::from_m3s(2.8);
        assert!((fr.to_mgd() - 64.0).abs() < 1.0, "tunnel MGD = {:.1}", fr.to_mgd());
    }

    // --- Bernoulli ---

    #[test]
    fn test_bernoulli_static_column() {
        // Static water column: v1=v2=0, h1=10m, h2=0m, P1=0 gauge
        // P2 = P1 + ρg(h1-h2) = 0 + 1000*9.80665*10 = 98066.5 Pa
        let p2 = bernoulli_pressure(0.0, 0.0, 10.0, 0.0, 0.0, 1000.0);
        assert!((p2 - 98_066.5).abs() < 1.0, "P2 = {p2:.1}");
    }

    #[test]
    fn test_bernoulli_venturi() {
        // Narrowing pipe: P drops as v increases
        // P1=200kPa, v1=2 m/s, h1=h2=0, v2=4 m/s, ρ=1000
        // P2 = 200000 + 0.5*1000*(4 - 16) = 200000 - 6000 = 194000
        let p2 = bernoulli_pressure(200_000.0, 2.0, 0.0, 4.0, 0.0, 1000.0);
        assert!((p2 - 194_000.0).abs() < 1.0, "P2 = {p2:.1}");
    }

    // --- Pump power (IGBWP §7) ---

    #[test]
    fn test_pump_power_igbwp_elevation() {
        // P = ρ·g·Q·H / η = 1000 × 9.80665 × 595 × 300 / 0.85 ≈ 2.06 GW
        // Note: IGBWP §7 quotes 4.1 GW using ρ≈1025 (seawater) and possibly different η.
        // Our calculation with fresh water ρ=1000 is correct at ~2.06 GW.
        let pp = pump_power(1000.0, 595.0, 300.0, 0.85);
        assert!((pp.shaft_gw() - 2.06).abs() < 0.1, "elevation pump = {:.2} GW", pp.shaft_gw());
    }

    #[test]
    fn test_pump_power_california() {
        // README: Q=1172, H=500m, η=0.85 → ~6.76 GW
        let pp = pump_power(1000.0, 1172.0, 500.0, 0.85);
        assert!((pp.shaft_gw() - 6.76).abs() < 0.5, "California pump = {:.2} GW", pp.shaft_gw());
    }

    // --- Well meter sensor ---

    #[test]
    fn test_pressure_to_depth_well_meter() {
        // WELL_METER §5: P=49050 Pa → d = 49050/(1000×9.81) ≈ 5.0 m
        let d = pressure_to_depth(49_050.0, 1000.0);
        assert!((d - 5.0).abs() < 0.05, "well depth = {d:.2} m");
    }

    #[test]
    fn test_depth_to_pressure_roundtrip() {
        let depth = 25.0;
        let p = depth_to_pressure(depth, 997.0);
        let d2 = pressure_to_depth(p, 997.0);
        assert!((d2 - depth).abs() < 1e-10);
    }

    // --- Flow rate from well meter DN40 ---

    #[test]
    fn test_well_meter_dn40_flow() {
        // WELL_METER §5: DN40 pipe, A=π×0.02², v=2 m/s → Q = 2.51 L/s = 151 L/min
        let q = flow_rate(0.04, 2.0);
        let fr = FlowRate::from_m3s(q);
        assert!((fr.to_lpm() - 151.0).abs() < 1.0, "DN40 flow = {:.1} L/min", fr.to_lpm());
    }

    // --- Required replenishment flow ---

    #[test]
    fn test_required_flow_igbwp() {
        // IGBWP: 15 km³/year overdraft, 80% recharge efficiency → ~595 m³/s
        let q = required_replenishment_flow(15e9, 0.80);
        assert!((q - 595.0).abs() < 5.0, "IGBWP required Q = {q:.1} m³/s");
    }

    #[test]
    fn test_required_flow_california() {
        // README: V₀=3.7e12, R=1%/yr → annual_loss = 3.7e10 → Q ≈ 1172 m³/s
        // README assumes 100% efficiency (no recharge factor mentioned)
        let q = required_replenishment_flow(3.7e10, 1.0);
        assert!((q - 1172.0).abs() < 5.0, "California required Q = {q:.1} m³/s");
    }
}
