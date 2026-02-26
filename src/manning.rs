//! # manning
//!
//! ## Purpose
//! Manning equation for open-channel and gravity-fed tunnel/pipe flow.
//! Used by the Boring Company tunnel hydraulics calculations (TUNNEL.md §3).
//!
//! ## Algorithms
//! - Manning equation (full pipe): Q = (1/n)·A·R^(2/3)·S^(1/2)
//! - Hydraulic radius for full circular pipe: R = D/4
//! - Slope from head difference and length: S = Δh/L
//!
//! ## Data Structures
//! - [`ManningFlow`] — input parameters for Manning calculation
//! - [`ManningResult`] — computed flow velocity, rate, and volume conversions

use core::f64::consts::PI;

/// Input parameters for Manning equation calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ManningFlow {
    /// Pipe/tunnel inner diameter (m)
    pub diameter_m: f64,
    /// Pipe/tunnel length (m)
    pub length_m: f64,
    /// Manning roughness coefficient (dimensionless)
    /// Typical values: concrete 0.012, steel 0.011, HDPE 0.009
    pub roughness_n: f64,
    /// Head difference between inlet and outlet (m)
    /// Positive = downhill (gravity flow), negative = uphill (pumped)
    pub head_difference_m: f64,
}

/// Result of Manning equation calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ManningResult {
    /// Cross-sectional area (m²)
    pub area_m2: f64,
    /// Hydraulic radius (m) — D/4 for full circular pipe
    pub hydraulic_radius_m: f64,
    /// Slope (dimensionless) — head/length
    pub slope: f64,
    /// Flow velocity (m/s)
    pub velocity_m_s: f64,
    /// Volume flow rate (m³/s)
    pub flow_rate_m3_s: f64,
    /// Flow rate in million gallons per day
    pub flow_rate_mgd: f64,
}

impl ManningFlow {
    /// Create a new Manning flow calculation.
    pub fn new(diameter_m: f64, length_m: f64, roughness_n: f64, head_difference_m: f64) -> Self {
        Self {
            diameter_m,
            length_m,
            roughness_n,
            head_difference_m,
        }
    }

    /// Compute Manning equation for full-pipe flow.
    ///
    /// Q = (1/n) · A · R^(2/3) · S^(1/2)
    ///
    /// where:
    /// - n = Manning roughness coefficient
    /// - A = cross-sectional area = π·(D/2)²
    /// - R = hydraulic radius = D/4 (full circular pipe)
    /// - S = slope = |head_difference| / length
    pub fn calculate(&self) -> ManningResult {
        let r = self.diameter_m / 2.0;
        let area = PI * r * r;
        let hydraulic_radius = self.diameter_m / 4.0;
        let slope = self.head_difference_m.abs() / self.length_m;

        // Manning velocity: v = (1/n) · R^(2/3) · S^(1/2)
        let velocity = (1.0 / self.roughness_n)
            * hydraulic_radius.powf(2.0 / 3.0)
            * slope.sqrt();

        let flow_rate = area * velocity;

        // m³/s → MGD: 1 m³/s = 22.8245 MGD
        let flow_rate_mgd = flow_rate * 22.824_5;

        ManningResult {
            area_m2: area,
            hydraulic_radius_m: hydraulic_radius,
            slope,
            velocity_m_s: velocity,
            flow_rate_m3_s: flow_rate,
            flow_rate_mgd,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tunnel_gentle_gravity() {
        // D=3.66m, L=1609m, n=0.012, head=5m (S≈0.31%)
        // Manning full-pipe: v = (1/0.012) × (0.915)^(2/3) × (0.0031)^(1/2) ≈ 4.3 m/s
        // Q = 10.52 × 4.3 ≈ 45.5 m³/s
        let mf = ManningFlow::new(3.66, 1609.0, 0.012, 5.0);
        let r = mf.calculate();
        assert!((r.slope - 0.0031).abs() < 0.001, "slope = {:.4}", r.slope);
        assert!((r.velocity_m_s - 4.3).abs() < 0.5, "v = {:.2} m/s", r.velocity_m_s);
        assert!(r.flow_rate_m3_s > 30.0 && r.flow_rate_m3_s < 60.0,
            "Q = {:.1} m³/s", r.flow_rate_m3_s);
    }

    #[test]
    fn test_tunnel_pumped_high() {
        // D=3.66m, L=1609m, n=0.012, head=50m (S≈3.1%)
        // Manning full-pipe: v = (1/0.012) × (0.915)^(2/3) × (0.031)^(1/2) ≈ 13.8 m/s
        // Q = 10.52 × 13.8 ≈ 145 m³/s
        let mf = ManningFlow::new(3.66, 1609.0, 0.012, 50.0);
        let r = mf.calculate();
        assert!((r.slope - 0.031).abs() < 0.002, "slope = {:.4}", r.slope);
        assert!((r.velocity_m_s - 13.8).abs() < 1.0, "v = {:.2} m/s", r.velocity_m_s);
        assert!(r.flow_rate_m3_s > 100.0 && r.flow_rate_m3_s < 200.0,
            "Q = {:.1} m³/s", r.flow_rate_m3_s);
    }

    #[test]
    fn test_hydraulic_radius_full_pipe() {
        // Full circular pipe: R = D/4
        let mf = ManningFlow::new(4.0, 100.0, 0.012, 1.0);
        let r = mf.calculate();
        assert!((r.hydraulic_radius_m - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zero_slope_zero_flow() {
        let mf = ManningFlow::new(3.66, 1609.0, 0.012, 0.0);
        let r = mf.calculate();
        assert!(r.velocity_m_s.abs() < 1e-10, "zero slope → zero velocity");
        assert!(r.flow_rate_m3_s.abs() < 1e-10, "zero slope → zero flow");
    }

    #[test]
    fn test_rougher_pipe_slower_flow() {
        let smooth = ManningFlow::new(3.66, 1609.0, 0.009, 10.0).calculate();
        let rough = ManningFlow::new(3.66, 1609.0, 0.020, 10.0).calculate();
        assert!(smooth.velocity_m_s > rough.velocity_m_s, "smoother pipe → faster flow");
    }

    #[test]
    fn test_tucson_tunnel_velocity() {
        // Tucson tunnel: D=3.66m, L=1287m, head=42.6m (S≈3.3%)
        // Manning full-pipe gives high velocity (~14 m/s) — the tunnel would be
        // throttled or partially filled. Verify Manning math is internally consistent:
        // v = (1/n) × R^(2/3) × S^(1/2)
        let mf = ManningFlow::new(3.66, 1287.0, 0.012, 42.6);
        let r = mf.calculate();
        // Verify: v = Q / A
        let expected_v = r.flow_rate_m3_s / r.area_m2;
        assert!((r.velocity_m_s - expected_v).abs() < 1e-10,
            "Manning velocity consistent: {:.4} vs {:.4}", r.velocity_m_s, expected_v);
        // At this slope, full-pipe velocity is high (>10 m/s)
        assert!(r.velocity_m_s > 10.0, "high slope → high Manning velocity = {:.2}", r.velocity_m_s);
    }
}
