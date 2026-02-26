//! # garbongus
//!
//! Fluid mechanics library for pipeline hydraulics, desalination water transport,
//! tunnel flow, suction/vacuum lift, capillary rise, and pipe flow calculations.
//! Supports fresh water and seawater/brine. Pure Rust, zero dependencies.
//!
//! ## Modules
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`fluid`] | Fluid properties: fresh water + seawater/brine (density, viscosity, surface tension, vapor pressure) |
//! | [`flow`] | Flow rate, pipe sizing, unit conversions, Bernoulli, pump power, pressure↔depth |
//! | [`manning`] | Manning equation for tunnel and open-channel gravity flow |
//! | [`pipeline`] | Multi-segment pipeline with elevation profile, friction losses, pump power |
//! | [`pipe`] | Darcy-Weisbach friction loss, Colebrook-White, Reynolds number, flow regime |
//! | [`vacuum`] | Suction/vacuum lift over arbitrary distances |
//! | [`capillary`] | Capillary rise height calculations |
//!
//! ## Quick Example
//!
//! ```rust
//! use garbongus::fluid::Fluid;
//! use garbongus::flow::{required_diameter, pump_power, FlowRate};
//!
//! // Seawater at 35 ppt salinity, 20°C
//! let fluid = Fluid::seawater(20.0, 35.0);
//! println!("Density: {:.1} kg/m³", fluid.density_kg_m3);
//!
//! // Size a pipe: Q = 10 m³/s at 2 m/s → D ≈ 2.52 m
//! let d = required_diameter(10.0, 2.0);
//! println!("Required diameter: {:.2} m", d);
//!
//! // Pump power for 50m total head, 85% efficiency
//! let pp = pump_power(fluid.density_kg_m3, 10.0, 50.0, 0.85);
//! println!("Pump power: {:.1} kW", pp.shaft_kw());
//!
//! // Flow rate conversions
//! let fr = FlowRate::from_m3s(1.0);
//! println!("{:.1} MGD = {:.0} L/min", fr.to_mgd(), fr.to_lpm());
//! ```

pub mod capillary;
pub mod fluid;
pub mod flow;
pub mod manning;
pub mod pipe;
pub mod pipeline;
pub mod vacuum;
