//! # garbongus
//!
//! Fluid mechanics library for suction/vacuum lift, capillary rise,
//! and pipe flow calculations. Pure Rust, `no_std`-compatible, zero dependencies.
//!
//! ## Modules
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`fluid`] | Fluid property lookups (density, viscosity, surface tension by temperature) |
//! | [`capillary`] | Capillary rise height calculations |
//! | [`vacuum`] | Suction/vacuum lift over arbitrary distances |
//! | [`pipe`] | Darcy-Weisbach friction loss, Reynolds number, flow regime |
//!
//! ## Quick Example
//!
//! ```rust
//! use garbongus::{fluid::Fluid, vacuum::VacuumLift};
//!
//! let fluid = Fluid::water(20.0);
//! let lift = VacuumLift::new(fluid, 0.05, 100.0);
//! let result = lift.calculate();
//! println!("Atm max lift: {:.2} m", result.atmospheric_max_lift_m);
//! println!("Required pump: {:.2} Pa", result.required_pump_pressure_pa);
//! ```

pub mod capillary;
pub mod fluid;
pub mod pipe;
pub mod vacuum;
