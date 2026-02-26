//! Benchmarks for garbongus fluid mechanics calculations.
//!
//! Run with: cargo bench

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use garbongus::{
    capillary::CapillaryRise,
    fluid::Fluid,
    pipe::DarcyWeisbach,
    vacuum::VacuumLift,
};

// ── Fluid property benchmarks ──────────────────────────────────────────────

fn bench_fluid_water(c: &mut Criterion) {
    c.bench_function("fluid::water_properties_20c", |b| {
        b.iter(|| Fluid::water(black_box(20.0)))
    });

    let mut group = c.benchmark_group("fluid::water_by_temperature");
    for temp in [0.0, 20.0, 40.0, 60.0, 80.0, 100.0_f64] {
        group.bench_with_input(BenchmarkId::new("temp_c", temp as i32), &temp, |b, &t| {
            b.iter(|| Fluid::water(black_box(t)))
        });
    }
    group.finish();
}

// ── Capillary rise benchmarks ──────────────────────────────────────────────

fn bench_capillary(c: &mut Criterion) {
    let fluid = Fluid::water(20.0);

    c.bench_function("capillary::jurin_1mm_pipe", |b| {
        let calc = CapillaryRise::new(fluid.clone(), 0.001, 0.0);
        b.iter(|| black_box(calc.calculate()))
    });

    let mut group = c.benchmark_group("capillary::rise_by_radius");
    for radius_mm in [0.1, 0.5, 1.0, 5.0, 10.0_f64] {
        group.bench_with_input(
            BenchmarkId::new("radius_mm", (radius_mm * 10.0) as i32),
            &radius_mm,
            |b, &r| {
                let calc = CapillaryRise::new(fluid.clone(), r / 1000.0, 0.0);
                b.iter(|| black_box(calc.calculate()))
            },
        );
    }
    group.finish();
}

// ── Vacuum lift benchmarks ─────────────────────────────────────────────────

fn bench_vacuum_lift(c: &mut Criterion) {
    let fluid = Fluid::water(20.0);

    c.bench_function("vacuum::lift_5m_static", |b| {
        let lift = VacuumLift::new(fluid.clone(), 0.05, 5.0);
        b.iter(|| black_box(lift.calculate()))
    });

    c.bench_function("vacuum::lift_50m_static", |b| {
        let lift = VacuumLift::new(fluid.clone(), 0.05, 50.0);
        b.iter(|| black_box(lift.calculate()))
    });

    c.bench_function("vacuum::lift_1000m_static", |b| {
        let lift = VacuumLift::new(fluid.clone(), 0.05, 1000.0);
        b.iter(|| black_box(lift.calculate()))
    });

    // With flow velocity (includes friction calculation)
    c.bench_function("vacuum::lift_50m_with_flow_1ms", |b| {
        let lift = VacuumLift::new(fluid.clone(), 0.05, 50.0)
            .flow_velocity(1.0)
            .roughness(1.5e-6);
        b.iter(|| black_box(lift.calculate()))
    });

    // Sweep over distances — demonstrates O(1) per calculation
    let mut group = c.benchmark_group("vacuum::lift_by_distance");
    for height_m in [1.0, 10.0, 100.0, 1_000.0, 10_000.0_f64] {
        group.bench_with_input(
            BenchmarkId::new("height_m", height_m as u64),
            &height_m,
            |b, &h| {
                let lift = VacuumLift::new(fluid.clone(), 0.05, h);
                b.iter(|| black_box(lift.calculate()))
            },
        );
    }
    group.finish();
}

// ── Pipe flow benchmarks ───────────────────────────────────────────────────

fn bench_pipe_flow(c: &mut Criterion) {
    let fluid = Fluid::water(20.0);

    c.bench_function("pipe::darcy_weisbach_laminar", |b| {
        // Re ≈ 100 (laminar)
        let dw = DarcyWeisbach::new(&fluid, 0.05, 100.0, 0.002, 1.5e-6);
        b.iter(|| black_box(dw.calculate()))
    });

    c.bench_function("pipe::darcy_weisbach_turbulent", |b| {
        // Re ≈ 50000 (turbulent)
        let dw = DarcyWeisbach::new(&fluid, 0.05, 100.0, 1.0, 1.5e-6);
        b.iter(|| black_box(dw.calculate()))
    });

    c.bench_function("pipe::darcy_weisbach_smooth_turbulent", |b| {
        let dw = DarcyWeisbach::new(&fluid, 0.05, 100.0, 1.0, 0.0);
        b.iter(|| black_box(dw.calculate()))
    });

    let mut group = c.benchmark_group("pipe::friction_factor_by_reynolds");
    for re in [100.0, 2300.0, 10_000.0, 100_000.0, 1_000_000.0_f64] {
        group.bench_with_input(
            BenchmarkId::new("Re", re as u64),
            &re,
            |b, &r| {
                let dw = DarcyWeisbach::new(&fluid, 0.05, 100.0, r * fluid.dynamic_viscosity_pa_s / (fluid.density_kg_m3 * 0.05), 1.5e-6);
                b.iter(|| black_box(dw.friction_factor(r)))
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_fluid_water,
    bench_capillary,
    bench_vacuum_lift,
    bench_pipe_flow,
);
criterion_main!(benches);
