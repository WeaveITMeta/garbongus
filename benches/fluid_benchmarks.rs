//! Benchmarks for garbongus fluid mechanics calculations.
//!
//! Run with: cargo bench

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use garbongus::{
    capillary::CapillaryRise,
    fluid::Fluid,
    flow,
    manning::ManningFlow,
    pipe::DarcyWeisbach,
    pipeline::{Pipeline, PipelineSegment},
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

// ── Seawater property benchmarks ─────────────────────────────────────────

fn bench_fluid_seawater(c: &mut Criterion) {
    c.bench_function("fluid::seawater_35ppt_20c", |b| {
        b.iter(|| Fluid::seawater(black_box(20.0), black_box(35.0)))
    });

    let mut group = c.benchmark_group("fluid::seawater_by_salinity");
    for s in [0.0, 10.0, 32.0, 35.0, 50.0, 70.0_f64] {
        group.bench_with_input(BenchmarkId::new("ppt", s as i32), &s, |b, &sal| {
            b.iter(|| Fluid::seawater(black_box(25.0), black_box(sal)))
        });
    }
    group.finish();
}

// ── Flow & pump power benchmarks ─────────────────────────────────────────

fn bench_flow(c: &mut Criterion) {
    c.bench_function("flow::required_diameter_large", |b| {
        b.iter(|| flow::required_diameter(black_box(595.0), black_box(3.0)))
    });

    c.bench_function("flow::pump_power_high_flow", |b| {
        b.iter(|| flow::pump_power(black_box(1025.0), black_box(595.0), black_box(300.0), black_box(0.85)))
    });

    c.bench_function("flow::bernoulli_pressure", |b| {
        b.iter(|| flow::bernoulli_pressure(
            black_box(101_325.0), black_box(2.0), black_box(10.0),
            black_box(4.0), black_box(0.0), black_box(1000.0),
        ))
    });

    c.bench_function("flow::pressure_to_depth", |b| {
        b.iter(|| flow::pressure_to_depth(black_box(49_050.0), black_box(1000.0)))
    });
}

// ── Manning equation benchmarks ──────────────────────────────────────────

fn bench_manning(c: &mut Criterion) {
    c.bench_function("manning::3_66m_gentle_slope", |b| {
        let mf = ManningFlow::new(3.66, 1609.0, 0.012, 5.0);
        b.iter(|| black_box(mf.calculate()))
    });

    c.bench_function("manning::3_66m_steep_slope", |b| {
        let mf = ManningFlow::new(3.66, 1609.0, 0.012, 50.0);
        b.iter(|| black_box(mf.calculate()))
    });
}

// ── Pipeline benchmarks ──────────────────────────────────────────────────

fn bench_pipeline(c: &mut Criterion) {
    c.bench_function("pipeline::4_segment_750km", |b| {
        let fluid = Fluid::water(25.0);
        let mut pl = Pipeline::new(fluid, 595.0, 0.85);
        pl.add_segment(PipelineSegment {
            name: "Coast→Foothills".into(),
            horizontal_distance_m: 200_000.0,
            start_elevation_m: 0.0, end_elevation_m: 170.0,
            terrain_multiplier: 1.10, diameter_m: 16.0, roughness_m: 0.000_046,
        });
        pl.add_segment(PipelineSegment {
            name: "Foothills→Ridge".into(),
            horizontal_distance_m: 180_000.0,
            start_elevation_m: 170.0, end_elevation_m: 300.0,
            terrain_multiplier: 1.20, diameter_m: 16.0, roughness_m: 0.000_046,
        });
        pl.add_segment(PipelineSegment {
            name: "Ridge→Valley".into(),
            horizontal_distance_m: 200_000.0,
            start_elevation_m: 300.0, end_elevation_m: 80.0,
            terrain_multiplier: 1.15, diameter_m: 16.0, roughness_m: 0.000_046,
        });
        pl.add_segment(PipelineSegment {
            name: "Valley→Destination".into(),
            horizontal_distance_m: 80_000.0,
            start_elevation_m: 80.0, end_elevation_m: 98.0,
            terrain_multiplier: 1.05, diameter_m: 16.0, roughness_m: 0.000_046,
        });
        b.iter(|| black_box(pl.analyze()))
    });
}

criterion_group!(
    benches,
    bench_fluid_water,
    bench_fluid_seawater,
    bench_capillary,
    bench_vacuum_lift,
    bench_pipe_flow,
    bench_flow,
    bench_manning,
    bench_pipeline,
);
criterion_main!(benches);
