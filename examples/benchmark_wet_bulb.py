"""Reproducible benchmark of baseline and Monteith-assisted ASAE wet bulb."""

from argparse import ArgumentParser
from pathlib import Path
import statistics
import sys
from time import perf_counter


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from alternative_methods import (  # noqa: E402
    diagnose_experimental_asae_wet_bulb,
    diagnose_wet_bulb_asae_monteith_assisted,
)
from psychrometrics import vapor_pressure  # noqa: E402


TEMPERATURES_C = (-5.0, 0.0, 5.0, 15.0, 25.0, 35.0, 40.0)
RELATIVE_HUMIDITIES = (0.20, 0.40, 0.60, 0.80, 0.95)
PRESSURES_PA = (101_325.0, 80_000.0)


def build_cases() -> list[tuple[float, float, float]]:
    """Return 70 valid ``(dry_bulb_K, vapor_pressure_Pa, pressure_Pa)`` cases."""
    return [
        (
            temperature_c + 273.15,
            vapor_pressure(temperature_c, relative_humidity),
            pressure_pa,
        )
        for pressure_pa in PRESSURES_PA
        for temperature_c in TEMPERATURES_C
        for relative_humidity in RELATIVE_HUMIDITIES
    ]


def _time_solver(solver, cases: list[tuple[float, float, float]], repetitions: int) -> float:
    started = perf_counter()
    for _ in range(repetitions):
        for case in cases:
            solver(*case)
    return perf_counter() - started


def benchmark(repetitions: int, rounds: int) -> None:
    cases = build_cases()
    baseline_diagnostics = [
        diagnose_experimental_asae_wet_bulb(*case) for case in cases
    ]
    assisted_diagnostics = [
        diagnose_wet_bulb_asae_monteith_assisted(*case) for case in cases
    ]

    # Warm both paths before measuring them. Timings include validation,
    # diagnostics construction, Monteith estimation, and bracket verification.
    for case in cases:
        diagnose_experimental_asae_wet_bulb(*case)
        diagnose_wet_bulb_asae_monteith_assisted(*case)

    baseline_samples = []
    assisted_samples = []
    for round_index in range(rounds):
        if round_index % 2 == 0:
            baseline_samples.append(
                _time_solver(
                    diagnose_experimental_asae_wet_bulb, cases, repetitions
                )
            )
            assisted_samples.append(
                _time_solver(
                    diagnose_wet_bulb_asae_monteith_assisted, cases, repetitions
                )
            )
        else:
            assisted_samples.append(
                _time_solver(
                    diagnose_wet_bulb_asae_monteith_assisted, cases, repetitions
                )
            )
            baseline_samples.append(
                _time_solver(
                    diagnose_experimental_asae_wet_bulb, cases, repetitions
                )
            )
    baseline_seconds = statistics.median(baseline_samples)
    assisted_seconds = statistics.median(assisted_samples)

    baseline_iterations = [item.iterations for item in baseline_diagnostics]
    assisted_iterations = [item.iterations for item in assisted_diagnostics]
    reductions = [
        100.0 * (baseline - assisted) / baseline
        for baseline, assisted in zip(
            baseline_iterations, assisted_iterations, strict=True
        )
    ]
    fallbacks = sum(item.fallback_used for item in assisted_diagnostics)
    worse_cases = sum(
        assisted > baseline
        for baseline, assisted in zip(
            baseline_iterations, assisted_iterations, strict=True
        )
    )
    maximum_difference = max(
        abs(base.wet_bulb_temperature_k - assisted.wet_bulb_temperature_k)
        for base, assisted in zip(
            baseline_diagnostics, assisted_diagnostics, strict=True
        )
    )

    print("ASAE wet-bulb benchmark")
    print(f"Cases: {len(cases)}")
    print(f"Repetitions per case: {repetitions}")
    print(f"Timing rounds: {rounds} (alternating method order)")
    print(f"Solves per method per round: {len(cases) * repetitions}")
    print(f"Baseline median round time: {baseline_seconds:.6f} s")
    print(f"Monteith-assisted median round time: {assisted_seconds:.6f} s")
    print(
        "Measured time change: "
        f"{100.0 * (baseline_seconds - assisted_seconds) / baseline_seconds:.2f}%"
    )
    print(f"Baseline iterations, mean: {statistics.mean(baseline_iterations):.3f}")
    print(f"Baseline iterations, median: {statistics.median(baseline_iterations):.1f}")
    print(f"Assisted iterations, mean: {statistics.mean(assisted_iterations):.3f}")
    print(f"Assisted iterations, median: {statistics.median(assisted_iterations):.1f}")
    print(f"Mean per-case iteration reduction: {statistics.mean(reductions):.2f}%")
    print(f"Fallback cases: {fallbacks} ({100.0 * fallbacks / len(cases):.2f}%)")
    print(f"Cases with more assisted iterations: {worse_cases}")
    print(f"Maximum baseline/assisted difference: {maximum_difference:.9f} K")


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repetitions",
        type=int,
        default=200,
        help="timed repetitions of the complete 70-case matrix (default: 200)",
    )
    parser.add_argument(
        "--rounds",
        type=int,
        default=5,
        help="timing rounds with alternating method order (default: 5)",
    )
    args = parser.parse_args()
    if args.repetitions <= 0:
        parser.error("--repetitions must be greater than zero")
    if args.rounds <= 0:
        parser.error("--rounds must be greater than zero")
    benchmark(args.repetitions, args.rounds)


if __name__ == "__main__":
    main()
