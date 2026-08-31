"""Reconstruct one psychrometric state from several supported input pairs."""

from pathlib import Path
import sys


# Allow direct execution with: python examples/flexible_state_solver.py
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from psychrometrics import solve_state  # noqa: E402


def main() -> None:
    pressure_pa = 101_325.0
    states = {
        "Tdb + RH": solve_state(
            25.0,
            pressure_pa=pressure_pa,
            relative_humidity=0.50,
        ),
        "Tdb + Twb": solve_state(
            25.0,
            pressure_pa=pressure_pa,
            wet_bulb_temperature_c=17.8894,
        ),
        "h + W": solve_state(
            pressure_pa=pressure_pa,
            enthalpy_kj_kg_dry_air=50.32196,
            humidity_ratio_kg_kg_dry_air=0.00988115,
        ),
        "Tdp + RH": solve_state(
            pressure_pa=pressure_pa,
            dew_point_temperature_c=13.86397,
            relative_humidity=0.50,
        ),
    }

    print("Flexible psychrometric state solver (SI)")
    print("Input pair     Tdb (°C)       RH      W (kg/kg)   h (kJ/kg)")
    for input_pair, state in states.items():
        print(
            f"{input_pair:<12} "
            f"{state.dry_bulb_temperature_c:9.4f} "
            f"{state.relative_humidity:8.5f} "
            f"{state.humidity_ratio_kg_kg_dry_air:14.8f} "
            f"{state.enthalpy_kj_kg_dry_air:11.4f}"
        )


if __name__ == "__main__":
    main()
