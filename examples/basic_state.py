"""Minimal executable example for the core psychrometric calculations."""

from pathlib import Path
import sys


# Allow direct execution with: python examples/basic_state.py
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from psychrometrics import calculate_state  # noqa: E402


def main() -> None:
    state = calculate_state(
        dry_bulb_temperature_c=25.0,
        relative_humidity=0.50,
        pressure_pa=101_325.0,
    )

    print("Basic psychrometric state (SI)")
    print(f"Dry-bulb temperature: {state.dry_bulb_temperature_c:.2f} °C")
    print(f"Relative humidity: {state.relative_humidity:.1%}")
    print(f"Pressure: {state.pressure_pa:.0f} Pa")
    print(f"Humidity ratio: {state.humidity_ratio_kg_kg_dry_air:.6f} kg/kg dry air")
    print(f"Enthalpy: {state.enthalpy_kj_kg_dry_air:.2f} kJ/kg dry air")
    print(f"Specific volume: {state.specific_volume_m3_kg_dry_air:.4f} m³/kg dry air")
    print(f"Dew-point temperature: {state.dew_point_temperature_c:.2f} °C")


if __name__ == "__main__":
    main()
