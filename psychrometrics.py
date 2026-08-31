"""Core psychrometric calculations in SI units.

The equations implemented here follow Chapter 1, *Psychrometrics*, of the
ASHRAE Handbook—Fundamentals (2021).  The saturation-pressure coefficients
and equation numbering are also available in the open-source PsychroLib
implementation (ASHRAE, SI unit system):
https://github.com/psychrometrics/psychrolib

Temperatures are dry-bulb temperatures in degrees Celsius, pressures are in
pascal, and properties expressed per mass of air use kilograms of dry air.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


STANDARD_ATMOSPHERIC_PRESSURE_PA = 101_325.0
"""Standard atmospheric pressure used only as a convenient public default."""

_TRIPLE_POINT_WATER_C = 0.01
_MIN_TEMPERATURE_C = -100.0
_MAX_TEMPERATURE_C = 200.0
_MOLECULAR_WEIGHT_RATIO = 0.621945
_DRY_AIR_GAS_CONSTANT_J_KG_K = 287.042


def _finite(name: str, value: float) -> float:
    """Return *value* as float, rejecting booleans and non-finite values."""
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a real number, not bool")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _temperature_c(value: float) -> float:
    temperature = _finite("temperature_c", value)
    if not _MIN_TEMPERATURE_C <= temperature <= _MAX_TEMPERATURE_C:
        raise ValueError("temperature_c must be between -100 and 200 °C")
    return temperature


def _pressure_pa(value: float) -> float:
    pressure = _finite("pressure_pa", value)
    if pressure <= 0.0:
        raise ValueError("pressure_pa must be greater than 0 Pa")
    return pressure


def saturation_vapor_pressure(temperature_c: float) -> float:
    """Return saturation vapor pressure in Pa at ``temperature_c``.

    Uses ASHRAE Handbook—Fundamentals (2021), Chapter 1, equations 5 and 6.
    The ice equation is used from -100 °C through the water triple point
    (0.01 °C); the liquid-water equation is used above it through 200 °C.
    Temperature is converted to kelvin before evaluating the correlation.
    """
    temperature_c = _temperature_c(temperature_c)
    temperature_k = temperature_c + 273.15

    if temperature_c <= _TRIPLE_POINT_WATER_C:
        log_pressure = (
            -5.6745359e3 / temperature_k
            + 6.3925247
            - 9.677843e-3 * temperature_k
            + 6.2215701e-7 * temperature_k**2
            + 2.0747825e-9 * temperature_k**3
            - 9.484024e-13 * temperature_k**4
            + 4.1635019 * math.log(temperature_k)
        )
    else:
        log_pressure = (
            -5.8002206e3 / temperature_k
            + 1.3914993
            - 4.8640239e-2 * temperature_k
            + 4.1764768e-5 * temperature_k**2
            - 1.4452093e-8 * temperature_k**3
            + 6.5459673 * math.log(temperature_k)
        )
    return math.exp(log_pressure)


def vapor_pressure(
    temperature_c: float,
    relative_humidity: float,
) -> float:
    """Return water-vapor partial pressure in Pa.

    ``relative_humidity`` is a fraction from 0 to 1, not a percentage.
    This evaluates ``p_w = phi * p_ws`` (ASHRAE 2021, Chapter 1, eq. 12).
    """
    temperature_c = _temperature_c(temperature_c)
    relative_humidity = _finite("relative_humidity", relative_humidity)
    if not 0.0 <= relative_humidity <= 1.0:
        raise ValueError("relative_humidity must be between 0 and 1")
    return relative_humidity * saturation_vapor_pressure(temperature_c)


def humidity_ratio(vapor_pressure_pa: float, pressure_pa: float) -> float:
    """Return humidity ratio in kg water/kg dry air.

    Uses ``W = 0.621945 p_w / (p - p_w)`` from ASHRAE 2021, Chapter 1,
    eq. 20.  Both input pressures are in Pa.
    """
    vapor_pressure_pa = _finite("vapor_pressure_pa", vapor_pressure_pa)
    pressure_pa = _pressure_pa(pressure_pa)
    if vapor_pressure_pa < 0.0:
        raise ValueError("vapor_pressure_pa must be at least 0 Pa")
    if vapor_pressure_pa >= pressure_pa:
        raise ValueError("vapor_pressure_pa must be less than pressure_pa")
    return _MOLECULAR_WEIGHT_RATIO * vapor_pressure_pa / (
        pressure_pa - vapor_pressure_pa
    )


def moist_air_enthalpy(temperature_c: float, humidity_ratio_kg_kg: float) -> float:
    """Return moist-air enthalpy in kJ/kg dry air.

    Uses ``h = 1.006 t + W (2501 + 1.86 t)`` from ASHRAE 2021,
    Chapter 1, eq. 30, with ``t`` in °C and ``W`` in kg/kg dry air.
    """
    temperature_c = _temperature_c(temperature_c)
    humidity_ratio_kg_kg = _finite(
        "humidity_ratio_kg_kg", humidity_ratio_kg_kg
    )
    if humidity_ratio_kg_kg < 0.0:
        raise ValueError("humidity_ratio_kg_kg must be at least 0 kg/kg")
    return 1.006 * temperature_c + humidity_ratio_kg_kg * (
        2501.0 + 1.86 * temperature_c
    )


def specific_volume(
    temperature_c: float,
    humidity_ratio_kg_kg: float,
    pressure_pa: float,
) -> float:
    """Return specific volume in m³/kg dry air.

    Uses the ideal-gas moist-air relation from ASHRAE 2021, Chapter 1,
    eq. 28: ``v = R_da T (1 + 1.607858 W) / p`` in SI units.
    """
    temperature_c = _temperature_c(temperature_c)
    humidity_ratio_kg_kg = _finite(
        "humidity_ratio_kg_kg", humidity_ratio_kg_kg
    )
    pressure_pa = _pressure_pa(pressure_pa)
    if humidity_ratio_kg_kg < 0.0:
        raise ValueError("humidity_ratio_kg_kg must be at least 0 kg/kg")
    temperature_k = temperature_c + 273.15
    return (
        _DRY_AIR_GAS_CONSTANT_J_KG_K
        * temperature_k
        * (1.0 + 1.607858 * humidity_ratio_kg_kg)
        / pressure_pa
    )


def dew_point_temperature(vapor_pressure_pa: float) -> float:
    """Return dew-point temperature in °C from vapor pressure in Pa.

    Numerically inverts the ASHRAE saturation-pressure equations with
    bisection.  The supported pressure range corresponds to -100 to 200 °C.
    A strictly positive pressure is required because dry air has no finite
    dew point in this model.
    """
    vapor_pressure_pa = _finite("vapor_pressure_pa", vapor_pressure_pa)
    minimum_pressure = saturation_vapor_pressure(_MIN_TEMPERATURE_C)
    maximum_pressure = saturation_vapor_pressure(_MAX_TEMPERATURE_C)
    if not minimum_pressure <= vapor_pressure_pa <= maximum_pressure:
        raise ValueError(
            "vapor_pressure_pa is outside the saturation-pressure range "
            "for -100 to 200 °C"
        )

    low = _MIN_TEMPERATURE_C
    high = _MAX_TEMPERATURE_C
    for _ in range(80):
        midpoint = (low + high) / 2.0
        if saturation_vapor_pressure(midpoint) < vapor_pressure_pa:
            low = midpoint
        else:
            high = midpoint
    return (low + high) / 2.0


@dataclass(frozen=True)
class PsychrometricState:
    """Basic moist-air state, with all quantities in the documented SI units."""

    dry_bulb_temperature_c: float
    relative_humidity: float
    pressure_pa: float
    saturation_vapor_pressure_pa: float
    vapor_pressure_pa: float
    humidity_ratio_kg_kg_dry_air: float
    enthalpy_kj_kg_dry_air: float
    specific_volume_m3_kg_dry_air: float
    dew_point_temperature_c: float | None


def calculate_state(
    dry_bulb_temperature_c: float,
    relative_humidity: float,
    pressure_pa: float = STANDARD_ATMOSPHERIC_PRESSURE_PA,
) -> PsychrometricState:
    """Calculate a basic state from dry-bulb temperature, RH, and pressure.

    Relative humidity is supplied as a fraction.  At zero relative humidity,
    ``dew_point_temperature_c`` is ``None`` because the model has no finite
    dew point for zero vapor pressure.
    """
    dry_bulb_temperature_c = _temperature_c(dry_bulb_temperature_c)
    pressure_pa = _pressure_pa(pressure_pa)
    saturation_pressure = saturation_vapor_pressure(dry_bulb_temperature_c)
    partial_pressure = vapor_pressure(
        dry_bulb_temperature_c, relative_humidity
    )
    ratio = humidity_ratio(partial_pressure, pressure_pa)
    dew_point = (
        dew_point_temperature(partial_pressure)
        if partial_pressure > 0.0
        else None
    )
    return PsychrometricState(
        dry_bulb_temperature_c=dry_bulb_temperature_c,
        relative_humidity=float(relative_humidity),
        pressure_pa=pressure_pa,
        saturation_vapor_pressure_pa=saturation_pressure,
        vapor_pressure_pa=partial_pressure,
        humidity_ratio_kg_kg_dry_air=ratio,
        enthalpy_kj_kg_dry_air=moist_air_enthalpy(
            dry_bulb_temperature_c, ratio
        ),
        specific_volume_m3_kg_dry_air=specific_volume(
            dry_bulb_temperature_c, ratio, pressure_pa
        ),
        dew_point_temperature_c=dew_point,
    )
