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
_DEW_POINT_TOLERANCE_C = 0.001
_DEW_POINT_MAX_ITERATIONS = 100
_STATE_TEMPERATURE_TOLERANCE_C = 0.001
_STATE_ROOT_TOLERANCE_C = 0.001
_STATE_ROOT_MAX_ITERATIONS = 100


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
    The ASHRAE correlation and this API both support -100 to 200 °C. The ice
    equation is used through the water triple point (0.01 °C), inclusive; the
    liquid-water equation is used above it. Temperature is converted to kelvin
    before evaluating the correlation.
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
    The API restricts temperature to the -100 to 200 °C range used throughout
    this module. This equation does not introduce a separate ice-phase term.
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
    The API accepts temperatures from -100 to 200 °C and pressure in Pa.
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

    Numerically inverts the ASHRAE ice/liquid saturation-pressure correlation
    with bisection. The pressure range corresponds to the correlation and API
    limits of -100 to 200 °C. Iteration stops when the enclosing temperature
    interval is no wider than 0.001 °C, with at most 100 iterations. A pressure
    below the supported range (including zero) has no representable dew point
    in this API.
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
    for _ in range(_DEW_POINT_MAX_ITERATIONS):
        midpoint = (low + high) / 2.0
        if saturation_vapor_pressure(midpoint) < vapor_pressure_pa:
            low = midpoint
        else:
            high = midpoint
        if high - low <= _DEW_POINT_TOLERANCE_C:
            return (low + high) / 2.0

    raise RuntimeError(
        "dew-point bisection did not converge within "
        f"{_DEW_POINT_MAX_ITERATIONS} iterations"
    )


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


def _vapor_pressure_from_humidity_ratio(
    humidity_ratio_kg_kg_dry_air: float,
    pressure_pa: float,
) -> float:
    """Invert ASHRAE equation 20 and return vapor pressure in Pa."""
    humidity_ratio_kg_kg_dry_air = _finite(
        "humidity_ratio_kg_kg_dry_air", humidity_ratio_kg_kg_dry_air
    )
    pressure_pa = _pressure_pa(pressure_pa)
    if humidity_ratio_kg_kg_dry_air < 0.0:
        raise ValueError(
            "humidity_ratio_kg_kg_dry_air must be at least 0 kg/kg"
        )
    return (
        pressure_pa
        * humidity_ratio_kg_kg_dry_air
        / (_MOLECULAR_WEIGHT_RATIO + humidity_ratio_kg_kg_dry_air)
    )


def _humidity_ratio_from_wet_bulb(
    dry_bulb_temperature_c: float,
    wet_bulb_temperature_c: float,
    pressure_pa: float,
) -> float:
    """Return humidity ratio from dry- and wet-bulb temperatures.

    Uses the SI forms of ASHRAE Handbook—Fundamentals (2017), Chapter 1,
    equations 33 and 35, as implemented by PsychroLib 2.5.0. The liquid-water
    branch is used at and above 0 °C wet bulb; the ice branch is used below
    0 °C. Temperatures are in °C, pressure is in Pa, and the result is kg water
    per kg dry air.
    """
    dry_bulb_temperature_c = _temperature_c(dry_bulb_temperature_c)
    wet_bulb_temperature_c = _temperature_c(wet_bulb_temperature_c)
    pressure_pa = _pressure_pa(pressure_pa)
    if (
        wet_bulb_temperature_c
        > dry_bulb_temperature_c + _STATE_TEMPERATURE_TOLERANCE_C
    ):
        raise ValueError(
            "wet_bulb_temperature_c must not exceed dry_bulb_temperature_c"
        )
    wet_bulb_temperature_c = min(
        wet_bulb_temperature_c, dry_bulb_temperature_c
    )

    saturation_pressure = saturation_vapor_pressure(wet_bulb_temperature_c)
    saturated_ratio = humidity_ratio(saturation_pressure, pressure_pa)
    if wet_bulb_temperature_c >= 0.0:
        result = (
            (2501.0 - 2.326 * wet_bulb_temperature_c) * saturated_ratio
            - 1.006 * (dry_bulb_temperature_c - wet_bulb_temperature_c)
        ) / (
            2501.0
            + 1.86 * dry_bulb_temperature_c
            - 4.186 * wet_bulb_temperature_c
        )
    else:
        result = (
            (2830.0 - 0.24 * wet_bulb_temperature_c) * saturated_ratio
            - 1.006 * (dry_bulb_temperature_c - wet_bulb_temperature_c)
        ) / (
            2830.0
            + 1.86 * dry_bulb_temperature_c
            - 2.1 * wet_bulb_temperature_c
        )
    if result < 0.0:
        raise ValueError(
            "dry- and wet-bulb temperatures imply a negative humidity ratio"
        )
    return result


def _relative_humidity_from_vapor_pressure(
    dry_bulb_temperature_c: float,
    vapor_pressure_pa: float,
    pressure_pa: float,
) -> float:
    """Return RH after validating total pressure and saturation limits."""
    dry_bulb_temperature_c = _temperature_c(dry_bulb_temperature_c)
    vapor_pressure_pa = _finite("vapor_pressure_pa", vapor_pressure_pa)
    pressure_pa = _pressure_pa(pressure_pa)
    if vapor_pressure_pa < 0.0 or vapor_pressure_pa >= pressure_pa:
        raise ValueError("vapor_pressure_pa must be in the range [0, pressure_pa)")
    relative_humidity = vapor_pressure_pa / saturation_vapor_pressure(
        dry_bulb_temperature_c
    )
    if not 0.0 <= relative_humidity <= 1.0:
        raise ValueError(
            "inputs imply a physically impossible supersaturated state"
        )
    return relative_humidity


def _bisect_monotonic_temperature(
    value_at_temperature,
    target: float,
    low_c: float,
    high_c: float,
    property_name: str,
) -> float:
    """Invert a continuous increasing property over a bounded °C interval."""
    low_value = value_at_temperature(low_c)
    high_value = value_at_temperature(high_c)
    if target < low_value or target > high_value:
        raise ValueError(
            f"{property_name} has no solution in the supported physical range"
        )
    if target == low_value:
        return low_c
    if target == high_value:
        return high_c

    for _ in range(_STATE_ROOT_MAX_ITERATIONS):
        midpoint = (low_c + high_c) / 2.0
        if value_at_temperature(midpoint) < target:
            low_c = midpoint
        else:
            high_c = midpoint
        if high_c - low_c <= _STATE_ROOT_TOLERANCE_C:
            return (low_c + high_c) / 2.0
    raise RuntimeError(
        f"{property_name} temperature solve did not converge within "
        f"{_STATE_ROOT_MAX_ITERATIONS} iterations"
    )


def _physical_temperature_upper_bound(
    relative_humidity: float,
    pressure_pa: float,
) -> float:
    """Return the highest usable Tdb where RH * saturation pressure < P."""
    if relative_humidity == 0.0:
        return _MAX_TEMPERATURE_C
    limiting_saturation_pressure = pressure_pa / relative_humidity
    maximum_saturation_pressure = saturation_vapor_pressure(_MAX_TEMPERATURE_C)
    if limiting_saturation_pressure > maximum_saturation_pressure:
        return _MAX_TEMPERATURE_C
    minimum_saturation_pressure = saturation_vapor_pressure(_MIN_TEMPERATURE_C)
    if limiting_saturation_pressure <= minimum_saturation_pressure:
        raise ValueError(
            "relative_humidity and pressure leave no physical temperature "
            "in the supported range"
        )
    phase_boundary = dew_point_temperature(limiting_saturation_pressure)
    upper_bound = phase_boundary - _STATE_ROOT_TOLERANCE_C
    while (
        relative_humidity * saturation_vapor_pressure(upper_bound)
        >= pressure_pa
    ):
        upper_bound -= _STATE_ROOT_TOLERANCE_C
    return max(upper_bound, _MIN_TEMPERATURE_C)


def _resolve_tdb_rh(
    values: dict[str, float], pressure_pa: float
) -> tuple[float, float]:
    return (
        _temperature_c(values["dry_bulb_temperature_c"]),
        _finite("relative_humidity", values["relative_humidity"]),
    )


def _resolve_tdb_tdp(
    values: dict[str, float], pressure_pa: float
) -> tuple[float, float]:
    dry_bulb = _temperature_c(values["dry_bulb_temperature_c"])
    dew_point = _temperature_c(values["dew_point_temperature_c"])
    if dew_point > dry_bulb + _STATE_TEMPERATURE_TOLERANCE_C:
        raise ValueError(
            "dew_point_temperature_c must not exceed dry_bulb_temperature_c"
        )
    partial_pressure = saturation_vapor_pressure(min(dew_point, dry_bulb))
    return dry_bulb, _relative_humidity_from_vapor_pressure(
        dry_bulb, partial_pressure, pressure_pa
    )


def _resolve_tdb_w(
    values: dict[str, float], pressure_pa: float
) -> tuple[float, float]:
    dry_bulb = _temperature_c(values["dry_bulb_temperature_c"])
    partial_pressure = _vapor_pressure_from_humidity_ratio(
        values["humidity_ratio_kg_kg_dry_air"], pressure_pa
    )
    return dry_bulb, _relative_humidity_from_vapor_pressure(
        dry_bulb, partial_pressure, pressure_pa
    )


def _resolve_tdb_twb(
    values: dict[str, float], pressure_pa: float
) -> tuple[float, float]:
    dry_bulb = _temperature_c(values["dry_bulb_temperature_c"])
    ratio = _humidity_ratio_from_wet_bulb(
        dry_bulb, values["wet_bulb_temperature_c"], pressure_pa
    )
    partial_pressure = _vapor_pressure_from_humidity_ratio(ratio, pressure_pa)
    return dry_bulb, _relative_humidity_from_vapor_pressure(
        dry_bulb, partial_pressure, pressure_pa
    )


def _resolve_h_w(
    values: dict[str, float], pressure_pa: float
) -> tuple[float, float]:
    enthalpy = _finite(
        "enthalpy_kj_kg_dry_air", values["enthalpy_kj_kg_dry_air"]
    )
    ratio = _finite(
        "humidity_ratio_kg_kg_dry_air",
        values["humidity_ratio_kg_kg_dry_air"],
    )
    if ratio < 0.0:
        raise ValueError("humidity_ratio_kg_kg_dry_air must be at least 0 kg/kg")
    dry_bulb = _temperature_c(
        (enthalpy - 2501.0 * ratio) / (1.006 + 1.86 * ratio)
    )
    partial_pressure = _vapor_pressure_from_humidity_ratio(ratio, pressure_pa)
    return dry_bulb, _relative_humidity_from_vapor_pressure(
        dry_bulb, partial_pressure, pressure_pa
    )


def _resolve_v_w(
    values: dict[str, float], pressure_pa: float
) -> tuple[float, float]:
    volume = _finite(
        "specific_volume_m3_kg_dry_air",
        values["specific_volume_m3_kg_dry_air"],
    )
    ratio = _finite(
        "humidity_ratio_kg_kg_dry_air",
        values["humidity_ratio_kg_kg_dry_air"],
    )
    if volume <= 0.0:
        raise ValueError("specific_volume_m3_kg_dry_air must be greater than 0")
    if ratio < 0.0:
        raise ValueError("humidity_ratio_kg_kg_dry_air must be at least 0 kg/kg")
    temperature_k = (
        volume
        * pressure_pa
        / (_DRY_AIR_GAS_CONSTANT_J_KG_K * (1.0 + 1.607858 * ratio))
    )
    dry_bulb = _temperature_c(temperature_k - 273.15)
    partial_pressure = _vapor_pressure_from_humidity_ratio(ratio, pressure_pa)
    return dry_bulb, _relative_humidity_from_vapor_pressure(
        dry_bulb, partial_pressure, pressure_pa
    )


def _resolve_tdp_rh(
    values: dict[str, float], pressure_pa: float
) -> tuple[float, float]:
    dew_point = _temperature_c(values["dew_point_temperature_c"])
    relative_humidity = _finite("relative_humidity", values["relative_humidity"])
    if not 0.0 < relative_humidity <= 1.0:
        raise ValueError(
            "relative_humidity must be greater than 0 and at most 1 for "
            "the dew-point pair"
        )
    partial_pressure = saturation_vapor_pressure(dew_point)
    if partial_pressure >= pressure_pa:
        raise ValueError("dew point implies vapor_pressure_pa >= pressure_pa")
    target_saturation_pressure = partial_pressure / relative_humidity
    dry_bulb = dew_point_temperature(target_saturation_pressure)
    if dry_bulb < dew_point - _STATE_TEMPERATURE_TOLERANCE_C:
        raise ValueError("dew-point and relative humidity imply Tdb < Tdp")
    return max(dry_bulb, dew_point), relative_humidity


def _resolve_h_rh(
    values: dict[str, float], pressure_pa: float
) -> tuple[float, float]:
    enthalpy = _finite(
        "enthalpy_kj_kg_dry_air", values["enthalpy_kj_kg_dry_air"]
    )
    relative_humidity = _finite("relative_humidity", values["relative_humidity"])
    if not 0.0 <= relative_humidity <= 1.0:
        raise ValueError("relative_humidity must be between 0 and 1")
    upper_bound = _physical_temperature_upper_bound(
        relative_humidity, pressure_pa
    )

    def enthalpy_at_temperature(temperature_c: float) -> float:
        partial_pressure = (
            relative_humidity * saturation_vapor_pressure(temperature_c)
        )
        ratio = humidity_ratio(partial_pressure, pressure_pa)
        return moist_air_enthalpy(temperature_c, ratio)

    dry_bulb = _bisect_monotonic_temperature(
        enthalpy_at_temperature,
        enthalpy,
        _MIN_TEMPERATURE_C,
        upper_bound,
        "enthalpy_kj_kg_dry_air",
    )
    return dry_bulb, relative_humidity


_STATE_PAIR_SOLVERS = {
    frozenset(("dry_bulb_temperature_c", "relative_humidity")): _resolve_tdb_rh,
    frozenset(
        ("dry_bulb_temperature_c", "dew_point_temperature_c")
    ): _resolve_tdb_tdp,
    frozenset(
        ("dry_bulb_temperature_c", "humidity_ratio_kg_kg_dry_air")
    ): _resolve_tdb_w,
    frozenset(
        ("dry_bulb_temperature_c", "wet_bulb_temperature_c")
    ): _resolve_tdb_twb,
    frozenset(
        ("enthalpy_kj_kg_dry_air", "humidity_ratio_kg_kg_dry_air")
    ): _resolve_h_w,
    frozenset(
        ("specific_volume_m3_kg_dry_air", "humidity_ratio_kg_kg_dry_air")
    ): _resolve_v_w,
    frozenset(("dew_point_temperature_c", "relative_humidity")): _resolve_tdp_rh,
    frozenset(("enthalpy_kj_kg_dry_air", "relative_humidity")): _resolve_h_rh,
}


def solve_state(
    dry_bulb_temperature_c: float | None = None,
    *,
    pressure_pa: float,
    relative_humidity: float | None = None,
    dew_point_temperature_c: float | None = None,
    humidity_ratio_kg_kg_dry_air: float | None = None,
    wet_bulb_temperature_c: float | None = None,
    enthalpy_kj_kg_dry_air: float | None = None,
    specific_volume_m3_kg_dry_air: float | None = None,
) -> PsychrometricState:
    """Solve a state from exactly one documented pair and explicit pressure.

    Supported pairs are Tdb+RH, Tdb+Tdp, Tdb+W, Tdb+Twb, h+W, v+W,
    Tdp+RH, and h+RH. Algebraic and bounded numerical inversions use the same
    ASHRAE equations as the stable core. Every route finishes through
    :func:`calculate_state` and returns :class:`PsychrometricState`.
    """
    pressure_pa = _pressure_pa(pressure_pa)
    inputs = {
        "dry_bulb_temperature_c": dry_bulb_temperature_c,
        "relative_humidity": relative_humidity,
        "dew_point_temperature_c": dew_point_temperature_c,
        "humidity_ratio_kg_kg_dry_air": humidity_ratio_kg_kg_dry_air,
        "wet_bulb_temperature_c": wet_bulb_temperature_c,
        "enthalpy_kj_kg_dry_air": enthalpy_kj_kg_dry_air,
        "specific_volume_m3_kg_dry_air": specific_volume_m3_kg_dry_air,
    }
    provided = {name: value for name, value in inputs.items() if value is not None}
    if len(provided) != 2:
        raise ValueError(
            "provide exactly two psychrometric properties from a supported pair"
        )
    resolver = _STATE_PAIR_SOLVERS.get(frozenset(provided))
    if resolver is None:
        raise ValueError(
            "unsupported psychrometric property pair: "
            + " + ".join(sorted(provided))
        )
    resolved_dry_bulb, resolved_relative_humidity = resolver(
        provided, pressure_pa
    )
    return calculate_state(
        resolved_dry_bulb,
        resolved_relative_humidity,
        pressure_pa,
    )
