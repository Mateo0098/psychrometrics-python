"""Experimental psychrometric relations from the original university model.

This module is deliberately separate from :mod:`psychrometrics`. Its
temperature inputs are in kelvin, matching the notation in the source report,
and it must not be treated as a second official API.

Primary local source: ``docs/Copia de Proyecto A.docx``. The report transcribes
ASAE D271.2 APR1979 (R2014) equations and Monteith & Unsworth (2013), equation
2.2.10. Only relations that are complete and dimensionally reproducible in the
report are implemented here.
"""

from __future__ import annotations

import math


ASAE_MIN_TEMPERATURE_K = 255.38
ASAE_PHASE_TRANSITION_K = 273.16
ASAE_LINEAR_VAPORIZATION_MAX_K = 338.72
ASAE_MAX_TEMPERATURE_K = 533.16

_ASAE_PRESSURE_SCALE_PA = 22_105_649.25
_ASAE_WET_BULB_MASS_RATIO = 0.62194
_ASAE_DRY_AIR_HEAT_CAPACITY_J_KG_K = 1006.9254
_WET_BULB_TOLERANCE_K = 0.001
_WET_BULB_MAX_ITERATIONS = 100


def _finite(name: str, value: float) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a real number, not bool")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _in_range(name: str, value: float, minimum: float, maximum: float) -> float:
    value = _finite(name, value)
    if not minimum <= value <= maximum:
        raise ValueError(f"{name} must be between {minimum} and {maximum}")
    return value


def _positive_pressure(name: str, value: float) -> float:
    pressure = _finite(name, value)
    if pressure <= 0.0:
        raise ValueError(f"{name} must be greater than 0 Pa")
    return pressure


def asae_saturation_vapor_pressure(temperature_k: float) -> float:
    """Return ASAE saturation vapor pressure in Pa for temperature in K.

    The source report gives an ice correlation for 255.38--273.16 K and a
    liquid-water correlation for 273.16--533.16 K. At their shared endpoint,
    this implementation selects the ice relation, mirroring the report's
    piecewise ordering. ``ln`` is interpreted as the natural logarithm.
    """
    temperature_k = _in_range(
        "temperature_k",
        temperature_k,
        ASAE_MIN_TEMPERATURE_K,
        ASAE_MAX_TEMPERATURE_K,
    )
    if temperature_k <= ASAE_PHASE_TRANSITION_K:
        log_pressure = (
            31.9602
            - 6270.3605 / temperature_k
            - 0.46057 * math.log(temperature_k)
        )
        return math.exp(log_pressure)

    a = -27_405.526
    b = 97.5413
    c = -0.146244
    d = 0.12558e-3
    e = -0.48502e-7
    f = 4.34903
    g = 0.39381e-2
    numerator = (
        a
        + b * temperature_k
        + c * temperature_k**2
        + d * temperature_k**3
        + e * temperature_k**4
    )
    denominator = f * temperature_k - g * temperature_k**2
    return _ASAE_PRESSURE_SCALE_PA * math.exp(numerator / denominator)


def asae_saturation_temperature(vapor_pressure_pa: float) -> float:
    """Return ASAE saturation temperature in K from pressure in Pa.

    Implements the report's eighth-degree inverse polynomial. Its documented
    pressure range is 620.52--4,688,396 Pa, narrower than the forward
    correlation. No extrapolation is allowed.
    """
    vapor_pressure_pa = _in_range(
        "vapor_pressure_pa", vapor_pressure_pa, 620.52, 4_688_396.0
    )
    coefficients = (
        19.5322,
        13.6626,
        1.17678,
        -0.189693,
        0.087453,
        -0.0174053,
        0.00214768,
        -0.138343e-3,
        0.38e-5,
    )
    logarithm = math.log(0.00145 * vapor_pressure_pa)
    return 255.38 + sum(
        coefficient * logarithm**index
        for index, coefficient in enumerate(coefficients)
    )


def monteith_vapor_pressure(
    temperature_k: float,
    reference_temperature_k: float,
    reference_vapor_pressure_pa: float,
    coefficient: float,
) -> float:
    """Evaluate the Monteith vapor-pressure relation in Pa.

    The documented equation is ``e(T) = e(T*) exp[A (T - T*) / T]``.
    Temperatures are K, both pressures use the same units (Pa here), and ``A``
    is dimensionless. The report states A=19.65 for reference temperatures
    273--293 K and A=18 for 293--313 K, but does not supply reference
    pressures; callers must therefore provide the anchor explicitly.
    """
    temperature_k = _finite("temperature_k", temperature_k)
    reference_temperature_k = _finite(
        "reference_temperature_k", reference_temperature_k
    )
    reference_vapor_pressure_pa = _positive_pressure(
        "reference_vapor_pressure_pa", reference_vapor_pressure_pa
    )
    coefficient = _finite("coefficient", coefficient)
    if temperature_k <= 0.0 or reference_temperature_k <= 0.0:
        raise ValueError("absolute temperatures must be greater than 0 K")
    if coefficient <= 0.0:
        raise ValueError("coefficient must be greater than 0")
    return reference_vapor_pressure_pa * math.exp(
        coefficient
        * (temperature_k - reference_temperature_k)
        / temperature_k
    )


def monteith_temperature(
    vapor_pressure_pa: float,
    reference_temperature_k: float,
    reference_vapor_pressure_pa: float,
    coefficient: float,
) -> float:
    """Algebraically invert the documented Monteith relation, returning K.

    This is an exact inversion of :func:`monteith_vapor_pressure`, not an
    additional empirical correlation. A singular or nonphysical result is
    rejected.
    """
    vapor_pressure_pa = _positive_pressure("vapor_pressure_pa", vapor_pressure_pa)
    reference_temperature_k = _finite(
        "reference_temperature_k", reference_temperature_k
    )
    reference_vapor_pressure_pa = _positive_pressure(
        "reference_vapor_pressure_pa", reference_vapor_pressure_pa
    )
    coefficient = _finite("coefficient", coefficient)
    if reference_temperature_k <= 0.0:
        raise ValueError("reference_temperature_k must be greater than 0 K")
    if coefficient <= 0.0:
        raise ValueError("coefficient must be greater than 0")

    denominator = 1.0 - math.log(
        vapor_pressure_pa / reference_vapor_pressure_pa
    ) / coefficient
    if denominator <= 0.0:
        raise ValueError("inputs imply a singular or nonphysical temperature")
    return reference_temperature_k / denominator


def asae_latent_heat_sublimation(temperature_k: float) -> float:
    """Return ASAE saturated latent heat of sublimation in J/kg.

    Documented range: 255.38--273.16 K.
    """
    temperature_k = _in_range(
        "temperature_k",
        temperature_k,
        ASAE_MIN_TEMPERATURE_K,
        ASAE_PHASE_TRANSITION_K,
    )
    return 2_839_683.144 - 212.56384 * (temperature_k - 255.38)


def asae_latent_heat_vaporization(temperature_k: float) -> float:
    """Return the reproducible ASAE latent heat of vaporization in J/kg.

    Only the report's complete linear branch, 273.16--338.72 K, is
    implemented. The reported 338.72--533.16 K expression is dimensionally
    incomplete (an apparent radical or outer exponent is absent), so this
    function deliberately refuses that range rather than guessing.
    """
    temperature_k = _in_range(
        "temperature_k",
        temperature_k,
        ASAE_PHASE_TRANSITION_K,
        ASAE_LINEAR_VAPORIZATION_MAX_K,
    )
    return 2_502_535.259 - 2_385.76424 * (
        temperature_k - ASAE_PHASE_TRANSITION_K
    )


def _asae_phase_change_latent_heat(temperature_k: float) -> float:
    if temperature_k <= ASAE_PHASE_TRANSITION_K:
        return asae_latent_heat_sublimation(temperature_k)
    return asae_latent_heat_vaporization(temperature_k)


def asae_wet_bulb_line_residual(
    dry_bulb_temperature_k: float,
    wet_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> float:
    """Return the ASAE wet-bulb-line residual in Pa.

    Zero satisfies report equation 5::

        P_swb - P_v = B' (T_wb - T)

    where ``B'`` is evaluated exactly as transcribed. Temperatures are K,
    pressures are Pa, and latent heat is J/kg. The usable wet-bulb range is
    255.38--338.72 K because the higher-temperature latent-heat branch is
    ambiguous in the report.
    """
    dry_bulb_temperature_k = _in_range(
        "dry_bulb_temperature_k",
        dry_bulb_temperature_k,
        ASAE_MIN_TEMPERATURE_K,
        ASAE_LINEAR_VAPORIZATION_MAX_K,
    )
    wet_bulb_temperature_k = _in_range(
        "wet_bulb_temperature_k",
        wet_bulb_temperature_k,
        ASAE_MIN_TEMPERATURE_K,
        ASAE_LINEAR_VAPORIZATION_MAX_K,
    )
    vapor_pressure_pa = _finite("vapor_pressure_pa", vapor_pressure_pa)
    atmospheric_pressure_pa = _positive_pressure(
        "atmospheric_pressure_pa", atmospheric_pressure_pa
    )
    if not 0.0 <= vapor_pressure_pa < atmospheric_pressure_pa:
        raise ValueError(
            "vapor_pressure_pa must be between 0 and atmospheric pressure"
        )

    wet_bulb_saturation_pressure = asae_saturation_vapor_pressure(
        wet_bulb_temperature_k
    )
    latent_heat = _asae_phase_change_latent_heat(wet_bulb_temperature_k)
    coefficient = (
        _ASAE_DRY_AIR_HEAT_CAPACITY_J_KG_K
        * (wet_bulb_saturation_pressure - atmospheric_pressure_pa)
        * (1.0 + 0.15577 * vapor_pressure_pa / atmospheric_pressure_pa)
        / (_ASAE_WET_BULB_MASS_RATIO * latent_heat)
    )
    return (
        wet_bulb_saturation_pressure
        - vapor_pressure_pa
        - coefficient
        * (wet_bulb_temperature_k - dry_bulb_temperature_k)
    )


def experimental_asae_wet_bulb_temperature(
    dry_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> float:
    """Solve the reconstructed ASAE wet-bulb line by bisection, returning K.

    This experimental solver is restricted to 255.38--338.72 K. It does not
    use the incomplete Dossat proposal, and it is not part of the stable
    ASHRAE state API.
    """
    dry_bulb_temperature_k = _in_range(
        "dry_bulb_temperature_k",
        dry_bulb_temperature_k,
        ASAE_MIN_TEMPERATURE_K,
        ASAE_LINEAR_VAPORIZATION_MAX_K,
    )
    atmospheric_pressure_pa = _positive_pressure(
        "atmospheric_pressure_pa", atmospheric_pressure_pa
    )
    vapor_pressure_pa = _finite("vapor_pressure_pa", vapor_pressure_pa)
    dry_bulb_saturation_pressure = asae_saturation_vapor_pressure(
        dry_bulb_temperature_k
    )
    if not 0.0 < vapor_pressure_pa <= dry_bulb_saturation_pressure:
        raise ValueError(
            "vapor_pressure_pa must be positive and no greater than saturation "
            "pressure at the dry-bulb temperature"
        )
    if vapor_pressure_pa >= atmospheric_pressure_pa:
        raise ValueError("vapor_pressure_pa must be less than atmospheric pressure")

    low = ASAE_MIN_TEMPERATURE_K
    high = dry_bulb_temperature_k
    low_residual = asae_wet_bulb_line_residual(
        dry_bulb_temperature_k, low, vapor_pressure_pa, atmospheric_pressure_pa
    )
    high_residual = asae_wet_bulb_line_residual(
        dry_bulb_temperature_k, high, vapor_pressure_pa, atmospheric_pressure_pa
    )
    if low_residual == 0.0:
        return low
    if high_residual == 0.0:
        return high
    if low_residual * high_residual > 0.0:
        raise ValueError("wet-bulb solution is outside the documented ASAE range")

    for _ in range(_WET_BULB_MAX_ITERATIONS):
        midpoint = (low + high) / 2.0
        midpoint_residual = asae_wet_bulb_line_residual(
            dry_bulb_temperature_k,
            midpoint,
            vapor_pressure_pa,
            atmospheric_pressure_pa,
        )
        if midpoint_residual == 0.0:
            return midpoint
        if low_residual * midpoint_residual <= 0.0:
            high = midpoint
        else:
            low = midpoint
            low_residual = midpoint_residual
        if high - low <= _WET_BULB_TOLERANCE_K:
            return (low + high) / 2.0

    raise RuntimeError(
        "experimental ASAE wet-bulb bisection did not converge within "
        f"{_WET_BULB_MAX_ITERATIONS} iterations"
    )
