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

from dataclasses import dataclass
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
_MONTEITH_INITIAL_BRACKET_FRACTION = 0.25
_MONTEITH_MINIMUM_HALF_WIDTH_K = 1.0


@dataclass(frozen=True)
class WetBulbSolveDiagnostics:
    """Diagnostics for an experimental ASAE wet-bulb solve.

    Widths are in K. ``proposed_bracket_width_k`` is the first Monteith-based
    proposal (or the full baseline width), while ``initial_bracket_width_k``
    is the verified bracket actually passed to bisection.
    """

    wet_bulb_temperature_k: float
    iterations: int
    used_monteith_bracket: bool
    fallback_used: bool
    bracket_expansions: int
    proposed_bracket_width_k: float
    initial_bracket_width_k: float
    final_bracket_width_k: float
    monteith_dew_point_estimate_k: float | None = None
    wet_bulb_estimate_k: float | None = None


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


def _validate_wet_bulb_solver_inputs(
    dry_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> tuple[float, float, float]:
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
    return dry_bulb_temperature_k, vapor_pressure_pa, atmospheric_pressure_pa


def _wet_bulb_residual(
    dry_bulb_temperature_k: float,
    wet_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> float:
    return asae_wet_bulb_line_residual(
        dry_bulb_temperature_k,
        wet_bulb_temperature_k,
        vapor_pressure_pa,
        atmospheric_pressure_pa,
    )


def _has_sign_change(low_residual: float, high_residual: float) -> bool:
    return (
        low_residual == 0.0
        or high_residual == 0.0
        or low_residual * high_residual < 0.0
    )


def _bisect_wet_bulb_bracket(
    dry_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
    low: float,
    high: float,
    *,
    used_monteith_bracket: bool,
    fallback_used: bool,
    bracket_expansions: int,
    proposed_bracket_width_k: float,
    monteith_dew_point_estimate_k: float | None = None,
    wet_bulb_estimate_k: float | None = None,
) -> WetBulbSolveDiagnostics:
    initial_bracket_width = high - low
    low_residual = _wet_bulb_residual(
        dry_bulb_temperature_k, low, vapor_pressure_pa, atmospheric_pressure_pa
    )
    high_residual = _wet_bulb_residual(
        dry_bulb_temperature_k, high, vapor_pressure_pa, atmospheric_pressure_pa
    )
    if not _has_sign_change(low_residual, high_residual):
        raise ValueError("wet-bulb bracket does not contain an ASAE root")

    if low_residual == 0.0:
        return WetBulbSolveDiagnostics(
            low,
            0,
            used_monteith_bracket,
            fallback_used,
            bracket_expansions,
            proposed_bracket_width_k,
            initial_bracket_width,
            initial_bracket_width,
            monteith_dew_point_estimate_k,
            wet_bulb_estimate_k,
        )
    if high_residual == 0.0:
        return WetBulbSolveDiagnostics(
            high,
            0,
            used_monteith_bracket,
            fallback_used,
            bracket_expansions,
            proposed_bracket_width_k,
            initial_bracket_width,
            initial_bracket_width,
            monteith_dew_point_estimate_k,
            wet_bulb_estimate_k,
        )

    for iteration in range(1, _WET_BULB_MAX_ITERATIONS + 1):
        midpoint = (low + high) / 2.0
        midpoint_residual = _wet_bulb_residual(
            dry_bulb_temperature_k,
            midpoint,
            vapor_pressure_pa,
            atmospheric_pressure_pa,
        )
        if midpoint_residual == 0.0:
            return WetBulbSolveDiagnostics(
                midpoint,
                iteration,
                used_monteith_bracket,
                fallback_used,
                bracket_expansions,
                proposed_bracket_width_k,
                initial_bracket_width,
                0.0,
                monteith_dew_point_estimate_k,
                wet_bulb_estimate_k,
            )
        if low_residual * midpoint_residual <= 0.0:
            high = midpoint
        else:
            low = midpoint
            low_residual = midpoint_residual
        if high - low <= _WET_BULB_TOLERANCE_K:
            return WetBulbSolveDiagnostics(
                (low + high) / 2.0,
                iteration,
                used_monteith_bracket,
                fallback_used,
                bracket_expansions,
                proposed_bracket_width_k,
                initial_bracket_width,
                high - low,
                monteith_dew_point_estimate_k,
                wet_bulb_estimate_k,
            )

    raise RuntimeError(
        "experimental ASAE wet-bulb bisection did not converge within "
        f"{_WET_BULB_MAX_ITERATIONS} iterations"
    )


def diagnose_experimental_asae_wet_bulb(
    dry_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> WetBulbSolveDiagnostics:
    """Run the full-bracket ASAE baseline and return solve diagnostics."""
    (
        dry_bulb_temperature_k,
        vapor_pressure_pa,
        atmospheric_pressure_pa,
    ) = _validate_wet_bulb_solver_inputs(
        dry_bulb_temperature_k, vapor_pressure_pa, atmospheric_pressure_pa
    )
    full_width = dry_bulb_temperature_k - ASAE_MIN_TEMPERATURE_K
    return _bisect_wet_bulb_bracket(
        dry_bulb_temperature_k,
        vapor_pressure_pa,
        atmospheric_pressure_pa,
        ASAE_MIN_TEMPERATURE_K,
        dry_bulb_temperature_k,
        used_monteith_bracket=False,
        fallback_used=False,
        bracket_expansions=0,
        proposed_bracket_width_k=full_width,
    )


def experimental_asae_wet_bulb_temperature(
    dry_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> float:
    """Solve the baseline ASAE wet-bulb line by full-bracket bisection in K.

    This existing API remains the experimental baseline. It is restricted to
    255.38--338.72 K and does not use Monteith or Dossat.
    """
    return diagnose_experimental_asae_wet_bulb(
        dry_bulb_temperature_k, vapor_pressure_pa, atmospheric_pressure_pa
    ).wet_bulb_temperature_k


def _monteith_reference(dry_bulb_temperature_k: float) -> tuple[float, float, float]:
    if dry_bulb_temperature_k <= 293.15:
        reference_temperature_k = ASAE_PHASE_TRANSITION_K
        coefficient = 19.65
    else:
        reference_temperature_k = 293.15
        coefficient = 18.0
    reference_pressure_pa = asae_saturation_vapor_pressure(reference_temperature_k)
    return reference_temperature_k, reference_pressure_pa, coefficient


def _monteith_wet_bulb_estimates(
    dry_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> tuple[float, float]:
    """Return Monteith dew-point and secant wet-bulb estimates in K."""
    full_low = ASAE_MIN_TEMPERATURE_K
    full_high = dry_bulb_temperature_k
    reference_temperature, reference_pressure, coefficient = _monteith_reference(
        dry_bulb_temperature_k
    )
    monteith_dew_point = monteith_temperature(
        vapor_pressure_pa,
        reference_temperature,
        reference_pressure,
        coefficient,
    )
    monteith_dew_point = min(max(monteith_dew_point, full_low), full_high)

    dew_residual = _wet_bulb_residual(
        dry_bulb_temperature_k,
        monteith_dew_point,
        vapor_pressure_pa,
        atmospheric_pressure_pa,
    )
    dry_residual = _wet_bulb_residual(
        dry_bulb_temperature_k,
        full_high,
        vapor_pressure_pa,
        atmospheric_pressure_pa,
    )
    residual_difference = dry_residual - dew_residual
    if residual_difference == 0.0:
        wet_bulb_estimate = (monteith_dew_point + full_high) / 2.0
    else:
        wet_bulb_estimate = monteith_dew_point - dew_residual * (
            full_high - monteith_dew_point
        ) / residual_difference
    wet_bulb_estimate = min(max(wet_bulb_estimate, full_low), full_high)
    return monteith_dew_point, wet_bulb_estimate


def diagnose_wet_bulb_asae_monteith_assisted(
    dry_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> WetBulbSolveDiagnostics:
    """Solve the ASAE root using a verified, expandable Monteith bracket.

    Monteith estimates a dew point from an unmodified documented coefficient.
    A secant interpolation of ASAE residuals then estimates wet bulb. The first
    bracket spans at least ±1 K and otherwise ±25% of the estimated dew-to-dry
    interval. It is doubled until a sign change is proven. Failure to obtain a
    narrower valid bracket triggers the exact baseline bracket.
    """
    (
        dry_bulb_temperature_k,
        vapor_pressure_pa,
        atmospheric_pressure_pa,
    ) = _validate_wet_bulb_solver_inputs(
        dry_bulb_temperature_k, vapor_pressure_pa, atmospheric_pressure_pa
    )
    full_low = ASAE_MIN_TEMPERATURE_K
    full_high = dry_bulb_temperature_k
    monteith_dew_point, wet_bulb_estimate = _monteith_wet_bulb_estimates(
        dry_bulb_temperature_k,
        vapor_pressure_pa,
        atmospheric_pressure_pa,
    )

    estimated_span = max(full_high - monteith_dew_point, 0.0)
    half_width = max(
        _MONTEITH_MINIMUM_HALF_WIDTH_K,
        _MONTEITH_INITIAL_BRACKET_FRACTION * estimated_span,
    )
    low = max(full_low, wet_bulb_estimate - half_width)
    high = min(full_high, wet_bulb_estimate + half_width)
    proposed_width = high - low
    expansions = 0

    while True:
        low_residual = _wet_bulb_residual(
            dry_bulb_temperature_k,
            low,
            vapor_pressure_pa,
            atmospheric_pressure_pa,
        )
        high_residual = _wet_bulb_residual(
            dry_bulb_temperature_k,
            high,
            vapor_pressure_pa,
            atmospheric_pressure_pa,
        )
        if _has_sign_change(low_residual, high_residual):
            break
        if low == full_low and high == full_high:
            break
        expansions += 1
        half_width *= 2.0
        low = max(full_low, wet_bulb_estimate - half_width)
        high = min(full_high, wet_bulb_estimate + half_width)

    fallback_used = low == full_low and high == full_high
    used_monteith_bracket = not fallback_used
    if not _has_sign_change(low_residual, high_residual):
        low = full_low
        high = full_high
        fallback_used = True
        used_monteith_bracket = False

    return _bisect_wet_bulb_bracket(
        dry_bulb_temperature_k,
        vapor_pressure_pa,
        atmospheric_pressure_pa,
        low,
        high,
        used_monteith_bracket=used_monteith_bracket,
        fallback_used=fallback_used,
        bracket_expansions=expansions,
        proposed_bracket_width_k=proposed_width,
        monteith_dew_point_estimate_k=monteith_dew_point,
        wet_bulb_estimate_k=wet_bulb_estimate,
    )


def solve_wet_bulb_asae_monteith_assisted(
    dry_bulb_temperature_k: float,
    vapor_pressure_pa: float,
    atmospheric_pressure_pa: float,
) -> float:
    """Return the ASAE wet-bulb root using safeguarded Monteith assistance."""
    return diagnose_wet_bulb_asae_monteith_assisted(
        dry_bulb_temperature_k, vapor_pressure_pa, atmospheric_pressure_pa
    ).wet_bulb_temperature_k
