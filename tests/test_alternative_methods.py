"""Validation of isolated historical and experimental psychrometric methods."""

import pytest

from alternative_methods import (
    asae_latent_heat_sublimation,
    asae_latent_heat_vaporization,
    asae_saturation_temperature,
    asae_saturation_vapor_pressure,
    asae_wet_bulb_line_residual,
    experimental_asae_wet_bulb_temperature,
    monteith_temperature,
    monteith_vapor_pressure,
)
from psychrometrics import saturation_vapor_pressure


@pytest.mark.parametrize(
    ("temperature_c", "expected_asae_pressure_pa", "relative_error_percent"),
    [
        (-10.0, 261.2520544650352, 0.5191129821152218),
        (0.01, 614.9098456687155, 0.5318047775349832),
        (25.0, 3165.119918891563, -0.12926069552702257),
        (40.0, 7370.951868410325, -0.16940757531795914),
    ],
)
def test_asae_saturation_pressure_and_documented_ashrae_difference(
    temperature_c: float,
    expected_asae_pressure_pa: float,
    relative_error_percent: float,
) -> None:
    # ASAE values reproduce equation 1 and its SI coefficients in the original
    # project DOCX. The comparison uses the stable ASHRAE API, in Pa. Expected
    # relative errors are reported explicitly rather than fitted away.
    asae_pressure = asae_saturation_vapor_pressure(temperature_c + 273.15)
    ashrae_pressure = saturation_vapor_pressure(temperature_c)
    assert asae_pressure == pytest.approx(expected_asae_pressure_pa, rel=1e-12)
    assert 100.0 * (asae_pressure / ashrae_pressure - 1.0) == pytest.approx(
        relative_error_percent, abs=1e-10
    )


@pytest.mark.parametrize("temperature_c", [5.0, 20.0, 40.0, 100.0])
def test_asae_inverse_saturation_polynomial(temperature_c: float) -> None:
    # Original DOCX equation 2: pressure in Pa and result in K, valid only for
    # 620.52 <= Ps <= 4,688,396 Pa. Its rounded polynomial is allowed 0.04 K.
    pressure = asae_saturation_vapor_pressure(temperature_c + 273.15)
    assert asae_saturation_temperature(pressure) == pytest.approx(
        temperature_c + 273.15, abs=0.04
    )


def test_asae_inverse_rejects_undocumented_low_pressure() -> None:
    with pytest.raises(ValueError, match="620.52"):
        asae_saturation_temperature(500.0)


@pytest.mark.parametrize(
    (
        "temperature_c",
        "reference_temperature_k",
        "reference_pressure_pa",
        "coefficient",
        "expected_pressure_pa",
        "absolute_error_pa",
        "relative_error_percent",
    ),
    [
        # Subzero use is an explicit extrapolation; the DOCX gives no target
        # range. Its large error is retained to demonstrate that limitation.
        (-10.0, 273.16, 611.657, 19.65, 289.6584713044197, 29.755606352240648, 11.448741189411482),
        (0.01, 273.16, 611.657, 19.65, 611.6569999999975, -0.000024390879048, -0.000003987672515),
        (20.0, 273.16, 611.657, 19.65, 2335.805250032521, -2.9984500414602735, -0.1282044338037136),
        (25.0, 293.15, 2338.8, 18.0, 3162.932034556652, -6.284435586975633, -0.1982961923295501),
        (35.0, 293.15, 2338.8, 18.0, 5617.207390439753, -10.61205610048728, -0.1885642601240689),
        (40.0, 293.15, 2338.8, 18.0, 7383.492690746618, 0.032681760498235235, 0.0004426347601116376),
    ],
)
def test_monteith_pressure_and_explicit_error_against_ashrae(
    temperature_c: float,
    reference_temperature_k: float,
    reference_pressure_pa: float,
    coefficient: float,
    expected_pressure_pa: float,
    absolute_error_pa: float,
    relative_error_percent: float,
) -> None:
    # Monteith & Unsworth (2013), equation 2.2.10, as transcribed in the DOCX.
    # Inputs are K and Pa; A is dimensionless. Anchors are ASHRAE table values.
    monteith_pressure = monteith_vapor_pressure(
        temperature_c + 273.15,
        reference_temperature_k,
        reference_pressure_pa,
        coefficient,
    )
    ashrae_pressure = saturation_vapor_pressure(temperature_c)
    assert monteith_pressure == pytest.approx(expected_pressure_pa, rel=1e-12)
    assert monteith_pressure - ashrae_pressure == pytest.approx(
        absolute_error_pa, abs=1e-9
    )
    assert 100.0 * (monteith_pressure / ashrae_pressure - 1.0) == pytest.approx(
        relative_error_percent, abs=1e-10
    )


@pytest.mark.parametrize("temperature_k", [263.15, 273.16, 293.15, 313.15])
def test_monteith_algebraic_inverse(temperature_k: float) -> None:
    coefficient = 19.65 if temperature_k <= 293.15 else 18.0
    reference_temperature_k = 273.16 if coefficient == 19.65 else 293.15
    reference_pressure_pa = 611.657 if coefficient == 19.65 else 2338.8
    pressure = monteith_vapor_pressure(
        temperature_k,
        reference_temperature_k,
        reference_pressure_pa,
        coefficient,
    )
    assert monteith_temperature(
        pressure,
        reference_temperature_k,
        reference_pressure_pa,
        coefficient,
    ) == pytest.approx(temperature_k, abs=1e-12)


@pytest.mark.parametrize(
    ("temperature_k", "expected_j_kg"),
    [
        (263.15, 2_838_031.5229632),
        (273.16, 2_835_903.7589248),
    ],
)
def test_asae_latent_heat_of_sublimation(
    temperature_k: float, expected_j_kg: float
) -> None:
    # ASAE equation 3 in the DOCX; T is K, output is J/kg, valid 255.38--273.16 K.
    assert asae_latent_heat_sublimation(temperature_k) == pytest.approx(
        expected_j_kg, rel=1e-12
    )


@pytest.mark.parametrize(
    ("temperature_k", "expected_j_kg"),
    [
        (293.15, 2_454_843.8318424),
        (313.15, 2_407_128.5470424),
        (338.15, 2_347_484.4410424),
    ],
)
def test_asae_linear_latent_heat_of_vaporization(
    temperature_k: float, expected_j_kg: float
) -> None:
    # Reproducible branch of ASAE equation 4 in the DOCX; output J/kg and
    # documented range 273.16--338.72 K.
    assert asae_latent_heat_vaporization(temperature_k) == pytest.approx(
        expected_j_kg, rel=1e-12
    )


def test_ambiguous_high_temperature_latent_branch_is_not_extrapolated() -> None:
    with pytest.raises(ValueError, match="338.72"):
        asae_latent_heat_vaporization(350.0)


@pytest.mark.parametrize(
    (
        "dry_bulb_c",
        "vapor_pressure_pa",
        "psychrolib_wet_bulb_c",
        "expected_alternative_wet_bulb_c",
    ),
    [
        (-5.0, 281.23488573516084, -6.334791830147736, -6.356578674316438),
        (0.0, 366.69214253446074, -2.3559407491807716, -2.3842103576660065),
        (25.0, 1584.6082350718139, 17.889432148552928, 17.901795883178806),
        (40.0, 2953.384003594448, 27.831560606135824, 27.853367385864317),
    ],
)
def test_experimental_wet_bulb_against_psychrolib(
    dry_bulb_c: float,
    vapor_pressure_pa: float,
    psychrolib_wet_bulb_c: float,
    expected_alternative_wet_bulb_c: float,
) -> None:
    # Independent references were generated with official PsychroLib 2.5.0
    # (ASHRAE SI), p=101325 Pa, at RH=(0.70, 0.60, 0.50, 0.40), respectively.
    # Temperatures are °C here; the alternative API itself consumes/returns K.
    wet_bulb_c = experimental_asae_wet_bulb_temperature(
        dry_bulb_c + 273.15, vapor_pressure_pa, 101_325.0
    ) - 273.15
    assert wet_bulb_c == pytest.approx(expected_alternative_wet_bulb_c, abs=0.001)
    assert wet_bulb_c == pytest.approx(psychrolib_wet_bulb_c, abs=0.04)
    assert asae_wet_bulb_line_residual(
        dry_bulb_c + 273.15,
        wet_bulb_c + 273.15,
        vapor_pressure_pa,
        101_325.0,
    ) == pytest.approx(0.0, abs=0.1)
