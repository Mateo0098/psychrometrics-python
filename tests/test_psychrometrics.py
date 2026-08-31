"""Tests for the SI psychrometric calculation core."""

import math

import pytest

import psychrometrics
from psychrometrics import (
    calculate_state,
    dew_point_temperature,
    humidity_ratio,
    moist_air_enthalpy,
    saturation_vapor_pressure,
    specific_volume,
    vapor_pressure,
)


@pytest.mark.parametrize(
    ("temperature_c", "reference_pressure_pa", "tolerance"),
    [
        (-20.0, 103.24, 0.0003),
        (25.0, 3169.7, 0.0003),
        (150.0, 476_101.4, 0.0003),
    ],
)
def test_saturation_pressure_against_ashrae_table_3(
    temperature_c: float,
    reference_pressure_pa: float,
    tolerance: float,
) -> None:
    # Independent values: ASHRAE Handbook—Fundamentals (2017), ch. 1,
    # Table 3, as transcribed in the official PsychroLib 2.5.0 SI tests.
    # Units are °C and Pa. ASHRAE states 300 ppm accuracy, hence rel=0.0003.
    assert saturation_vapor_pressure(temperature_c) == pytest.approx(
        reference_pressure_pa, rel=tolerance
    )


def test_saturation_pressure_at_water_triple_point() -> None:
    # ASHRAE Handbook—Fundamentals gives approximately 611.657 Pa.
    assert saturation_vapor_pressure(0.01) == pytest.approx(611.657, rel=2e-6)


def test_saturation_pressure_at_25_c() -> None:
    assert saturation_vapor_pressure(25.0) == pytest.approx(3169.216, rel=2e-6)


def test_saturation_pressure_is_continuous_at_triple_point() -> None:
    below = saturation_vapor_pressure(0.009)
    at_triple_point = saturation_vapor_pressure(0.01)
    above = saturation_vapor_pressure(0.011)

    assert below < at_triple_point < above
    # Across 0.001 °C on either side, each branch changes by about 0.05 Pa.
    # A 0.06 Pa bound detects a branch discontinuity without requiring equal
    # slopes for saturation over ice and over liquid water.
    assert at_triple_point - below < 0.06
    assert above - at_triple_point < 0.06


def test_humidity_ratio_against_psychrolib_excel_reference() -> None:
    # Independent PsychroLib 2.5.0 SI test case, calculated in Excel:
    # Pv=3169.7 Pa, total pressure=95461 Pa, W=kg water/kg dry air.
    # The source comparison uses one part per million relative tolerance.
    assert humidity_ratio(3169.7, 95_461.0) == pytest.approx(
        0.0213603998047487, rel=1e-6
    )


def test_moist_air_enthalpy_against_psychrolib_excel_reference() -> None:
    # Independent PsychroLib 2.5.0 SI test case, calculated in Excel:
    # t=30 °C, W=0.02 kg/kg, h=81316 J/kg dry air = 81.316 kJ/kg.
    # PsychroLib uses the ASHRAE correlation accuracy, rel=0.0003.
    assert moist_air_enthalpy(30.0, 0.02) == pytest.approx(
        81.316, rel=0.0003
    )


def test_saturated_enthalpy_against_ashrae_table_2() -> None:
    # Independent ASHRAE Fundamentals (2017), ch. 1, Table 2 values as
    # transcribed in PsychroLib 2.5.0: at 25 °C, W=0.020173 kg/kg and
    # h=76504 J/kg dry air = 76.504 kJ/kg. The source test uses rel=0.01
    # because agreement between the correlation and rounded table is ~1%.
    assert moist_air_enthalpy(25.0, 0.020173) == pytest.approx(
        76.504, rel=0.01
    )


def test_specific_volume_against_psychrolib_excel_reference() -> None:
    # Independent PsychroLib 2.5.0 SI test case, calculated in Excel:
    # t=30 °C, W=0.02 kg/kg, p=95461 Pa, v in m³/kg dry air.
    assert specific_volume(30.0, 0.02, 95_461.0) == pytest.approx(
        0.940855374352943, rel=0.0003
    )


@pytest.mark.parametrize(
    ("reference_vapor_pressure_pa", "reference_dew_point_c"),
    [
        (103.24, -20.0),
        (872.6, 5.0),
        (12_351.3, 50.0),
    ],
)
def test_dew_point_against_ashrae_table_3_pressures(
    reference_vapor_pressure_pa: float,
    reference_dew_point_c: float,
) -> None:
    # Inputs are independent saturation pressures from ASHRAE Fundamentals
    # (2017), ch. 1, Table 3, via the official PsychroLib 2.5.0 SI tests.
    # abs=0.01 °C allows for the table pressures' published rounding.
    assert dew_point_temperature(reference_vapor_pressure_pa) == pytest.approx(
        reference_dew_point_c, abs=0.01
    )


def test_basic_state_at_25_c_and_50_percent_rh() -> None:
    state = calculate_state(25.0, 0.50, 101_325.0)

    assert state.vapor_pressure_pa == pytest.approx(1584.608, rel=2e-6)
    assert state.humidity_ratio_kg_kg_dry_air == pytest.approx(
        0.009881, rel=5e-5
    )
    assert state.enthalpy_kj_kg_dry_air == pytest.approx(50.32, abs=0.02)
    assert state.specific_volume_m3_kg_dry_air == pytest.approx(
        0.8580, abs=0.0002
    )
    assert state.dew_point_temperature_c == pytest.approx(13.864, abs=0.002)


@pytest.mark.parametrize("temperature_c", [-80.0, -10.0, 0.01, 25.0, 100.0])
def test_dew_point_inverts_saturation_pressure(temperature_c: float) -> None:
    # Internal round-trip check at the public dew-point tolerance. Independent
    # absolute reference cases are tested separately above.
    pressure = saturation_vapor_pressure(temperature_c)
    assert dew_point_temperature(pressure) == pytest.approx(
        temperature_c, abs=0.001
    )


def test_zero_relative_humidity_has_no_finite_dew_point() -> None:
    state = calculate_state(20.0, 0.0)
    assert state.humidity_ratio_kg_kg_dry_air == 0.0
    assert state.dew_point_temperature_c is None


@pytest.mark.parametrize("relative_humidity", [-0.01, 1.01, math.inf])
def test_relative_humidity_validation(relative_humidity: float) -> None:
    with pytest.raises(ValueError):
        vapor_pressure(25.0, relative_humidity)


@pytest.mark.parametrize("temperature_c", [-100.01, 200.01, math.nan])
def test_temperature_validation(temperature_c: float) -> None:
    with pytest.raises(ValueError):
        saturation_vapor_pressure(temperature_c)


def test_pressure_and_humidity_ratio_validation() -> None:
    with pytest.raises(ValueError, match="less than"):
        humidity_ratio(101_325.0, 101_325.0)
    with pytest.raises(ValueError, match="greater than"):
        specific_volume(25.0, 0.01, 0.0)
    with pytest.raises(ValueError, match="at least"):
        moist_air_enthalpy(25.0, -0.01)


def test_boolean_is_not_accepted_as_a_number() -> None:
    with pytest.raises(TypeError):
        saturation_vapor_pressure(True)


def test_dew_point_reports_non_convergence(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(psychrometrics, "_DEW_POINT_MAX_ITERATIONS", 0)
    with pytest.raises(RuntimeError, match="did not converge"):
        dew_point_temperature(872.6)
