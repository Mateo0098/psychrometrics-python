"""Tests for the SI psychrometric calculation core."""

import math

import pytest

from psychrometrics import (
    calculate_state,
    dew_point_temperature,
    humidity_ratio,
    moist_air_enthalpy,
    saturation_vapor_pressure,
    specific_volume,
    vapor_pressure,
)


def test_saturation_pressure_at_water_triple_point() -> None:
    # ASHRAE Handbook—Fundamentals gives approximately 611.657 Pa.
    assert saturation_vapor_pressure(0.01) == pytest.approx(611.657, rel=2e-6)


def test_saturation_pressure_at_25_c() -> None:
    assert saturation_vapor_pressure(25.0) == pytest.approx(3169.216, rel=2e-6)


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
    pressure = saturation_vapor_pressure(temperature_c)
    assert dew_point_temperature(pressure) == pytest.approx(
        temperature_c, abs=1e-9
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
