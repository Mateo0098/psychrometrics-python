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
    solve_state,
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


@pytest.mark.parametrize(
    ("dry_bulb_c", "relative_humidity", "pressure_pa"),
    [
        (-5.0, 0.60, 101_325.0),
        (15.0, 0.80, 80_000.0),
        (25.0, 0.50, 101_325.0),
        (35.0, 0.20, 101_325.0),
        (5.0, 0.95, 80_000.0),
        (25.0, 0.999, 101_325.0),
    ],
)
def test_solve_state_reconstructs_rh_state_from_dew_point_and_humidity_ratio(
    dry_bulb_c: float,
    relative_humidity: float,
    pressure_pa: float,
) -> None:
    reference = calculate_state(dry_bulb_c, relative_humidity, pressure_pa)
    from_rh = solve_state(
        dry_bulb_c,
        pressure_pa=pressure_pa,
        relative_humidity=relative_humidity,
    )
    from_dew_point = solve_state(
        dry_bulb_c,
        pressure_pa=pressure_pa,
        dew_point_temperature_c=reference.dew_point_temperature_c,
    )
    from_humidity_ratio = solve_state(
        dry_bulb_c,
        pressure_pa=pressure_pa,
        humidity_ratio_kg_kg_dry_air=reference.humidity_ratio_kg_kg_dry_air,
    )

    assert from_rh == reference
    assert from_dew_point.relative_humidity == pytest.approx(
        reference.relative_humidity, abs=5e-5
    )
    assert from_dew_point.humidity_ratio_kg_kg_dry_air == pytest.approx(
        reference.humidity_ratio_kg_kg_dry_air, rel=1e-4
    )
    assert from_dew_point.dew_point_temperature_c == pytest.approx(
        reference.dew_point_temperature_c, abs=0.001
    )
    assert from_humidity_ratio.relative_humidity == pytest.approx(
        reference.relative_humidity, abs=1e-12
    )
    assert from_humidity_ratio.humidity_ratio_kg_kg_dry_air == pytest.approx(
        reference.humidity_ratio_kg_kg_dry_air, rel=1e-12
    )


@pytest.mark.parametrize(
    (
        "dry_bulb_c",
        "wet_bulb_c",
        "pressure_pa",
        "psychrolib_relative_humidity",
    ),
    [
        (-5.0, -7.720009955743, 101_325.0, 0.40),
        (0.0, -2.355940749181, 101_325.0, 0.60),
        (25.0, 17.889432148553, 101_325.0, 0.50),
        (40.0, 39.180340157543, 101_325.0, 0.95),
        (15.0, 4.848862086905, 80_000.0, 0.20),
        (5.0, 4.609370025723, 80_000.0, 0.95),
    ],
)
def test_solve_state_from_wet_bulb_against_psychrolib(
    dry_bulb_c: float,
    wet_bulb_c: float,
    pressure_pa: float,
    psychrolib_relative_humidity: float,
) -> None:
    # Independent wet-bulb values are frozen from PsychroLib 2.5.0 SI
    # GetTWetBulbFromRelHum. PsychroLib solves to a 0.001 °C bracket, so the
    # inverse relative humidity is allowed abs=5e-5.
    state = solve_state(
        dry_bulb_c,
        pressure_pa=pressure_pa,
        wet_bulb_temperature_c=wet_bulb_c,
    )
    assert state.relative_humidity == pytest.approx(
        psychrolib_relative_humidity, abs=5e-5
    )


def test_solve_state_accepts_zero_humidity_ratio() -> None:
    state = solve_state(
        20.0,
        pressure_pa=101_325.0,
        humidity_ratio_kg_kg_dry_air=0.0,
    )
    assert state.relative_humidity == 0.0
    assert state.dew_point_temperature_c is None


def test_solve_state_requires_explicit_pressure() -> None:
    with pytest.raises(TypeError, match="pressure_pa"):
        solve_state(25.0, relative_humidity=0.50)


def test_solve_state_rejects_insufficient_or_ambiguous_inputs() -> None:
    with pytest.raises(ValueError, match="exactly two"):
        solve_state(25.0, pressure_pa=101_325.0)
    with pytest.raises(ValueError, match="exactly two"):
        solve_state(pressure_pa=101_325.0, relative_humidity=0.50)
    with pytest.raises(ValueError, match="exactly two"):
        solve_state(
            25.0,
            pressure_pa=101_325.0,
            relative_humidity=0.50,
            dew_point_temperature_c=15.0,
        )
    with pytest.raises(ValueError, match="exactly two"):
        solve_state(
            25.0,
            pressure_pa=101_325.0,
            relative_humidity=0.50,
            dew_point_temperature_c=15.0,
            humidity_ratio_kg_kg_dry_air=0.01,
        )
    with pytest.raises(TypeError, match="unexpected keyword"):
        solve_state(25.0, pressure_pa=101_325.0, enthalpy_kj_kg=50.0)


@pytest.mark.parametrize(
    "arguments",
    [
        {"relative_humidity": 1.01},
        {"dew_point_temperature_c": 25.1},
        {"wet_bulb_temperature_c": 25.1},
        {"humidity_ratio_kg_kg_dry_air": -0.001},
    ],
)
def test_solve_state_rejects_invalid_humidity_inputs(arguments: dict) -> None:
    with pytest.raises(ValueError):
        solve_state(25.0, pressure_pa=101_325.0, **arguments)


def test_solve_state_rejects_invalid_pressure_and_impossible_states() -> None:
    with pytest.raises(ValueError, match="greater than"):
        solve_state(25.0, pressure_pa=0.0, relative_humidity=0.50)
    with pytest.raises(ValueError, match="physically impossible"):
        solve_state(
            25.0,
            pressure_pa=101_325.0,
            humidity_ratio_kg_kg_dry_air=1.0,
        )
    with pytest.raises(ValueError, match="negative humidity ratio"):
        solve_state(
            25.0,
            pressure_pa=101_325.0,
            wet_bulb_temperature_c=-20.0,
        )


_PSYCHROLIB_STATE_REFERENCES = (
    # Tdb °C, RH, P Pa, W kg/kg dry air, Tdp °C, h kJ/kg dry air,
    # v m³/kg dry air, Twb °C. Generated with official PsychroLib 2.5.0 SI
    # CalcPsychrometricsFromRelHum; values are frozen to avoid a test dependency.
    (25.0, 0.50, 101_325.0, 0.009881043690750, 13.863973265965,
     50.321958802185, 0.858043263853, 17.889432148553),
    (35.0, 0.20, 101_325.0, 0.006986454817593, 8.706690674789,
     53.137941707427, 0.882759374606, 18.870393145730),
    (10.0, 0.90, 101_325.0, 0.006858634124814, 8.437213810535,
     27.341014540880, 0.810976854686, 9.157179030144),
    (-5.0, 0.60, 101_325.0, 0.001483174379916, -10.845082425052,
     -1.334374397563, 0.761449454772, -6.790555950399),
    (15.0, 0.80, 80_000.0, 0.010790981847048, 11.581588542887,
     42.379313993001, 1.051827763816, 12.796309703978),
    (40.0, 0.40, 101_325.0, 0.018672483881760, 23.822151902538,
     88.329114989084, 0.913751384567, 27.831560606136),
    (25.0, 0.999, 101_325.0, 0.020060393923260, 24.983218219432,
     76.253853519504, 0.871867189731, 24.987675879895),
)


@pytest.mark.parametrize(
    (
        "dry_bulb_c",
        "relative_humidity",
        "pressure_pa",
        "ratio",
        "dew_point_c",
        "enthalpy",
        "volume",
        "wet_bulb_c",
    ),
    _PSYCHROLIB_STATE_REFERENCES,
)
def test_all_supported_pairs_reconstruct_psychrolib_states(
    dry_bulb_c: float,
    relative_humidity: float,
    pressure_pa: float,
    ratio: float,
    dew_point_c: float,
    enthalpy: float,
    volume: float,
    wet_bulb_c: float,
) -> None:
    routes = (
        solve_state(
            dry_bulb_c,
            pressure_pa=pressure_pa,
            relative_humidity=relative_humidity,
        ),
        solve_state(
            dry_bulb_c,
            pressure_pa=pressure_pa,
            dew_point_temperature_c=dew_point_c,
        ),
        solve_state(
            dry_bulb_c,
            pressure_pa=pressure_pa,
            humidity_ratio_kg_kg_dry_air=ratio,
        ),
        solve_state(
            dry_bulb_c,
            pressure_pa=pressure_pa,
            wet_bulb_temperature_c=wet_bulb_c,
        ),
        solve_state(
            pressure_pa=pressure_pa,
            enthalpy_kj_kg_dry_air=enthalpy,
            humidity_ratio_kg_kg_dry_air=ratio,
        ),
        solve_state(
            pressure_pa=pressure_pa,
            specific_volume_m3_kg_dry_air=volume,
            humidity_ratio_kg_kg_dry_air=ratio,
        ),
        solve_state(
            pressure_pa=pressure_pa,
            dew_point_temperature_c=dew_point_c,
            relative_humidity=relative_humidity,
        ),
        solve_state(
            pressure_pa=pressure_pa,
            enthalpy_kj_kg_dry_air=enthalpy,
            relative_humidity=relative_humidity,
        ),
    )

    # Algebraic routes are effectively exact; 0.001 °C covers the documented
    # temperature width of both bounded inversions and PsychroLib's wet-bulb
    # reference solve. Derived-property tolerances include that temperature
    # uncertainty and rounded frozen reference digits.
    for state in routes:
        assert state.dry_bulb_temperature_c == pytest.approx(
            dry_bulb_c, abs=0.001
        )
        assert state.relative_humidity == pytest.approx(
            relative_humidity, abs=5e-5
        )
        assert state.humidity_ratio_kg_kg_dry_air == pytest.approx(
            ratio, rel=2e-4
        )
        assert state.enthalpy_kj_kg_dry_air == pytest.approx(
            enthalpy, abs=0.003
        )
        assert state.specific_volume_m3_kg_dry_air == pytest.approx(
            volume, abs=5e-6
        )
        assert state.dew_point_temperature_c == pytest.approx(
            dew_point_c, abs=0.002
        )


@pytest.mark.parametrize("relative_humidity", [0.0, -0.1, 1.01])
def test_dew_point_relative_humidity_pair_requires_positive_valid_rh(
    relative_humidity: float,
) -> None:
    with pytest.raises(ValueError, match="relative_humidity"):
        solve_state(
            pressure_pa=101_325.0,
            dew_point_temperature_c=10.0,
            relative_humidity=relative_humidity,
        )


def test_new_pairs_reject_invalid_and_unsupported_inputs() -> None:
    with pytest.raises(ValueError):
        solve_state(
            pressure_pa=101_325.0,
            enthalpy_kj_kg_dry_air=-1000.0,
            humidity_ratio_kg_kg_dry_air=0.0,
        )
    with pytest.raises(ValueError, match="greater than"):
        solve_state(
            pressure_pa=101_325.0,
            specific_volume_m3_kg_dry_air=0.0,
            humidity_ratio_kg_kg_dry_air=0.01,
        )
    with pytest.raises(ValueError, match="at least"):
        solve_state(
            pressure_pa=101_325.0,
            specific_volume_m3_kg_dry_air=0.85,
            humidity_ratio_kg_kg_dry_air=-0.01,
        )
    with pytest.raises(ValueError, match="unsupported"):
        solve_state(
            pressure_pa=101_325.0,
            enthalpy_kj_kg_dry_air=50.0,
            specific_volume_m3_kg_dry_air=0.85,
        )


def test_new_pairs_reject_pressure_conflicts_and_missing_solutions() -> None:
    with pytest.raises(ValueError, match=">= pressure_pa"):
        solve_state(
            pressure_pa=80_000.0,
            dew_point_temperature_c=100.0,
            relative_humidity=0.50,
        )
    with pytest.raises(ValueError, match="outside"):
        solve_state(
            pressure_pa=1_000_000.0,
            dew_point_temperature_c=150.0,
            relative_humidity=0.01,
        )
    with pytest.raises(ValueError, match="no solution"):
        solve_state(
            pressure_pa=101_325.0,
            enthalpy_kj_kg_dry_air=-1000.0,
            relative_humidity=0.50,
        )


def test_enthalpy_relative_humidity_reports_non_convergence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(psychrometrics, "_STATE_ROOT_MAX_ITERATIONS", 0)
    with pytest.raises(RuntimeError, match="did not converge"):
        solve_state(
            pressure_pa=101_325.0,
            enthalpy_kj_kg_dry_air=50.321958802185,
            relative_humidity=0.50,
        )
