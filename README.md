# psychrometrics-python

`psychrometrics-python` is a small, script-oriented Python library for
calculating fundamental moist-air properties in SI units. It favors explicit,
traceable equations and a simple module-based structure; it is not a GUI or web
application.

## Stable ASHRAE core

The primary implementation in `psychrometrics.py` uses formulations from the
ASHRAE Handbook—Fundamentals. Its current public calculations are:

- saturation vapor pressure;
- vapor pressure from dry-bulb temperature and relative humidity;
- humidity ratio;
- moist-air enthalpy;
- specific volume;
- dew-point temperature.

`calculate_state()` combines these properties in an immutable
`PsychrometricState`. Temperatures are expressed in °C at the public API,
pressures in Pa, relative humidity as a fraction from `0` to `1`, humidity ratio
in kg water/kg dry air, enthalpy in kJ/kg dry air, and specific volume in
m³/kg dry air. Function docstrings document applicable ranges and equations.

## Use as a psychrometric calculator

The stable API has two levels. `calculate_state()` is the simplest option: it
accepts dry-bulb temperature in °C, relative humidity as a fraction from `0`
to `1` (not as a percentage), and atmospheric pressure in Pa:

```python
from psychrometrics import calculate_state

state = calculate_state(
    dry_bulb_temperature_c=25.0,
    relative_humidity=0.50,
    pressure_pa=101325.0,
)

print(f"Saturation vapor pressure: {state.saturation_vapor_pressure_pa:.2f} Pa")
print(f"Vapor pressure: {state.vapor_pressure_pa:.2f} Pa")
print(
    "Humidity ratio: "
    f"{state.humidity_ratio_kg_kg_dry_air:.6f} kg/kg dry air"
)
print(f"Moist-air enthalpy: {state.enthalpy_kj_kg_dry_air:.2f} kJ/kg dry air")
print(
    "Specific volume: "
    f"{state.specific_volume_m3_kg_dry_air:.4f} m³/kg dry air"
)
print(f"Dew-point temperature: {state.dew_point_temperature_c:.2f} °C")
```

Output for 25 °C, 50% relative humidity, and 101325 Pa:

```text
Saturation vapor pressure: 3169.22 Pa
Vapor pressure: 1584.61 Pa
Humidity ratio: 0.009881 kg/kg dry air
Moist-air enthalpy: 50.32 kJ/kg dry air
Specific volume: 0.8580 m³/kg dry air
Dew-point temperature: 13.86 °C
```

The returned `PsychrometricState` also retains the three inputs as
`dry_bulb_temperature_c`, `relative_humidity`, and `pressure_pa`. Its calculated
properties are `saturation_vapor_pressure_pa`, `vapor_pressure_pa`,
`humidity_ratio_kg_kg_dry_air`, `enthalpy_kj_kg_dry_air`,
`specific_volume_m3_kg_dry_air`, and `dew_point_temperature_c`. The experimental
methods described below are separate from this basic calculator workflow.

For more flexible input, `solve_state()` supports multiple independent input
pairs, always with explicit atmospheric pressure. Dry-bulb temperature does
not have to be one of the known properties. Every route returns the same
`PsychrometricState` and uses the stable ASHRAE core.

| Input pair | Supported |
|---|---|
| Tdb + RH | Yes |
| Tdb + Tdp | Yes |
| Tdb + W | Yes |
| Tdb + Twb | Yes |
| h + W | Yes |
| v + W | Yes |
| Tdp + RH | Yes |
| h + RH | Yes |

```python
from psychrometrics import solve_state

from_rh = solve_state(
    25.0, pressure_pa=101325.0, relative_humidity=0.50
)
from_wet_bulb = solve_state(
    25.0, pressure_pa=101325.0, wet_bulb_temperature_c=17.8894
)
from_enthalpy_and_humidity_ratio = solve_state(
    pressure_pa=101325.0,
    enthalpy_kj_kg_dry_air=50.32196,
    humidity_ratio_kg_kg_dry_air=0.00988115,
)
from_dew_point_and_relative_humidity = solve_state(
    pressure_pa=101325.0,
    dew_point_temperature_c=13.86397,
    relative_humidity=0.50,
)
```

Here Tdb, Tdp, and Twb are dry-bulb, dew-point, and wet-bulb temperatures in
°C; RH is a fraction; W is kg water/kg dry air; h is kJ/kg dry air; and v is
m³/kg dry air. The h+W and v+W routes recover Tdb algebraically. Tdp+RH inverts
the ASHRAE saturation-pressure relation, while h+RH uses bounded bisection over
the physical temperature domain. The wet-bulb route uses the ASHRAE
liquid-water or ice equation as appropriate.

Only the eight documented pairs are supported; this is not an unrestricted
any-two-properties solver. `alternative_methods.py` remains separate from all
stable solver routes.

## Experimental methods

`alternative_methods.py` is a complementary research and validation module,
not part of the stable API. It contains:

- ASAE/ASABE saturation-pressure and latent-heat formulations;
- the local Monteith pressure–temperature relation and its inverse;
- an experimental ASAE wet-bulb-temperature solver;
- a safeguarded Monteith-assisted strategy that proposes a narrower search
  bracket while retaining the ASAE equation as the final solved root.

The wet-bulb solvers were compared against PsychroLib over 70 states spanning
two atmospheric pressures. In the published benchmark, the assisted strategy
reduced mean bisection iterations from **15.429 to 12.257**, a **20.41%**
reduction, and the maximum difference between the baseline and assisted roots
was **0.000643165 K**. Both methods remained within approximately **0.04 °C**
of the PsychroLib reference values.

Fewer iterations did not improve measured computational performance in this
case: the Monteith-assisted variant was approximately **6.65% slower** in the
published timing run because estimation and bracket-verification work added
overhead. These timings are platform-dependent and do not establish Monteith
as a faster solver.

See [docs/alternative_methods_validation.md](docs/alternative_methods_validation.md)
for the equations, sources, validation matrix, safeguards, fallback behavior,
and methodological limitations.

## Running the project

Run the minimal stable-core example:

```powershell
python examples/basic_state.py
```

Run the flexible solver example, which reconstructs one state from four input
pairs:

```powershell
python examples/flexible_state_solver.py
```

Run the reproducible experimental wet-bulb benchmark:

```powershell
python examples/benchmark_wet_bulb.py
```

Install the test dependency and run the complete test suite:

```powershell
python -m pip install -r requirements.txt
python -m pytest
```

The current suite contains 93 tests, including independent PsychroLib state and
wet-bulb references and validation of invalid or over-specified solver inputs.

## Repository structure

```text
psychrometrics.py                       Stable ASHRAE calculations and state API
alternative_methods.py                 Experimental ASAE/ASABE and Monteith methods
examples/basic_state.py                Minimal stable-core example
examples/flexible_state_solver.py      Stable flexible-solver example
examples/benchmark_wet_bulb.py         Reproducible experimental benchmark
tests/test_psychrometrics.py           Independent tests of the stable core
tests/test_alternative_methods.py      Experimental-method and fallback tests
docs/alternative_methods_validation.md Methodology and comparison results
requirements.txt                       Test dependencies
```

Air mixing, heating, cooling, humidification, and dehumidification processes
are not currently implemented.
