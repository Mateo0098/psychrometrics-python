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

Run the reproducible experimental wet-bulb benchmark:

```powershell
python examples/benchmark_wet_bulb.py
```

Install the test dependency and run the complete test suite:

```powershell
python -m pip install -r requirements.txt
python -m pytest
```

## Repository structure

```text
psychrometrics.py                       Stable ASHRAE calculations and state API
alternative_methods.py                 Experimental ASAE/ASABE and Monteith methods
examples/basic_state.py                Minimal stable-core example
examples/benchmark_wet_bulb.py         Reproducible experimental benchmark
tests/test_psychrometrics.py           Independent tests of the stable core
tests/test_alternative_methods.py      Experimental-method and fallback tests
docs/alternative_methods_validation.md Methodology and comparison results
requirements.txt                       Test dependencies
```

Air mixing, heating, cooling, humidification, and dehumidification processes
are not currently implemented.
