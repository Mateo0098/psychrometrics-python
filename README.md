# psychrometrics-python

A small, script-oriented Python implementation of fundamental moist-air
properties in SI units. The current phase intentionally contains no GUI, web
application, or air-treatment process models.

## Implemented calculations

- Saturation vapor pressure from -100 to 200 °C.
- Vapor partial pressure from dry-bulb temperature and relative humidity.
- Humidity ratio in kg water/kg dry air.
- Moist-air enthalpy in kJ/kg dry air.
- Specific volume in m³/kg dry air.
- Dew-point temperature obtained by numerically inverting saturation pressure.
- An immutable `PsychrometricState` assembled from dry-bulb temperature,
  relative humidity, and total pressure.

Relative humidity is always a fraction from `0` to `1`, pressure is in Pa, and
temperature is in °C. Functions validate finite numbers and their physical or
correlation ranges. At exactly zero relative humidity, the state's dew point is
`None` because zero vapor pressure has no finite dew point in this model.

## Example

```powershell
python examples/basic_state.py
```

```python
from psychrometrics import calculate_state

state = calculate_state(25.0, 0.50, 101_325.0)
print(state.humidity_ratio_kg_kg_dry_air)
print(state.enthalpy_kj_kg_dry_air)
```

## Tests

Install the test dependency and run:

```powershell
python -m pip install -r requirements.txt
python -m pytest
```

## Equation sources and scope

The equations and constants are from ASHRAE Handbook—Fundamentals (2021),
Chapter 1, *Psychrometrics*. Function docstrings identify the applicable
equation and units. The open-source [PsychroLib implementation](https://github.com/psychrometrics/psychrolib)
provides a publicly inspectable implementation of the same ASHRAE SI
formulation.

ASHRAE is the repository's primary formulation. The ASAE/ASABE D271 model used
by the earlier university project is retained as a historical and comparative
validation reference; it is not exposed as a second official calculation API,
and its results are not assumed to be interchangeable with ASHRAE results.

The test suite includes independent reference points transcribed by the
official PsychroLib 2.5.0 SI tests from Tables 2 and 3 of ASHRAE
Handbook—Fundamentals (2017), plus its independently calculated Excel cases.
Tests state source units and tolerances alongside each reference case.

Not yet implemented: wet-bulb temperature, air mixing, heating, cooling,
humidification, or dehumidification. Those are intentionally deferred until
the core has been independently validated over its intended operating range.
