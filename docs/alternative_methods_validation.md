# Validation of historical alternative psychrometric methods

## Scope and sources

The stable calculation core in `psychrometrics.py` remains based on ASHRAE.
The isolated `alternative_methods.py` module reconstructs selected equations
from the original university report, `docs/Copia de Proyecto A.docx`, for
comparison only. The report identifies its principal sources as ASAE D271.2
APR1979 (R2014) and Monteith and Unsworth (2013), equation 2.2.10. It also
mentions ASAE D271.3 in its methodology, although its bibliography lists
D271.2.

The report is the primary source for the transcribed constants, units, and
ranges below. PsychroLib 2.5.0 in SI mode, commit
`3066345dc8cf91bf59134147cf917f982c1fce13`, provides independent ASHRAE
wet-bulb comparison values. No constants were fitted to those values.

## Problem addressed by the original scheme

The original interactive MATLAB model attempted to determine a moist-air state
from any usable pair of properties. Its wet-bulb line depends on saturation
pressure and phase-change latent heat at the unknown wet-bulb temperature.
The report's conditional, equation-by-equation GUI could not select the latent
heat branch before that temperature was known, so it reported a circular
dependency.

Monteith was added to provide a direct, locally anchored relationship between
temperature and vapor pressure. Dossat was proposed later as a possible way to
obtain a latent-energy term more directly.

## Reconstructed equations

All ASAE temperatures are absolute temperatures in K and pressures are Pa.

### ASAE saturation pressure

For 255.38--273.16 K:

```text
ln(Ps) = 31.9602 - 6270.3605/T - 0.46057 ln(T)
```

For 273.16--533.16 K:

```text
ln(Ps/R) = (A + B T + C T² + D T³ + E T⁴) / (F T - G T²)
```

with `R=22105649.25`, `A=-27405.526`, `B=97.5413`,
`C=-0.146244`, `D=0.12558e-3`, `E=-0.48502e-7`, `F=4.34903`,
and `G=0.39381e-2`.

The inverse polynomial is implemented only over its stated pressure range,
620.52--4,688,396 Pa. The inverse is not used outside that narrower range.

Comparison with the stable ASHRAE implementation:

| T (°C) | ASAE (Pa) | ASHRAE (Pa) | Error (Pa) | Relative error |
|---:|---:|---:|---:|---:|
| -10.00 | 261.2521 | 259.9029 | +1.3492 | +0.5191% |
| 0.01 | 614.9098 | 611.6570 | +3.2528 | +0.5318% |
| 25.00 | 3165.1199 | 3169.2165 | -4.0966 | -0.1293% |
| 40.00 | 7370.9519 | 7383.4600 | -12.5081 | -0.1694% |

### Monteith relation and inverse

The report gives:

```text
e(T) = e(T*) exp[A (T - T*) / T]
```

It states `A=19.65` for reference temperatures from 273 to 293 K and
`A=18` for reference temperatures from 293 to 313 K. It does not supply
`e(T*)`, so the implementation requires both `T*` and `e(T*)` explicitly.
The inverse used here is algebra, not a new empirical equation:

```text
T = T* / (1 - ln(e/e*) / A)
```

Using ASHRAE tabular anchors of 611.657 Pa at 273.16 K and 2338.8 Pa at
293.15 K gives:

| T (°C) | Monteith (Pa) | ASHRAE (Pa) | Error (Pa) | Relative error |
|---:|---:|---:|---:|---:|
| -10.00 | 289.6585 | 259.9029 | +29.7556 | +11.4487% |
| 0.01 | 611.6570 | 611.6570 | -0.00002 | -0.000004% |
| 20.00 | 2335.8053 | 2338.8037 | -2.9985 | -0.1282% |
| 25.00 | 3162.9320 | 3169.2165 | -6.2844 | -0.1983% |
| 35.00 | 5617.2074 | 5627.8194 | -10.6121 | -0.1886% |
| 40.00 | 7383.4927 | 7383.4600 | +0.0327 | +0.00044% |

The subzero row extrapolates away from the 273.16 K anchor because the report
does not document a target-temperature range. Its 11.45% error demonstrates
that this use is not defensible as a general subzero saturation model.

### ASAE latent heat

The complete equations reconstructed from the report are:

```text
h_ig = 2839683.144 - 212.56384 (T - 255.38)       [J/kg]
```

for sublimation over 255.38--273.16 K, and:

```text
h_fg = 2502535.259 - 2385.76424 (T - 273.16)      [J/kg]
```

for vaporization over 273.16--338.72 K. No unit conversion is applied because
the report already specifies J/kg.

The displayed high-temperature branch, 338.72--533.16 K, contains
`7,329,155,978,000 - 15,995,964.08 T²` but no visible square root or outer
exponent. Taken literally it is not a plausible latent heat in J/kg. It is not
implemented.

### Wet-bulb line

The report gives:

```text
P_swb - P_v = B' (T_wb - T)

B' = 1006.9254 (P_swb - P_atm) (1 + 0.15577 P_v/P_atm)
     / (0.62194 h'_fg)
```

The implementation evaluates this residual and solves it by bisection within
255.38--338.72 K, selecting the documented sublimation or vaporization branch
at each candidate temperature. This numerical interpretation removes the GUI's
branch-selection circularity without adding constants.

Comparison against PsychroLib 2.5.0 at 101325 Pa:

| Dry bulb | RH | PsychroLib wet bulb | Alternative wet bulb | Absolute error | Relative error* |
|---:|---:|---:|---:|---:|---:|
| -5 °C | 70% | -6.3348 °C | -6.3566 °C | -0.0218 °C | -0.344% |
| 0 °C | 60% | -2.3559 °C | -2.3842 °C | -0.0283 °C | -1.200% |
| 25 °C | 50% | 17.8894 °C | 17.9018 °C | +0.0124 °C | +0.069% |
| 40 °C | 40% | 27.8316 °C | 27.8534 °C | +0.0218 °C | +0.078% |

\*Relative temperature error in Celsius is included only to reproduce the
requested comparison; absolute error is the meaningful metric near 0 °C.

## Monteith-assisted wet-bulb solver

### Method

The full-range ASAE bisection remains the baseline and its public experimental
function is unchanged. The assisted variant does not accept the Monteith value
as wet bulb and does not alter the ASAE residual, tolerance, or final root.
Instead it performs these safeguarded steps:

1. Select the documented Monteith coefficient without fitting it: `19.65`
   with the 273.16 K ASAE anchor for dry bulb up to 293.15 K, otherwise `18`
   with the 293.15 K ASAE anchor.
2. Invert Monteith to estimate dew point.
3. Evaluate the ASAE residual at estimated dew point and dry bulb, then use a
   secant interpolation only to propose a wet-bulb center.
4. Propose a bracket of at least ±1 K, or ±25% of the estimated dew-to-dry
   interval when larger. These are search-policy values, not modified
   psychrometric constants.
5. Prove that the proposed endpoints bracket an ASAE sign change. If not,
   double the half-width until a sign change is found or the full baseline
   interval is recovered.
6. Run the same ASAE bisection with the same 0.001 K stopping width. A forced
   bad-estimate test verifies fallback to the full interval and the baseline
   root.

`WetBulbSolveDiagnostics` records the temperature, iteration count, proposed,
initial, and final bracket widths, expansions, Monteith-bracket use, fallback,
and both intermediate estimates. The final root always comes from the ASAE
equation.

### Validation matrix and precision

The matrix contains all 70 combinations of:

- dry bulb: -5, 0, 5, 15, 25, 35, and 40 °C;
- relative humidity: 20, 40, 60, 80, and 95%;
- atmospheric pressure: 101325 and 80000 Pa.

No state was excluded. PsychroLib 2.5.0 SI reference values are frozen in the
test suite, so validation does not require a runtime dependency.

| Metric over 70 states | Baseline | Monteith-assisted |
|---|---:|---:|
| Maximum error vs PsychroLib | 0.035802 °C | 0.035986 °C |
| Mean bisection iterations | 15.429 | 12.257 |
| Median bisection iterations | 16 | 12 |
| Iteration range | 14--16 | 11--14 |
| Natural fallback cases | n/a | 0 (0%) |
| Cases with more assisted iterations | n/a | 0 |

The maximum baseline/assisted result difference is `0.000643165 K`, below the
shared 0.001 K convergence width. The mean per-case iteration reduction is
20.41%. In 37 of 70 states the assisted rounding point is marginally farther
from PsychroLib than the baseline point, but both remain within 0.04 °C. Thus
the assistance preserves practical accuracy but does not improve it.

### Timing

`examples/benchmark_wet_bulb.py` uses only the standard library and the local
modules. Its default run times 200 repetitions of all 70 cases in five rounds,
alternating method order and reporting the median round. On the development
Windows/Python 3.14.4 environment (14,000 solves per method per round), the
two complete observed runs were:

| Run | Baseline median | Assisted median | Assisted time change |
|---:|---:|---:|---:|
| 1 | 0.919919 s | 1.007807 s | 9.55% slower |
| 2 (final validation) | 1.305608 s | 1.374077 s | 5.24% slower |

The assisted method was 5.24--9.55% slower in these complete runs despite
20.41% fewer bisection iterations. Monteith inversion, two extra ASAE residual
evaluations, secant construction, bracket verification, and diagnostics cost
more than the saved iterations on this workload. Individual runs showed normal
timing variation, so absolute times are platform-specific; the benchmark must
be rerun on a target system before making performance claims.

The data therefore support an algorithmic reduction in bisection iterations,
but do **not** justify describing Monteith assistance as a faster solver. It
remains an experimental comparison path with safeguarded fallback.

## Dossat ambiguity

The report reproduces only the relation `h_L = w h_w` for one pound of dry
air, where `w` is lb water/lb dry air and `h_w` is saturated-vapor enthalpy at
the dew-point temperature in Btu/lb. This is an energy accounting identity,
not a standalone correlation for `h_fg` or `h_ig`.

The required function or table `h_w(T_dp)` is absent. The report itself says a
thermodynamic database would be needed and that the proposed SI result was not
validated. Although `1 Btu/lb = 2326 J/kg` is a valid unit conversion, applying
it cannot reconstruct the missing property values. Consequently no Dossat
function is implemented and no numerical accuracy claim is made.

## Technical assessment

The reproducible ASAE equations plus numerical branch selection can solve the
tested wet-bulb states without Monteith or Dossat, and agree with PsychroLib to
within 0.03 °C in the four sampled states. This supports the wet-bulb line as a
useful historical validation method over the implemented range, not yet as a
replacement for the ASHRAE method.

Monteith is accurate locally when supplied with a suitable explicit anchor,
but it is not a global saturation correlation and is poor in the tested
subzero extrapolation. Dossat cannot be evaluated reproducibly from the source
material. Therefore the complete proposed combination “ASAE + Monteith +
Dossat” is not presently a technically defensible general wet-bulb strategy.
The stable public API and `PsychrometricState` remain unchanged.
