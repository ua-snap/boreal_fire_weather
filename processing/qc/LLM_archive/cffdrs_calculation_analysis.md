# **CFFDRS Calculation Analysis**

**Date:** February 11, 2026  
**Analyst:** GitHub Copilot  
**Purpose:** This document explains the numeric/math rationale behind a defensive change proposed for `dmc_calc` and `dc_calc` in `04_calculate_cffdrs.py`. It contains (A) an explanation of the problem, (B) the exact code snippets (original vs proposed) with commentary, (C) a step-by-step demonstration of how small differences in `hursmin` can flip domain checks and produce NaNs, and (D) an investigation of how environment / package differences might cause different observed behaviors (NaN vs zero vs other sentinel values).

---

## **A. Problem summary (short)**

- Several operations in `dmc_calc` and `dc_calc` use mathematical functions that require a strictly positive argument: specifically, `np.log(Mr - 20.0)` and `np.log(800.0 / Qr)`.
- If the argument to `np.log()` is zero or negative the result is invalid and becomes `NaN` or `-inf`, which then propagates through downstream calculations (`dmc`, `dc`, `bui`, `fwi`, ...).
- Because the CFFDRS indices are computed sequentially (today depends on yesterday), a single NaN can invalidate many subsequent grid cells/time slices.

## **B. Code examples and the proposed change**

1) Original vulnerable lines (from `dmc_calc` and `dc_calc`) — these are the lines that can produce NaNs when their arguments are invalid:

```python
Pr = np.where(ro > 1.5, 244.72 - 43.43 * np.log(Mr - 20.0), Po)

Dr = 400.0 * np.log(800.0 / Qr)
Dr = np.where(ro <= 2.8, Do, Dr)
Dr[Dr < 0.0] = 0.0
```

Why these are risky:
- `np.log(Mr - 20.0)` requires `Mr - 20.0 > 0`.
- `np.log(800.0 / Qr)` requires `Qr > 0`.
- If either condition is not met for some grid cell, the `np.log` yields `NaN` or `-inf` and that value flows into `Pr` or `Dr` — and onward.

2) Proposed defensive patch (the change I recommended):

```python
# Guard against invalid log arguments
Mr_minus20 = Mr - 20.0
valid_Pr = (ro > 1.5) & (Mr_minus20 > 0.0)
Pr = np.where(valid_Pr, 244.72 - 43.43 * np.log(Mr_minus20), Po)
Pr = np.where(np.isnan(Pr), Po, Pr)
Pr[Pr < 0.0] = 0.0

# For Dr, only compute where Qr > 0
valid_Dr = Qr > 0.0
Dr = np.where(valid_Dr, 400.0 * np.log(800.0 / Qr), Do)
Dr = np.where(ro <= 2.8, Do, Dr)
Dr = np.where(np.isnan(Dr), Do, Dr)
Dr[Dr < 0.0] = 0.0
```

Key differences and rationale:
- We explicitly compute boolean masks for locations where the math is valid.
- We fall back to the previous-day value (`Po` / `Do`) at locations where the formula would be undefined.
- We sanitize any remaining `NaN` with a fallback and clamp negative values to zero afterwards.

Note: `np.where(condition, x, y)` will evaluate `x` and `y` prior to masking, so even with `np.where` alone you may still see warnings / intermediate NaNs during evaluation. The patch above therefore (a) uses a validity mask and (b) replaces NaNs afterwards — the most robust approach is to assign only at masked locations (see the optional alternative below).

Optional stronger pattern (avoid evaluating invalid expressions):

```python
Pr = Po.copy()
mask = (ro > 1.5) & ((Mr - 20.0) > 0)
if np.any(mask):
    computed = 244.72 - 43.43 * np.log(Mr - 20.0)
    Pr[mask] = computed[mask]
Pr[Pr < 0] = 0
```

This approach builds the result from the previous-day values and only computes the log where the inputs are valid, avoiding invalid intermediate results entirely.

## **C. Step-by-step demonstration (conceptual + numeric example)**

Goal: show how a `hursmin` change can lead to a NaN seed in `Pr` or `Dr` and hence in downstream indices.

Important note: the DMC/DC formulas are sequential and coupled. A tiny difference early in the chain (for example in `dmc`) can accumulate and change whether later expressions fall into a valid domain. The demonstration below shows a plausible chain of causation (a small change in `hursmin` modifies a dmc increment, which modifies the previous-day `Po`, which then changes `Mr - 20` or `Qr` sign in subsequent days).

1) Relevant formula pieces (simplified):

- Effective precipitation term for DMC:
  - $re = 0.92\times ro - 1.27$ when $ro > 1.5$, else $re = ro$.
- Mo (a function of previous-day DMC $Po$):
  - $Mo = 20 + \exp(5.6348 - Po/43.43)$
- Mr (only computed when $ro > 1.5$):
  - $Mr = Mo + \dfrac{1000 \cdot re}{48.77 + b \cdot re}$
  - where $b$ itself depends on $Po$.
- Then $Pr = 244.72 - 43.43\log(Mr - 20)$ (only valid if $Mr - 20>0$)

Because $Po$ appears inside an exponential to make $Mo$, even modest changes to $Po$ can significantly change $Mo$ (and thus $Mr$), which can flip whether $Mr - 20$ is positive.

2) Small numeric illustration (hand-worked):

- Start values (example cell / day):
  - $ro = 2.0$ mm (slightly above 1.5 so the precipitation branch uses the Mr formula)
  - Previous-day $Po = 20.0$ (DMC)
  - Compute:
    - $re = 0.92\times 2 - 1.27 = 0.57$
    - $Mo = 20 + e^{(5.6348 - 20/43.43)} \approx 20 + e^{5.1743} \approx 20 + 176.34 = 196.34$
    - For $Po=20$, the `b` formula (Po between 33 and 65 does not apply) gives $b = 100/(0.5 + 0.3\cdot Po) = 100/(0.5 + 6) = 15.38$.
    - $Mr = 196.34 + \dfrac{1000\cdot0.57}{48.77 + 15.38\cdot0.57} \approx 206.25$
    - $Mr - 20 \approx 186.25 > 0$, so the log is valid and $Pr$ gets a finite value.

- Now suppose a small change in `hursmin` (or prior-day `hursmin` chain) causes $Po$ to increase to $Po' = 220.0$ (this is a contrived but plausible cumulative change after several time steps). Then:
    - $Mo' = 20 + e^{(5.6348 - 220/43.43)} = 20 + e^{5.6348 - 5.066} = 20 + e^{0.569} \approx 20 + 1.77 = 21.77$
    - For $Po' = 220$ we use alternate formula for $b$ (Po > 65): $b' = 6.2 \log(Po') - 17.2$; numerically $b' \approx 6.2\cdot5.393 - 17.2 \approx 16.2$.
    - $Mr' = 21.77 + \dfrac{1000\cdot0.57}{48.77 + 16.2\cdot0.57} \approx 31.0$
    - $Mr' - 20 \approx 11.0 > 0$ — still valid here, but the magnitude fell a lot.

The point of these numbers is to show that a modest change in $Po$ can dramatically alter $Mo$ and hence $Mr - 20$. With different starting numbers or larger cumulative effect, a similar chain can push $Mr - 20$ to zero or negative, which will make `np.log(Mr - 20)` invalid and yield `NaN`.

Edge case that directly produces invalid log without large chains:

- If $ro \le 1.5$ then $Mr$ is set to $0.0$ in the original formula. Then $Mr - 20 = -20$ and the $\log(Mr - 20)$ is invalid. In other words, dry days (small precipitation) produce $Mr$=0 which directly causes the invalid log. That is handled in the proposed patch by only applying the log for locations/days where the full condition is mathematically valid.

3) How `hursmin` clipping plays in (your asked scenario):

- Suppose the _original_ dataset had `hursmin` values outside [0,100] (say negative numbers). Those un-clipped values may produce slightly different day-to-day increments of `dmc` and `dc` (the term $K$ in `dmc_calc` contains $(100 - hurs)$). In the sequential integrator, those small differences can accumulate and affect $Po$ or $Do$` on later days.
- Because the key domain checks depend on $Po$, $Do$, and precipitation, a small change in humidity sequences can flip a later-test like `Mr - 20 > 0` or `Qr > 0`. Thus, paradoxically, a dataset with out-of-range humidity could avoid producing a NaN at some location (because the chain resulted in $Mr - 20$ being positive), while the corrected/clamped humidity (0..100) leads the chain to a state where $Mr - 20 \le 0$ and a NaN occurs.

In short: the CFFDRS computations are path-dependent. Correcting one variable (clamping `hursmin`) can change the integrator path and cause a previously-valid computed term to fall into an invalid domain at some later time step.

## **D. Why environments or package versions might show different observed behavior (NaN vs zero vs other values)**

1) Numpy arithmetic rules and `np.seterr`:
- By default, `numpy` will perform floating-point operations and return `NaN` for invalid operations (e.g., `np.log(-1)` -> `nan`) and `inf` for overflow or division by zero (e.g., `1.0/0.0` -> `inf`), while emitting a runtime warning. The default behavior hasn't changed across recent numpy versions, but it is configurable via `np.seterr` and some codebases change the error handling.
- Example behaviors:
  - Default: invalid/log -> `nan`, division-by-zero -> `inf` (with a printed warning if warnings are not suppressed).
  - `np.seterr(invalid='ignore')` will suppress warnings but still produce `nan` values.
  - `np.errstate` context managers can temporarily change behavior.

2) Dask / lazy arrays: evaluation timing differences
- If input arrays are Dask arrays (lazy), evaluation may happen at a later point and some intermediate arrays might never be computed for masked regions. That can sometimes hide or postpone warnings or produce different memory/precision results.

3) Masked arrays and NetCDF `_FillValue` / encoding behavior
- NetCDF files use `_FillValue` and `missing_value` metadata to indicate missing data. When writing NetCDF, libraries sometimes map `NaN` to `_FillValue`. If `_FillValue` is set to 0.0, or if a reader interprets `_FillValue` incorrectly, a written `NaN` can end up read back as `0.0` (or other sentinel).
- Different readers (netCDF4, h5netcdf, xarray's engine options) and different conventions can therefore make a `NaN` look like `0` on some systems but `nan` on others.

4) Type casting and integer/floating behavior
- If arrays are forced into integer types (unlikely here, but possible due to upstream encoding), invalid values may be coerced, truncated, or replaced, producing surprising outcomes. The code in this repo converts to `float32` before writing, which is good, but encoding interactions remain important.

5) Summary of differences you might see across systems
- True math invalidity => `nan` or `inf` in numpy by default.
- If a later I/O layer maps `nan` to a `_FillValue` that a reader interprets as `0.0`, you may see zeros instead of `nan` in the loaded dataset.
- Dask/lazy evaluation may hide or postpone warnings or cause slightly different floating rounding depending on chunking and scheduler.

## **E. Recommendations**

- Keep the defensive guards in `dmc_calc` / `dc_calc` (the `valid_...` boolean masks and masked assignment). That prevents NaN seeds and reduces propagation.
- Prefer the masked-assignment pattern that avoids evaluating invalid expressions at all (see the "Optional stronger pattern" snippet). This both suppresses runtime warnings and avoids temporary NaNs in memory.
- Add diagnostic counters/logging to record how many gridpoints used the fallback (previous-day) values — that helps detect whether the fallback is rare (acceptable) or common (a sign of upstream data problems).
- When writing NetCDF, be explicit with `_FillValue` and attributes to avoid unintentional replacement of NaNs with zeros on read. Use `encoding={'_FillValue': np.nan}` or leave NaNs and allow CF conventions to handle them.