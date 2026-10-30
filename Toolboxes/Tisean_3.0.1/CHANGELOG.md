# Changes from TISEAN 3.0.1

Changes are listed by program.
**Output** says whether the program's results differ from 3.0.1.
Each entry corresponds to one or more commits after the tag `upstream-3.0.1`, whose messages give the details.

The changes come from three sources, credited in each entry:

- **authors**: changes the TISEAN authors made after 3.0.1 but never released.
  Rainer Hegger published them in 2015 in github.com/heggus/Tisean (since deleted), with fixes from Piotr Held, author of the Octave port, in preparation for a TISEAN 3.0.2 that was never released.
  They are merged here with their original commits, from the copy in github.com/matzeri/Tisean, and tagged `authors-2015`.
- **hctsa**: fixes made for [hctsa](https://github.com/benfulcher/hctsa); those commits carry an `hctsa-commit:` trailer.
- **here**: fixes made in this repository, including two first made by Stefan Caldararu in [StefanCaldararu/TISEAN](https://github.com/StefanCaldararu/TISEAN).

## Unreleased

### c1 (and its subroutine d1.f)

- **Fixed a hang** when the number of embedded points is a power of two for some requested embedding dimension (e.g. 520 points, `-d2`, `-M5`): the mass loop could reach 1, which asked for one more neighbor than exists.
  Also skips embedding dimensions longer than the series, and stops searching after 200 box enlargements as a safety net.
  Output: unchanged wherever c1 previously finished.
  *(hctsa)*
- **Fixed out-of-range memory use** when shuffling the centers: the shuffle partner `int(rand*nmax-(m-1)*d)+1` could be zero or negative for m ≥ 2, so the result depended on whatever was in memory.
  The number of centers is also capped at the number of embedded points.
  Output: changes.
  *(hctsa)*
- **Deterministic centers**: the centers are the first points of a fixed, evenly spread ordering (a golden-ratio lattice) instead of a random permutation, so the result no longer depends on the seed `-I`.
  Output: changes.
  *(hctsa)*
- **Fixed a bias in the dimension estimates**: the sum of log radii was divided by `ncmin-(m-1)*d` instead of by the number of centers, which shrank every estimate, more so at large m·d and with few centers (on iid noise at d=10, m=4: 2.35 / 3.17 / 3.32 for 100 / 500 / 2000 centers, now 3.36 / 3.37 / 3.37).
  Output: changes.
  *(hctsa)*
- **Repeated points**: a center that coincides with its k-th neighbor (to within 1e-12 of the data range) makes the mean log distance minus infinity, so that mass is skipped, instead of flooring the distance at 1e-20.
  Output: changes only for data with repeated points.
  *(hctsa)*

### c2g, c2t

- **Double precision** throughout.
  In single precision the interpolation prefactor in c2g overflowed on steep local slopes, giving Inf/NaN.
- **c2g no longer integrates a stale point**: the counter was incremented before the check that ends each block, so an uninitialized point, or one left from the previous embedding dimension, entered the integral.
- Output: changes.
  Inf/NaN results become finite; otherwise results agree closely.
  *(hctsa)*

### corr

- New option `-s`: compute the autocorrelation only every `s`-th lag.
  Output: unchanged by default.
  *(authors)*

### d2

- Writes correlation sums and slopes with `%.17g` instead of `%e`, so they are not rounded to six significant digits.
  Output: same values, more digits.
  *(hctsa)*

### delay

- **Fixed the file output for multivariate data**: with `-o`, every component was embedded with the first component's delays (`-D`), so the file differed from the standard output.
  Output: changes for `-o` with several components and component-specific delays.
  *(authors)*

### false_nearest

- **Fixed a segfault** when `(M+1)·d` is at least the data length (unsigned loop bounds wrapped around); the program now exits with status 54 and a message.
  The authors had added a similar check, `M·d+1 ≥ N`, which misses cases such as N=100, d=19, M=5; it is replaced by this one.
  Output: unchanged wherever it previously ran.
  *(hctsa)*
- **The delay `-d` now sets the lag between embedding coordinates for a scalar series.**
  In 3.0.1 the scalar case always embedded consecutive samples, and `-d` only changed which points were used.
  The authors made the same fix in 2009 but never released it.
  Output: changes for `-d` > 1; identical for `-d1`.
  *(authors, hctsa)*
- **Ties**: distances within 1e-12 (data rescaled to [0,1]) count as equal, pairs closer than that are repeated points rather than neighbors, and every neighbor tied with the nearest is tested, each point contributing the fraction of them that are false.
  The result no longer depends on rounding or on the order of the box search, and is unchanged by affine transforms and negation of the series.
  Output: unchanged when nearest neighbors are unique; changes with ties.
  *(hctsa)*
- **A fifth output column**: the number of points behind each fraction.
  Scripts that expect exactly four columns need updating.
  *(hctsa)*

### fsle

- No longer reads one element past the end of the data when following a trajectory.
  Output: unchanged.
  *(here; first fixed by Stefan Caldararu)*

### histogram

- New option `-D`: output densities instead of relative frequencies.
  Output: unchanged by default.
  *(authors)*

### ikeda

- The initial condition is set with `-R` and `-I`, as the help text always said; 3.0.1 read `-X` and `-Y`.
  *(authors: Piotr Held)*

### lyap_k

- With `-M1`, the number of reference points was computed before the dimension was raised to its minimum of 2, so the neighbor search read past the end of the data.
  Output: changes for `-M1` (previously affected by out-of-range memory); unchanged otherwise.
  *(here; first fixed by Stefan Caldararu)*

### lyap_r

- **Fixed an endless loop** on data where some reference point has no neighbor within the acceptance radius (coarsely quantized series, or a few outliers squashing the rest into one box): epsilon grew until it overflowed.
  The search now stops, and if any reference point is unresolved lyap_r writes no output and exits with status 54, rather than reporting an estimate from only some of the points.
  Output: unchanged wherever lyap_r previously finished, except that such data now give exit status 54 and no output.
  *(hctsa)*
- **Ties**: every neighbor tied with the nearest (to within 1e-12, repeated points excluded) contributes with equal weight, instead of whichever was found first, and pairs that merge within the steps followed are left out at every step.
  Output: unchanged when nearest neighbors are unique and never merge; changes with ties.
  *(hctsa)*

### lyap_spec

- **New option `-t`**, a Theiler window: neighbors closer in time than this to the reference point are not used in the local linear fits.
  Output: unchanged by default (`-t0`).
  *(hctsa)*

### lzo-run, lzo-test

- lzo-run: fixing the number of neighbors (`-K`) is off by default, as documented; in 3.0.1 it was always on.
  Output: changes by default.
  *(authors)*
- lzo-test: new option `-K`, as in lzo-run.
  *(authors)*

### mem_spec

- **The file and standard output are normalized the same way**, as in `spectrum` (divided by N); in 3.0.1 the file output was divided by √N and the standard output not at all.
  Output: changes.
  *(authors)*
- No longer reads one element past the end of the data in the Burg recursion.
  Output: unchanged.
  *(here)*

### polynom

- With `-L`, no longer closes its output file twice (which aborts under glibc).
  *(authors: Piotr Held)*
- The in-sample and out-of-sample errors are given in units of the data.
  Output: changes.
  *(authors)*

### polynomp

- Passed a pointer, cast to `char`, where `get_multi_series` expects a 0/1 flag, so its behavior depended on the address.
  *(authors)*

### randomize_\*

- The final write called `writecfile` with the argument list of `xwritecfile`, so the final surrogate was not written (the output file kept the last intermediate one) and the program could segfault on exit.
  Output: changes (the final surrogate is written).
  *(authors; also fixed by Stefan Caldararu)*

### rank.f (used by c1)

- No longer reads `x(0)` at the end of a list (Fortran does not short-circuit `.and.`), which stopped bounds-checked builds.
  Output: unchanged.
  *(hctsa)*

### xcor

- Type fix in the lag loop.
  Output: unchanged.
  *(authors)*

### Build, tests and documentation

- `configure` finds C compilers again: its test program declared `main` without a return type, which clang 16+ and GCC 14+ reject, so no C programs were built.
  *(here)*
- `configure` tries `gfortran -std=legacy` and `gfortran`.
  *(authors, hctsa)*
- `cluster.f` compiles with gfortran (a format with no field width).
  *(here; the authors fixed it differently)*
- `make install` creates the install directory if needed, and no longer stops before installing the C programs when no Fortran compiler is found.
  *(here, hctsa)*
- Hollerith constants in format strings replaced by quoted strings, and two C warnings fixed.
  *(hctsa)*
- The `unsupported/` directory and three stray files in `source_c/routines/` are removed.
  *(authors)*
- Documentation updated for corr, lfo-ar, lyap_spec, lzo-test and the tutorial.
  *(authors)*
- Regression tests (`tests/`) and continuous integration on Linux and macOS, including a build with AddressSanitizer and Fortran bounds checking.
  *(here)*
