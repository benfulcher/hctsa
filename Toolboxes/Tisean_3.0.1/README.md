# TISEAN, maintained

[TISEAN](https://www.pks.mpg.de/tisean/) is a package of programs for nonlinear time-series analysis by Rainer Hegger, Holger Kantz and Thomas Schreiber.
Its last release, 3.0.1, dates from 2007.
This repository is TISEAN 3.0.1 with bug fixes and build fixes, maintained by the [Dynamics and Neural Systems group](https://github.com/DynamicsAndNeuralSystems) at the University of Sydney.
It includes the authors' own fixes from 2007–2015, which were prepared for a TISEAN 3.0.2 that was never released, and further fixes developed for [hctsa](https://github.com/benfulcher/hctsa), which runs eleven of the TISEAN programs on large numbers of time series.

This is not an official TISEAN release.
The original is available from the [TISEAN website](https://www.pks.mpg.de/tisean/).

## What is different from 3.0.1

The first commit is the 3.0.1 tarball exactly as released, tagged `upstream-3.0.1`, and every later change is a separate commit.
`git diff upstream-3.0.1` shows everything that differs.
The authors' unreleased changes are merged with their original commits (Rainer Hegger's, with fixes from Piotr Held), and tagged `authors-2015`.
[CHANGELOG.md](CHANGELOG.md) lists the changes by program, and says for each whether it changes the program's output.

In brief:

- **Crashes and hangs**: `false_nearest` no longer segfaults when the embedding is too long for the data, `c1` no longer hangs for some data lengths, `lyap_r` no longer loops for ever on quantized data, `polynom -L` and `randomize_*` no longer crash, and `fsle`, `lyap_k` and `mem_spec` no longer read past the end of the data.
- **Results**: `false_nearest` now uses the delay `-d` for a scalar series (3.0.1 ignored it), `c1` no longer biases dimension estimates low, `c2g` no longer integrates a stale point or overflows, and `false_nearest`, `lyap_r` and `c1` treat tied distances and repeated points independently of rounding and search order.
  `delay -o` now embeds multivariate data correctly, and `mem_spec` normalizes its file and standard output the same way.
- **Build**: `./configure && make` works with current compilers (clang, GCC 14, gfortran), and all 80 programs build.

Of the programs hctsa uses, those whose output can differ from 3.0.1 are `c1`, `c2g`, `c2t`, `false_nearest` and `lyap_r`; of the others, `delay` (with `-o`), `lyap_k` (with `-M1`), `lzo-run`, `mem_spec`, `polynom` and `randomize_*` (see the changelog for when).
To reproduce published results computed with the original, build the `upstream-3.0.1` tag.

## Scope

All 80 programs build, but only those used by hctsa are tested and actively maintained:
`boxcount`, `c1`, `c2d`, `c2g`, `c2t`, `d2`, `false_nearest`, `lyap_r`, `lyap_spec`, `nrlazy` and `nstat_z`.
The others have the authors' fixes and those listed in the changelog, but are otherwise as released in 3.0.1.
Bug reports and fixes for any program are welcome.

## Installing

You need a C compiler, a Fortran compiler (gfortran) and `make`.
On macOS, `brew install gcc` provides gfortran.

```sh
./configure --prefix=$HOME/.local   # programs go in $HOME/.local/bin
make
make install
```

Windows is not supported directly; use WSL.

## Testing

```sh
python3 tests/run_tests.py                       # after make, using the programs in source_c/ and source_f/
python3 tests/run_tests.py --bindir ~/.local/bin # or the installed programs
```

The tests run the programs hctsa uses with hctsa's command lines, and six other programs on the bugs fixed in them, and compare the output with reference output in `tests/expected/` to a relative tolerance of 1e-6 (1e-5 for the single-precision c1).
Fourteen of the 22 cases fail on TISEAN 3.0.1; the fsle case catches its out-of-range read only in a build with AddressSanitizer, as in continuous integration.

## Documentation

The original HTML documentation is in [docs/](docs/) (open `index.html`), and also on the [TISEAN website](https://www.pks.mpg.de/tisean/Tisean_3.0.1/index.html).
Options added here are listed by each program's `-h`.

## Citing

If you use TISEAN, please cite:

R. Hegger, H. Kantz and T. Schreiber, *Practical implementation of nonlinear time series methods: The TISEAN package*, Chaos **9**, 413 (1999).
[doi:10.1063/1.166424](https://doi.org/10.1063/1.166424)

and, for the surrogate data programs, T. Schreiber and A. Schmitz, *Surrogate time series*, Physica D **142**, 346 (2000).

## License

TISEAN is licensed under the GNU General Public License, version 2 or later (see [COPYING](COPYING)).
The copyright to each program is held by the authors named in its source code.
Changes in this repository are under the same license.
