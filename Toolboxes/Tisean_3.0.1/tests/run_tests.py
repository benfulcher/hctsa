#!/usr/bin/env python3
"""Regression tests for the TISEAN programs used by hctsa.

Each case runs one or more TISEAN command lines in a fresh temporary directory
holding a copy of the test data, then compares the exit status, standard output
and every file the commands wrote with the reference copies in tests/expected/.
Numbers are compared with a relative tolerance (output formats and the last
digits of floating-point results can differ between compilers and platforms);
anything else must match exactly.

Usage (from the repository root, after ./configure && make):
    python3 tests/run_tests.py              # run all cases
    python3 tests/run_tests.py d2 c1        # run the cases whose name starts with these
    python3 tests/run_tests.py --update     # rewrite the reference outputs
    python3 tests/run_tests.py --bindir DIR # use installed binaries in DIR
"""

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
EXPECTED = ROOT / "tests" / "expected"

RTOL = 1e-6
ATOL = 1e-12
# c1 works in single precision, so its results differ more between platforms
RTOL_SINGLE = 1e-5

# name: (list of command lines, expected exit status of the last command,
#        timeout in seconds, relative tolerance)
# Command lines for the programs hctsa uses follow the way hctsa calls them
# (hctsa/Operations).
CASES = {
    # d2 and the programs that read its output (NL_d2, NL_GPCorrSum)
    "d2_henon": (["d2 -d1 -M1,4 -t10 -V0 henon.dat"], 0, 60, RTOL),
    "d2_c2g_c2t_henon": (["d2 -d1 -M1,4 -t10 -V0 henon.dat",
                          "c2g -V0 henon.dat.c2",
                          "c2t -V0 henon.dat.c2"], 0, 60, RTOL),
    "d2_gpcorrsum_lorenz": (["d2 -d5 -M1,3 -t20 -R0.5 -N1000 -#20 -V0 lorenz.dat"],
                            0, 60, RTOL),
    # c1 and c2d (NL_c1)
    "c1_c2d_henon": (["c1 -d1 -m1 -M4 -t10 -n500 -V0 -o henon.c1 henon.dat",
                      "c2d -a2 -V0 henon.c1"], 0, 60, RTOL_SINGLE),
    # used to hang: 520 points, delay 2, m = 5 gives 512 embedded points
    "c1_powerof2_nohang": (["c1 -d2 -m1 -M5 -t1 -n100 -l520 -V0 -o h520.c1 henon.dat"],
                           0, 60, RTOL_SINGLE),
    # false_nearest (NL_FNN)
    "fnn_lorenz_delay5": (["false_nearest -d5 -m1 -M1,6 -t10 -V0 lorenz.dat"], 0, 60, RTOL),
    "fnn_henon_ties": (["false_nearest -d1 -m1 -M1,5 -t1 -V0 henon_q1.dat"], 0, 60, RTOL),
    # used to segfault: (M+1)*d >= length
    "fnn_too_short": (["false_nearest -d19 -m1 -M1,5 -t1 -l100 -V0 henon.dat"], 54, 60, RTOL),
    # lyap_r (NL_LargestLyap)
    "lyap_r_lorenz": (["lyap_r -d5 -m3 -t20 -s30 -V0 -o lorenz.lyap lorenz.dat"], 0, 60, RTOL),
    "lyap_r_henon_ties": (["lyap_r -d1 -m2 -t1 -s10 -V0 -o henon_q1.lyap henon_q1.dat"],
                          0, 60, RTOL),
    # used to loop for ever: no neighbors within the acceptance radius
    "lyap_r_not_estimable": (["lyap_r -d1 -m3 -t1 -s20 -V0 -o binary.lyap binary.dat"],
                             54, 60, RTOL),
    # lyap_spec (NL_LyapSpec), with and without a Theiler window
    "lyap_spec_henon": (["lyap_spec -m1,2 -k30 -V0 -o henon.lyaps henon.dat"], 0, 120, 1e-4),
    "lyap_spec_henon_theiler": (["lyap_spec -m1,2 -k30 -t5 -V0 -o henon.lyaps henon.dat"],
                                0, 120, 1e-4),
    # boxcount (NL_BoxCountEntropyRate, NL_Dimensions)
    "boxcount_henon": (["boxcount -M1,4 -d1 -Q2.0 -#20 -V0 -o henon.box henon.dat"],
                       0, 60, RTOL),
    # nstat_z (SY_nstat_z)
    "nstat_z_lorenz": (["nstat_z -# 5 -d5 -m3 -V0 lorenz.dat"], 0, 60, RTOL),
    # nrlazy (PP_SchreiberDenoise)
    "nrlazy_lorenz": (["nrlazy -m1,5 -d1 -i2 -v0.5 -V0 -o lorenz.laz lorenz.dat"], 0, 60, RTOL),

    # Programs hctsa does not use, for bugs fixed since 3.0.1
    # delay: the file output (-o) used the first component's delays for every
    # component, so it differed from the stdout output
    "delay_multivariate": (["delay -M2 -F2,2 -D1,3 -V0 henon2.dat",
                            "delay -M2 -F2,2 -D1,3 -V0 -o henon2.del henon2.dat"], 0, 60, RTOL),
    # mem_spec: the file and stdout outputs were normalized differently
    "mem_spec_henon": (["mem_spec -p20 -P10 -V0 henon.dat",
                        "mem_spec -p20 -P10 -V0 -o henon.spec henon.dat"], 0, 60, RTOL),
    # polynom -L: closed its output file twice (aborts under glibc)
    "polynom_cast": (["polynom -m2 -p2 -L20 -V0 -o henon.pol henon.dat"], 0, 60, RTOL),
    # randomize: the final write called writecfile with the wrong arguments
    "randomize_auto": (["randomize_auto_exp_random -n1 -I1 -D5 -l200 -V0 -a0.5 -S200 -s50 -z40 "
                        "-o henon.rnd henon.dat"], 0, 60, RTOL),
    # fsle and lyap_k -M1: read past the end of the data (seen with AddressSanitizer)
    "fsle_henon": (["fsle -d1 -m2 -t0 -V0 -o henon.fsle henon.dat"], 0, 60, RTOL),
    "lyap_k_m1": (["lyap_k -M1 -m1 -d1 -s20 -V0 -o henon.lyk henon.dat"], 0, 60, RTOL),
}


def find_binaries(bindir):
    """Directories to search for the TISEAN programs, in order."""
    if bindir:
        return [Path(bindir).resolve()]
    return [ROOT / "source_c", ROOT / "source_f", ROOT / "source_f" / "randomize"]


def resolve(cmd, dirs):
    prog, _, rest = cmd.partition(" ")
    for d in dirs:
        exe = d / prog
        if exe.is_file() and os.access(exe, os.X_OK):
            return [str(exe)] + rest.split()
    raise FileNotFoundError(f"{prog} not found in {', '.join(map(str, dirs))} (run make?)")


def run_case(name, dirs):
    """Run a case; return (exit status of the last command, {output name: text})."""
    cmds, _, timeout, _ = CASES[name]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for f in DATA.iterdir():
            shutil.copy(f, tmp / f.name)
        inputs = {f.name for f in DATA.iterdir()}
        outputs = {}
        status = None
        for i, cmd in enumerate(cmds):
            try:
                proc = subprocess.run(resolve(cmd, dirs), cwd=tmp, capture_output=True,
                                      text=True, errors="replace", timeout=timeout)
            except subprocess.TimeoutExpired:
                return "timeout", {}
            status = proc.returncode
            if proc.stdout:
                outputs[f"stdout.{i + 1}"] = proc.stdout
        for f in sorted(tmp.iterdir()):
            if f.name not in inputs:
                outputs[f.name] = f.read_text(errors="replace")
        return status, outputs


def as_number(tok):
    try:
        return float(tok.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def compare(got, want, rtol):
    """Return a description of the first difference, or None."""
    gl, wl = got.splitlines(), want.splitlines()
    if len(gl) != len(wl):
        return f"{len(gl)} lines, expected {len(wl)}"
    for n, (g, w) in enumerate(zip(gl, wl), 1):
        gt, wt = g.split(), w.split()
        if len(gt) != len(wt):
            return f"line {n}: {g!r}, expected {w!r}"
        for a, b in zip(gt, wt):
            x, y = as_number(a), as_number(b)
            if x is None or y is None:
                if a != b:
                    return f"line {n}: {g!r}, expected {w!r}"
            elif not (x == y or (math.isnan(x) and math.isnan(y))
                      or abs(x - y) <= ATOL + rtol * max(abs(x), abs(y))):
                return f"line {n}: {a}, expected {b}"
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cases", nargs="*", help="run only cases whose name starts with these")
    ap.add_argument("--bindir", help="directory holding the TISEAN programs")
    ap.add_argument("--update", action="store_true", help="rewrite the reference outputs")
    args = ap.parse_args()

    dirs = find_binaries(args.bindir)
    names = [n for n in CASES if not args.cases or any(n.startswith(c) for c in args.cases)]
    failed = []
    for name in names:
        want_status, rtol = CASES[name][1], CASES[name][3]
        status, outputs = run_case(name, dirs)
        problems = []
        if status != want_status:
            problems.append(f"exit status {status}, expected {want_status}")
        ref = EXPECTED / name
        if args.update and not problems:
            shutil.rmtree(ref, ignore_errors=True)
            ref.mkdir(parents=True)
            for out, text in outputs.items():
                (ref / out).write_text(text)
        elif not problems:
            want = {f.name: f.read_text() for f in sorted(ref.iterdir())} if ref.is_dir() else {}
            if set(outputs) != set(want):
                problems.append(f"wrote {sorted(outputs)}, expected {sorted(want)}")
            for out in sorted(set(outputs) & set(want)):
                diff = compare(outputs[out], want[out], rtol)
                if diff:
                    problems.append(f"{out}: {diff}")
        print(f"{'FAIL' if problems else 'ok  '} {name}")
        for p in problems:
            print(f"       {p}")
        if problems:
            failed.append(name)
    print(f"\n{len(names) - len(failed)} of {len(names)} cases passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
