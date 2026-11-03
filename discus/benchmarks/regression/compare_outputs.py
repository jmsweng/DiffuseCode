#!/usr/bin/env python3
"""
Compare two DISCUS output files point by point, e.g. the same macro run with
the master build and with the perf/gemm-fourier build.

    python3 compare_outputs.py reference.h5 test.h5 [--tol 1e-5]
    python3 compare_outputs.py ref_dir/ test_dir/  [--tol 1e-5]   # all files with equal names

Formats
    HDF5 ('form common' / 'form hdf5' / 'form nexus'): every dataset whose name
         ends in data_values (or the first float dataset found) is compared.
    Text (gnu, xplor, vesta, ...): all numbers after the header are compared.
         Text formats round to ~5 significant digits, so use HDF5 when you can.

Reported per file:
    rel.max = max|a-b| / max|a|
    p99     = 99th percentile of |a-b|/|a| over points with |a| > 1e-6 max|a|
    corr    = Pearson correlation
Exit status 1 if any file exceeds --tol (default 1e-5 of the maximum).
"""
import argparse
import os
import re
import sys

import numpy as np

NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?')


def load(path):
    if path.endswith(('.h5', '.hdf5', '.nxs', '.nx')):
        import h5py
        found = {}

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and obj.dtype.kind == 'f' and obj.size > 1:
                found[name] = obj
        with h5py.File(path, 'r') as f:
            f.visititems(visit)
            keys = [k for k in found if k.endswith('data_values')] or list(found)[:1]
            if not keys:
                raise ValueError('no float dataset in %s' % path)
            return np.array(f[keys[0]], dtype=float).ravel()
    vals = []
    with open(path, errors='replace') as fh:
        for line in fh:
            t = line.split()
            if not t or not NUM.fullmatch(t[0].replace('D', 'E')):
                continue                     # header / title lines
            vals.extend(float(x.replace('D', 'E')) for x in t if NUM.fullmatch(x.replace('D', 'E')))
    return np.array(vals)


def compare(a, b):
    if a.shape != b.shape:
        return None
    m = np.max(np.abs(a))
    d = np.abs(a - b)
    sig = np.abs(a) > 1e-6 * m
    p99 = np.percentile(d[sig] / np.abs(a[sig]), 99) if sig.any() else 0.0
    corr = np.corrcoef(a, b)[0, 1] if a.size > 1 and a.std() > 0 and b.std() > 0 else 1.0
    return d.max() / m if m > 0 else d.max(), p99, corr


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ref')
    ap.add_argument('test')
    ap.add_argument('--tol', type=float, default=1e-5)
    a = ap.parse_args()
    if os.path.isdir(a.ref):
        names = sorted(f for f in os.listdir(a.ref) if os.path.isfile(os.path.join(a.ref, f))
                       and os.path.exists(os.path.join(a.test, f)) and not f.endswith(('.log', '.mac', '.stru', '.cll', '.cell')))
        pairs = [(os.path.join(a.ref, n), os.path.join(a.test, n)) for n in names]
    else:
        pairs = [(a.ref, a.test)]
    bad = 0
    print('%-40s %10s %10s %12s  %s' % ('file', 'rel.max', 'p99 rel', 'corr', 'result'))
    for fa, fb in pairs:
        try:
            r = compare(load(fa), load(fb))
        except Exception as e:                               # unreadable / not numeric
            print('%-40s  skipped (%s)' % (os.path.basename(fa), e))
            continue
        if r is None:
            print('%-40s  SIZE MISMATCH' % os.path.basename(fa))
            bad += 1
            continue
        ok = r[0] <= a.tol
        bad += 0 if ok else 1
        print('%-40s %10.2e %10.2e %12.9f  %s' % (os.path.basename(fa), r[0], r[1], r[2], 'ok' if ok else 'DIFFERS'))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
