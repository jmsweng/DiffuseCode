#!/usr/bin/env python3
"""
Full-volume regression tests for the DISCUS Fourier / 3D-PDF changes on
branch perf/gemm-fourier.

Runs the same DISCUS macros with a reference binary (built from master)
and a test binary (built from this branch), writes every result as an
HDF5 'common' file (full 3D volume, double precision) and compares all
points.

Usage
    python3 regress.py --ref  /path/to/master/discus_suite \
                       --test /path/to/branch/discus_suite \
                       [--work DIR] [--threads N] [--only PATTERN] [--keep]

Needs numpy and h5py. Exit status 0 if every check passes.

Check kinds
    exact   test with DISCUS_FOURIER_KERNEL=table vs ref       rel. max diff <= 1e-12
            (isolates every change except the GEMM kernel; must match bit for bit
             up to summation order)
    kernel  test default kernel vs ref                         rel. max diff <= 1e-2
            Gross-error guard only. With -DDISCUS_BLAS=ON the test binary uses the
            GEMM kernel (exact phases); ref uses the 2^16-entry sin/cos table, whose
            error reaches ~1e-3 of the maximum (with 'set aver') and ~1% at single
            weak diffuse points. Accuracy is judged by the oracle checks instead.
    oracle  test NUFFT vs test turbo (same binary)             rel. max diff <= 1e-5
            Two independent exact-phase methods (GEMM or table turbo vs FINUFFT,
            tol 1e-7); with the table kernel compiled in, expect ~1e-4 here.
    mask    NUFFT with 'set aver' vs plain NUFFT               must be 0 at integer
            hkl and identical elsewhere
    info    reported only (e.g. known-wrong behaviour of the old NUFFT)

"rel. max diff" = max|a-b| / max|ref|. Also reported: 99th percentile of the
pointwise relative difference over points above 1e-6 * max (the weak diffuse
signal, where a 3D-dPDF is most sensitive).
"""
import argparse
import fnmatch
import os
import re
import subprocess
import sys
import time

import numpy as np
import h5py

# ----------------------------------------------------------------------------
# Crystals (built once with the test binary, then read by both binaries)
# ----------------------------------------------------------------------------
CELL_CU = """title Cu test cell
spcgr P1
cell  3.6, 3.6, 3.6, 90.0, 90.0, 90.0
atoms
CU   0.00 0.00 0.00 0.5
CU   0.50 0.50 0.00 0.5
CU   0.50 0.00 0.50 0.5
CU   0.00 0.50 0.50 0.5
"""
CELL_CU2 = """title Cu test cell, two ADP types
spcgr P1
cell  3.6, 3.6, 3.6, 90.0, 90.0, 90.0
atoms
CU   0.00 0.00 0.00 0.5
CU   0.50 0.50 0.00 0.5
CU   0.50 0.00 0.50 0.8
CU   0.00 0.50 0.50 0.8
"""
CRYSTALS = {
    # name : (cell file contents, nx, ny, nz)
    'x1':   (CELL_CU, 16, 16, 16),     # 16384 atoms, Cu/Au random alloy, isotropic, 1 ADP type
    'x2':   (CELL_CU2, 12, 12, 12),    # 6912 atoms, 2 ADP types
    'thin': (CELL_CU, 16, 16, 1),      # one unit cell thick along c
    'x3':   (CELL_CU, 8, 8, 8),        # 2048 atoms, for lots over all cells
    'x4':   (CELL_CU, 12, 12, 12),     # 6912 atoms, random thermal displacements ('ther')
}
CRYSTAL_POST = {'x4': ['ther']}
ANIS = ['anis type:1, values:[0.010, 0.020, 0.005, 0.0, 0.0, 0.002]',
        'anis type:2, values:[0.015, 0.008, 0.012, 0.001, 0.0, 0.0]']
ANIS_CHANGED = ['anis type:1, values:[0.030, 0.010, 0.010, 0.0, 0.003, 0.0]']

# ----------------------------------------------------------------------------
# Reciprocal-space grids:  (ll, lr, ul, tl, nabs, nord, ntop)
# ----------------------------------------------------------------------------
G = {
    'odd3d':   ('-4,-4,-4', '4,-4,-4', '-4,4,-4', '-4,-4,4', 41, 41, 41),   # step 0.2, symmetric -> Friedel
    'even3d':  ('-3.9,-3.5,-3.1', '3.9,-3.5,-3.1', '-3.9,3.5,-3.1', '-3.9,-3.5,3.1', 40, 36, 32),  # no Friedel, 3 sizes
    'offc3d':  ('0.1,0.2,0.3', '3.1,3.2,0.3', '-1.9,2.2,0.3', '0.1,0.2,2.8', 31, 25, 21),  # off-centre, skewed axes
    'aver3d':  ('-3,-3,-3', '3,-3,-3', '-3,3,-3', '-3,-3,3', 31, 31, 31),   # step 0.2, Bragg on grid
    'layer':   ('-4,-4,0.5', '4,-4,0.5', '-4,4,0.5', '-4,-4,0.5', 81, 81, 1),
    'skewlay': ('-3,-3,0.5', '3,3,0.5', '-3,3,0.5', '-3,-3,0.5', 121, 61, 1),
    'line':    ('-6,0.5,0', '6,0.5,0', '-6,0.5,0', '-6,0.5,0', 241, 1, 1),
    'line_ap': ('-6,0,0', '6,0,0', '-6,0,0', '-6,0,0', 241, 1, 1),        # along a*, through origin
    'thin3d':  ('-3,-3,-3', '3,-3,-3', '-3,3,-3', '-3,-3,3', 31, 31, 31),
    'small3d': ('-2.5,-2.5,-2.5', '2.5,-2.5,-2.5', '-2.5,2.5,-2.5', '-2.5,-2.5,2.5', 21, 21, 21),
    'aver_ok': ('-3,-3,-3', '3,-3,-3', '-3,3,-3', '-3,-3,3', 25, 25, 25),  # step 0.25: 16 cells * 0.25 = 4
}


def fourier_block(grid, technique='turbo', radiation='xray', wvle='1.5406', temp='ignore',
                  lots='off', disp=None, aver=None, extra=()):
    ll, lr, ul, tl, na, no, nt = G[grid]
    lines = ['fourier', radiation]
    if radiation != 'neut':
        lines.append('wvle ' + wvle)
    if disp:
        lines.append('disp ' + disp)
    lines += ['temp ' + temp, 'lots ' + lots,
              'll ' + ll, 'lr ' + lr, 'ul ' + ul, 'tl ' + tl,
              'nabs %d' % na, 'nord %d' % no, 'ntop %d' % nt,
              'set technique:' + technique]
    if aver is not None:
        lines.append('set aver, %g' % aver)
    lines += list(extra)
    lines += ['run', 'exit']
    return lines


def output_block(fname, value='intensity'):
    return ['output', 'outfile ' + fname, 'form common', 'value ' + value, 'run', 'exit']


def session(crystal, steps, pre=()):
    """steps: list of (fourier_lines, [(outfile, value), ...]) run in ONE DISCUS session."""
    lines = ['discus', 'seed 4711', 'read', 'stru %s.stru' % crystal] + list(pre)
    for f_lines, outs in steps:
        lines += f_lines
        for fname, value in outs:
            lines += output_block(fname, value)
    lines += ['exit', 'exit']
    return '\n'.join(lines) + '\n'


# ----------------------------------------------------------------------------
# Cases.  Each case: name -> dict(macro=..., checks=[...])
#   A run is a (binary, environment, threads) combination, see RUNS in main().
#   check = (kind, run_a, file_a, run_b, file_b); run_a/file_a is the reference.
# ----------------------------------------------------------------------------
def build_cases():
    C = {}

    def simple(name, crystal, grid, value='intensity', pre=(), kinds=('exact', 'kernel'), **kw):
        out = name + '.h5'
        checks = []
        if 'exact' in kinds:
            checks.append(('exact', 'ref', out, 'test_table', out))
        if 'kernel' in kinds:
            checks.append(('kernel', 'ref', out, 'test', out))
        C[name] = dict(macro=session(crystal, [(fourier_block(grid, **kw), [(out, value)])], pre=pre),
                       checks=checks)

    # --- turbo kernel / form factor paths ---------------------------------
    simple('t_odd3d',     'x1', 'odd3d')                                    # Friedel half + mirror
    simple('t_even3d',    'x1', 'even3d')                                   # no Friedel, unequal sizes
    simple('t_offc3d',    'x1', 'offc3d')                                   # lmn/off_shift, skewed axes
    simple('t_layer',     'x1', 'layer')
    simple('t_skewlay',   'x1', 'skewlay')
    simple('t_line',      'x1', 'line')                                     # 1D: original kernel in both
    simple('t_anom',      'x1', 'odd3d', wvle='CUA1', disp='anom')          # anomalous -> Friedel off
    simple('t_neut',      'x1', 'odd3d', radiation='neut')
    simple('t_iso2',      'x2', 'odd3d', temp='use')                        # 2 isotropic ADP types
    simple('t_aniso2',    'x2', 'odd3d', temp='use', pre=ANIS)              # anisotropic, 2 types
    simple('t_aniso_off', 'x2', 'offc3d', temp='use', pre=ANIS)
    simple('t_aver',      'x1', 'aver3d', aver=100)                         # four_aver path
    # lots: 'all' places a lot at every unit cell, which is deterministic. Random lots are
    # not reproducible even in master ('seed' does not reach the lot placement).
    simple('t_lots_p',    'x3', 'small3d', lots='box, 3, 3, 3, all, yes')
    simple('t_lots_np',   'x3', 'small3d', lots='box, 3, 3, 3, all, no')
    simple('t_lots_ell',  'x3', 'small3d', lots='ell, 4, 4, 4, all, yes', aver=100)

    # displaced atoms: arbitrary phases (perfect-lattice crystals hit exact table entries)
    simple('t_disp',      'x4', 'odd3d')
    simple('t_disp_off',  'x4', 'offc3d', temp='use')
    simple('t_disp_lots', 'x4', 'small3d', lots='box, 3, 3, 3, all, yes')
    simple('t_disp_aver', 'x4', 'aver3d', aver=100)

    # --- 3D-PDF (r2c FFT) --------------------------------------------------
    simple('p_odd3d',  'x1', 'odd3d',  value='3DPDF')
    simple('p_even3d', 'x1', 'even3d', value='3DPDF')                       # even sizes, 3 different
    simple('p_offc3d', 'x1', 'offc3d', value='3DPDF')
    simple('p_layer',  'x1', 'layer',  value='3DPDF')                       # 2D path (unchanged code)
    simple('p_aver',   'x1', 'aver3d', value='3DPDF', aver=100)             # 3D-dPDF proper
    simple('p_disp',   'x4', 'aver3d', value='3DPDF', aver=100)             # 3D-dPDF, displaced atoms

    # --- serial kernel (1 thread, original table kernel in both) ------------
    for cr, g in (('x1', 'odd3d'), ('x1', 'even3d'), ('x4', 'offc3d')):
        n = 's_' + g + ('_disp' if cr == 'x4' else '')
        out = n + '.h5'
        C[n] = dict(macro=session(cr, [(fourier_block(g), [(out, 'intensity')])]),
                    checks=[('exact', 'ref_1t', out, 'test_table_1t', out)])

    # --- Debye-Waller cache: one session with many changes -------------------
    #     compared with the same session in master AND with fresh one-step sessions
    seq = [('d1', 'odd3d',  'use',    ANIS),             # first table
           ('d2', 'odd3d',  'use',    ()),               # unchanged -> reuse
           ('d3', 'even3d', 'use',    ()),               # grid size changes
           ('d4', 'offc3d', 'use',    ()),               # grid position / axes change
           ('d5', 'offc3d', 'ignore', ()),               # temp off
           ('d6', 'offc3d', 'use',    ()),               # temp on again, same grid
           ('d7', 'offc3d', 'use',    ANIS_CHANGED)]     # ADPs change
    steps, checks, current = [], [], []
    for tag, g, t, anis in seq:
        if anis is ANIS:
            current = list(ANIS)
        elif anis:
            current = list(current) + list(anis)
        steps.append((list(anis) + fourier_block(g, temp=t), [(tag + '.h5', 'intensity')]))
        C['dbw_fresh_' + tag] = dict(
            macro=session('x2', [(fourier_block(g, temp=t), [(tag + '.h5', 'intensity')])], pre=current),
            checks=[], hidden=True)
        checks += [('exact', 'ref',        tag + '.h5', 'test_table', tag + '.h5'),
                   ('exact', 'ref_fresh',  tag + '.h5', 'test_table', tag + '.h5'),
                   ('info',  'ref_fresh',  tag + '.h5', 'ref',        tag + '.h5')]
    C['dbw_session'] = dict(macro=session('x2', steps), checks=checks)

    # --- NUFFT vs turbo (test binary), and vs master -----------------------
    def nufft(name, crystal, grid, pre=(), info_ref=False, **kw):
        tn, tt = name + '_nufft.h5', name + '_turbo.h5'
        mac = session(crystal, [(fourier_block(grid, technique='nufft', **kw), [(tn, 'intensity')]),
                                (fourier_block(grid, technique='turbo', **kw), [(tt, 'intensity')])],
                      pre=pre)
        checks = [('oracle', 'test', tt, 'test', tn)]
        if info_ref:
            checks.append(('info', 'test', tt, 'ref', tn))      # old NUFFT vs new turbo
        C[name] = dict(macro=mac, checks=checks)

    nufft('n_odd3d',   'x1',   'odd3d', info_ref=True)     # step*cells = 3.2: old NUFFT was wrong
    nufft('n_even3d',  'x1',   'even3d')
    nufft('n_offc3d',  'x1',   'offc3d')
    nufft('n_layer',   'x1',   'layer')
    nufft('n_skewlay', 'x1',   'skewlay')
    nufft('n_line',    'x1',   'line_ap')
    nufft('n_thin',    'thin', 'thin3d')                   # cr_icc(3) = 1
    nufft('n_aniso2',  'x2',   'odd3d', pre=ANIS, temp='use')
    nufft('n_neut',    'x1',   'odd3d', radiation='neut')
    nufft('n_anom',    'x1',   'odd3d', wvle='CUA1', disp='anom')
    nufft('n_disp',    'x4',   'odd3d')
    nufft('n_disp_off','x4',   'offc3d', temp='use')

    # NUFFT with average subtraction. The NUFFT path removes the average structure by
    # zeroing the integer hkl points, so 'aver' must equal the plain NUFFT result with
    # exactly those points set to zero ('mask' check).
    #   aver_ok: step*cells = 4, the old NUFFT was correct here -> old and new must agree
    #   aver3d : step*cells = 3.2, the old NUFFT was wrong here -> reported only
    for g, kind in (('aver_ok', 'exact'), ('aver3d', 'info')):
        n = 'n_' + g
        fa, fp = n + '_aver.h5', n + '_plain.h5'
        C[n] = dict(macro=session('x1', [(fourier_block(g, technique='nufft', aver=100), [(fa, 'intensity')]),
                                         (fourier_block(g, technique='nufft'),           [(fp, 'intensity')])]),
                    checks=[(kind, 'ref', fa, 'test', fa),
                            ('mask', 'test', fp, 'test', fa)])
    return C


TOL = {'exact': 1e-12, 'kernel': 1e-2, 'oracle': 1e-5, 'mask': 1e-12, 'info': None}

ERR_RE = re.compile(r'\*\*\*[A-Z]{3,4} ?\*\*\*')


def run_discus(binary, macro_path, cwd, env_extra, threads, timeout=3600):
    env = dict(os.environ)
    env.update({'OMP_NUM_THREADS': str(threads), 'OPENBLAS_NUM_THREADS': str(threads),
                'MKL_NUM_THREADS': str(threads)})
    env.pop('DISCUS_FOURIER_KERNEL', None)
    env.update(env_extra)
    log = macro_path[:-4] + '.log'
    t0 = time.time()
    with open(log, 'w') as fh:
        subprocess.run([binary, '-macro', os.path.basename(macro_path)], cwd=cwd, env=env,
                       stdin=subprocess.DEVNULL, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)
    dt = time.time() - t0
    errs = [l.strip() for l in open(log, errors='replace') if ERR_RE.search(l)]
    return dt, errs, log


def load(path):
    with h5py.File(path, 'r') as f:
        return np.array(f['entry/data/data_values'], dtype=float)


def load_hkl(path):
    """hkl of every grid point (corner + index . increment vectors)."""
    with h5py.File(path, 'r') as f:
        g = f['entry/data']
        shape = g['data_values'].shape
        corner = np.array(g['data_corner'], dtype=float)
        inc = np.array(g['data_increment_vector'], dtype=float)
    idx = np.indices(shape).reshape(3, -1).T
    return (corner + idx @ inc).reshape(shape + (3,))


def compare_mask(plain, aver, hkl):
    """aver must be 0 at integer hkl and equal to plain elsewhere."""
    bragg = np.all(np.abs(hkl - np.round(hkl)) < 1e-6, axis=-1)
    m = np.max(np.abs(plain))
    rel = max(np.max(np.abs(aver[bragg])) / m if bragg.any() else 0.0,
              np.max(np.abs(plain[~bragg] - aver[~bragg])) / m)
    return dict(rel=rel, p99=0.0, n=int(bragg.sum()))


def compare(a, b):
    """a = reference."""
    if a.shape != b.shape:
        return dict(shape=(a.shape, b.shape))
    m = np.max(np.abs(a))
    d = np.abs(a - b)
    rel = d.max() / m if m > 0 else d.max()
    sig = np.abs(a) > 1e-6 * m
    p99 = np.percentile(d[sig] / np.abs(a[sig]), 99) if sig.any() else 0.0
    return dict(rel=rel, p99=p99, n=a.size)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True, help='discus_suite built from master')
    ap.add_argument('--test', required=True, help='discus_suite built from perf/gemm-fourier')
    ap.add_argument('--work', default='regress_work')
    ap.add_argument('--threads', type=int, default=os.cpu_count() or 2)
    ap.add_argument('--only', default='*', help='glob on case names')
    ap.add_argument('--keep', action='store_true', help='reuse existing outputs in --work')
    a = ap.parse_args()

    work = os.path.abspath(a.work)
    os.makedirs(work, exist_ok=True)
    ref, test = os.path.abspath(a.ref), os.path.abspath(a.test)

    RUNS = {
        'ref':           (ref,  {}, a.threads),
        'ref_1t':        (ref,  {}, 1),
        'ref_fresh':     (ref,  {}, a.threads),
        'test':          (test, {}, a.threads),
        'test_table':    (test, {'DISCUS_FOURIER_KERNEL': 'table'}, a.threads),
        'test_table_1t': (test, {'DISCUS_FOURIER_KERNEL': 'table'}, 1),
    }

    # crystals
    for name, (cell, nx, ny, nz) in CRYSTALS.items():
        if a.keep and os.path.exists(os.path.join(work, name + '.stru')):
            continue
        open(os.path.join(work, name + '.cll'), 'w').write(cell)
        mac = os.path.join(work, 'make_%s.mac' % name)
        post = ''.join(l + '\n' for l in CRYSTAL_POST.get(name, []))
        open(mac, 'w').write('discus\nseed 1234\nread\ncell %s.cll, %d, %d, %d\n'
                             'replace CU, AU, all, 0.3\n%ssave %s.stru\nexit\nexit\n'
                             % (name, nx, ny, nz, post, name))
        _, errs, log = run_discus(test, mac, work, {}, 1)
        if errs or not os.path.exists(os.path.join(work, name + '.stru')):
            sys.exit('crystal %s failed, see %s' % (name, log))

    cases = build_cases()
    selected = [n for n in cases if fnmatch.fnmatch(n, a.only)]
    if any(n == 'dbw_session' for n in selected):
        selected += [n for n in cases if n.startswith('dbw_fresh_') and n not in selected]
    failed = 0

    for name in selected:
        c = cases[name]
        need = {r for chk in c['checks'] for r in (chk[1], chk[3])}
        if c.get('hidden'):
            need = {'ref_fresh'}
        for run in sorted(need):
            d = os.path.join(work, run)
            os.makedirs(d, exist_ok=True)
            for f in os.listdir(work):
                if f.endswith(('.stru', '.cll')) and not os.path.exists(os.path.join(d, f)):
                    os.symlink(os.path.join(work, f), os.path.join(d, f))
            mac = os.path.join(d, name + '.mac')
            open(mac, 'w').write(c['macro'])
            if a.keep and os.path.exists(mac[:-4] + '.log'):
                continue
            binary, env, thr = RUNS[run]
            dt, errs, log = run_discus(binary, mac, d, env, thr)
            c.setdefault('time', {})[run] = dt
            if errs:
                print('  %-16s %-14s DISCUS ERRORS (%s):' % (name, run, log))
                for e in errs[:4]:
                    print('      ' + e)
                failed += 1

    # Without -DDISCUS_BLAS=ON the test binary's turbo IS the table kernel: then the
    # NUFFT oracle can only be as good as the table (~1e-4 .. 1e-3 of the maximum).
    fa, fb = os.path.join(work, 'test', 't_odd3d.h5'), os.path.join(work, 'test_table', 't_odd3d.h5')
    if os.path.exists(fa) and os.path.exists(fb):
        if np.array_equal(load(fa), load(fb)):
            TOL['oracle'] = 1e-3
            print('NOTE: test binary has no GEMM kernel (built without -DDISCUS_BLAS=ON); '
                  'oracle tolerance relaxed to 1e-3, kernel checks are identical to exact checks.')
        else:
            print('NOTE: test binary uses the GEMM kernel (-DDISCUS_BLAS=ON).')

    print('\n%-14s %-7s %-42s %10s %10s  %s' % ('case', 'kind', 'compare (reference vs test)',
                                                'rel.max', 'p99 rel', 'result'))
    for name in selected:
        for kind, ra, fa_, rb, fb_ in cases[name]['checks']:
            fa, fb = os.path.join(work, ra, fa_), os.path.join(work, rb, fb_)
            label = '%s/%s vs %s/%s' % (ra, fa_[:-3], rb, fb_[:-3])
            if not (os.path.exists(fa) and os.path.exists(fb)):
                print('%-14s %-7s %-42s %10s %10s  MISSING OUTPUT' % (name, kind, label, '', ''))
                failed += 1
                continue
            if kind == 'mask':
                r = compare_mask(load(fa), load(fb), load_hkl(fb))
            else:
                r = compare(load(fa), load(fb))
            if 'shape' in r:
                print('%-14s %-7s %-42s shape mismatch %s  FAIL' % (name, kind, label, r['shape']))
                failed += 1
                continue
            tol = TOL[kind]
            ok = tol is None or r['rel'] <= tol
            res = 'info' if tol is None else ('ok' if ok else 'FAIL (tol %.0e)' % tol)
            failed += 0 if ok else 1
            print('%-14s %-7s %-42s %10.2e %10.2e  %s' % (name, kind, label, r['rel'], r['p99'], res))

    print('\nTimings (s, includes DISCUS start-up):')
    for name in selected:
        t = cases[name].get('time', {})
        if t and not cases[name].get('hidden'):
            print('  %-16s ' % name + '  '.join('%s %.2f' % (k, v) for k, v in sorted(t.items())))
    print('\n%s' % ('ALL CHECKS PASSED' if failed == 0 else '%d CHECK(S) FAILED' % failed))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
