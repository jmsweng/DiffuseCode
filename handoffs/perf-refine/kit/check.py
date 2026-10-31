#!/usr/bin/env python3
"""Validation kit for the perf/refine branch of DiffuseCode.

Usage (inside WSL, from the kit directory):
    python3 check.py <item>      item = 1 .. 7, or "all" for the summary table

Binaries are expected as ../bins/c0 (upstream base) .. ../bins/c7 (last commit),
built by setup.sh. Item N compares bins/c(N-1) with bins/cN.
Every check prints PASS or FAIL; the last line is the overall result.
Raw outputs are kept in ../runs/<binary>/<case>/ for inspection.
"""
import os, re, shutil, subprocess, sys

KIT = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(KIT)
BINS = os.path.join(ROOT, 'bins')
RUNS = os.path.join(ROOT, 'runs')
SEED = 'seed 4711, 271828, 314159, 161803, 141421, 173205, 223606, 99991'

# name: (work macro, parameter file, extra option for 'run')
CASES = {
    'det':     ('work_det.mac', 'p_det.txt',     ''),
    'far':     ('work_det.mac', 'p_far.txt',     ''),
    'farfix':  ('work_det.mac', 'p_farfix.txt',  ''),
    'range':   ('work_det.mac', 'p_range.txt',   ''),
    'p5':      ('work_det.mac', 'p_p5.txt',      ''),
    'one':     ('work_det.mac', 'p_one.txt',     ''),
    'plot':    ('work_det.mac', 'p_det.txt',     ', plot:kplot.mac'),
    'ran':     ('work_ran.mac', 'p_det.txt',     ''),
    'ranfix':  ('work_ran.mac', 'p_ranfix.txt',  ''),
    'p3range': ('work_det.mac', 'p_p3range.txt', ''),
    'p2':      ('work_det.mac', 'p_p2.txt',      ''),
    'p2far':   ('work_det.mac', 'p_p2far.txt',   ''),
    'p2range': ('work_det.mac', 'p_p2range.txt', ''),
    'p2ran':   ('work_ran.mac', 'p_p2.txt',      ''),
}
OLD_CASES = ['det', 'far', 'farfix', 'range', 'p5', 'one', 'plot', 'ran', 'ranfix', 'p3range']
P2_CASES = ['p2', 'p2far', 'p2range', 'p2ran']
P2_REF = {'p2': 'det', 'p2far': 'far', 'p2range': 'p3range', 'p2ran': 'ran'}

failures = []


def check(ok, text):
    print(('PASS  ' if ok else 'FAIL  ') + text)
    if not ok:
        failures.append(text)


def binary(i):
    path = os.path.join(BINS, 'c%d' % i)
    if not os.path.isfile(path):
        sys.exit('missing binary %s - run setup.sh first' % path)
    return path


def run_refine(i, name):
    work, pfile, runopt = CASES[name]
    d = os.path.join(RUNS, 'c%d' % i, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    for f in (work, pfile, 'gauss.data', 'kplot.mac'):
        shutil.copy(os.path.join(KIT, 'refine', f), d)
    lines = ['set prompt,off', SEED, 'i[90] = 0', 'i[91] = 0', 'refine', 'data xy, gauss.data']
    lines += open(os.path.join(d, pfile)).read().splitlines()
    lines += ['set cycle, 25', 'run %s%s' % (work, runopt), 'exit',
              'eval i[90]', 'eval i[91]', 'exit']
    open(os.path.join(d, 'main.mac'), 'w').write('\n'.join(lines) + '\n')
    with open(os.path.join(d, 'out.txt'), 'w') as out:
        try:
            subprocess.run([binary(i), '-macro', 'main.mac'], cwd=d, stdout=out,
                           stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, timeout=600)
        except subprocess.TimeoutExpired:
            out.write('\nTIMEOUT\n')
    return parse_refine(open(os.path.join(d, 'out.txt'), errors='replace').read())


def parse_refine(t):
    r = {'table': [], 'pars': {}, 'esd': {}, 'chi2': None, 'runs': None, 'plots': None,
         'errors': [], 'shift': []}
    for line in t.splitlines():
        m = re.match(r'^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
        if m:
            r['table'].append(line.strip())
            r['shift'].append(float(m.group(3)))
        m = re.match(r'^\s+(P_\w+)\s+:\s+(\S+)\s+\+-\s+(\S+)', line)
        if m:
            r['pars'][m.group(1)] = float(m.group(2))
            r['esd'][m.group(1)] = float(m.group(3))
        m = re.search(r'Chi\^2\s+:\s+(\S+)', line)
        if m:
            r['chi2'] = float(m.group(1))
        m = re.search(r'Value of i\[90\] =\s*(\S+)', line)
        if m:
            r['runs'] = int(float(m.group(1)))
        m = re.search(r'Value of i\[91\] =\s*(\S+)', line)
        if m:
            r['plots'] = int(float(m.group(1)))
        if '***' in line and 'Version' not in line and '*****' not in line:
            r['errors'].append(re.sub(r'\x1b\[[0-9;]*m', '', line).strip())
        if 'TIMEOUT' in line:
            r['errors'].append('TIMEOUT')
    return r


def results_key(r, with_shift=True):
    # without shift: drop the 'MAX(dP/sig)' and 'Par' columns (3 and 4)
    table = r['table'] if with_shift else [' '.join(x.split()[:2] + x.split()[4:]) for x in r['table']]
    return (table, sorted(r['pars'].items()), sorted(r['esd'].items()), r['chi2'], r['errors'])


def ok_run(r):
    return not r['errors'] and r['runs'] is not None and r['pars']


def gauss_ok(r, tol_amp=0.05, tol_cen=0.01, tol_wid=0.02, tol_bg=0.03):
    p = r['pars']
    return (abs(p.get('P_amp', 0) - 8.0) < tol_amp and abs(p.get('P_cen', 0) - 0.7) < tol_cen
            and abs(p.get('P_wid', 0) - 1.3) < tol_wid and abs(p.get('P_bg', 9) - 0.5) < tol_bg)


def run_rmc(i, tag, lots, fourier_lots=False):
    d = os.path.join(RUNS, 'c%d' % i, 'rmc_' + tag)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    for f in ('cu.cll', 'exp.nipl'):
        shutil.copy(os.path.join(KIT, 'rmc', f), d)
    body = open(os.path.join(KIT, 'rmc', 'rmc_body.inc')).read()
    if lots:
        la, lb = 'set lots, box, 3,3,3, 4, yes', 'set lots, box, 3,3,3, 4, yes, lots.txt'
    else:
        la = lb = 'set lots, off'
    flots = 'fourier\nlots box, 3,3,3, 4, yes\nexit\n' if fourier_lots else ''
    open(os.path.join(d, 'a.inc'), 'w').write(
        body.replace('SIGMA', '-9999.').replace('NCYC', '300').replace('LOTSLINE', la))
    open(os.path.join(d, 'b.inc'), 'w').write(
        body.replace('SIGMA', '-9999.').replace('NCYC', '1').replace('LOTSLINE', lb))
    open(os.path.join(d, 'A.mac'), 'w').write(
        'set prompt,off\ndiscus\nseed 1234\nread\ncell cu.cll, 6,6,6\n' + flots +
        '@a.inc\nset move, all, 0.05,0.05,0.05\nrun\nsave lots, lots.txt\nexit\neval res[1]\n'
        'save\noutf final.stru\nrun\nexit\nexit\nexit\n')
    open(os.path.join(d, 'B.mac'), 'w').write(
        'set prompt,off\ndiscus\nseed 1234\nread\nstru final.stru\n' + flots +
        '@b.inc\nset move, all, 0.0,0.0,0.0\nrun\nexit\neval res[1]\nexit\nexit\n')
    vals = {}
    for m in ('A', 'B'):
        with open(os.path.join(d, m + '.out'), 'w') as out:
            try:
                subprocess.run([binary(i), '-macro', m + '.mac'], cwd=d, stdout=out,
                               stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, timeout=300)
            except subprocess.TimeoutExpired:
                out.write('\nTIMEOUT\n')
        t = open(os.path.join(d, m + '.out'), errors='replace').read()
        v = re.findall(r'Value of res\[1\] =\s*(\S+)', t)
        vals[m] = float(v[-1]) if v else None
        if 'No valid move' in t:
            vals[m + '_err'] = 'no valid move'
    return vals


def rel(a, b):
    if a is None or b is None:
        return None
    return abs(a - b) / max(abs(a), abs(b))


def fmt_runs(rs):
    return '/'.join(str(r['runs']) for r in rs)


# --------------------------------------------------------------------------- items

def item1():
    print('Item 01: MPI derivative step sizes (c0 -> c1)')
    src = subprocess.run(['git', '-C', os.path.join(ROOT, 'src'), 'show', sha(1) + ':refine/prog/refine_run.f90'],
                         capture_output=True, text=True).stdout
    body = src.split('subroutine refine_calc_deriv_mpi')[1].split('end subroutine refine_calc_deriv_mpi')[0]
    bad = [l for l in body.splitlines() if re.match(r'^\s*delta\s*=', l)]
    check(not bad, 'refine_calc_deriv_mpi has no whole-array "delta =" assignment (%d found)' % len(bad))
    for name in ['det', 'p5', 'ran']:
        a, b = run_refine(0, name), run_refine(1, name)
        check(ok_run(b) and results_key(a) == results_key(b),
              '%-7s serial refinement identical to base (MPI-only change), runs %s' % (name, fmt_runs([a, b])))


def item2():
    print('Item 02: seed restore reproduces Gaussian random numbers (c1 -> c2)')
    for name in ['det', 'far', 'p5', 'range']:
        a, b = run_refine(1, name), run_refine(2, name)
        check(ok_run(b) and results_key(a) == results_key(b), '%-7s deterministic fit unchanged' % name)
    a, b = run_refine(1, 'ran'), run_refine(2, 'ran')
    print('INFO  ran before: cycles %d chi2 %s pars %s' % (len(a['table']), a['chi2'], a['pars']))
    print('INFO  ran after : cycles %d chi2 %s pars %s' % (len(b['table']), b['chi2'], b['pars']))
    check(ok_run(b) and gauss_ok(b), 'ran     noisy fit reaches the true parameters (8, 0.7, 1.3, 0.5)')
    check(b['chi2'] is not None and b['chi2'] < 150, 'ran     chi^2 < 150 (101 points, noise sigma 0.02 on 0.05)')


def item3():
    print('Item 03: RMC with lots (c2 -> c3); incremental chi^2 vs fresh recalculation')
    for tag, lots, flots in [('nolots', False, False), ('lots', True, False), ('fourlots', True, True)]:
        a, b = run_rmc(2, tag, lots, flots), run_rmc(3, tag, lots, flots)
        print('INFO  %-8s before: incremental %s fresh %s %s' % (tag, a.get('A'), a.get('B'), a.get('A_err', '')))
        print('INFO  %-8s after : incremental %s fresh %s %s' % (tag, b.get('A'), b.get('B'), b.get('A_err', '')))
        r = rel(b.get('A'), b.get('B'))
        check(r is not None and r < 1e-6, '%-8s after: incremental and fresh chi^2 agree (rel. diff %s)' % (tag, r))
        if tag == 'nolots':
            ra = rel(a.get('A'), a.get('B'))
            check(ra is not None and ra < 1e-6, 'nolots   before: also consistent (control)')


def item4():
    print('Item 04: dP/sigma from parameter sigmas (c3 -> c4)')
    for name in OLD_CASES:
        a, b = run_refine(3, name), run_refine(4, name)
        check(ok_run(b) and results_key(a, False) == results_key(b, False) and a['runs'] == b['runs'],
              '%-7s cycles, chi^2, parameters, esd unchanged apart from the dP/sig column' % name)
    a, b = run_refine(3, 'det'), run_refine(4, 'det')
    print('INFO  det dP/sig before: %s' % a['shift'])
    print('INFO  det dP/sig after : %s' % b['shift'])
    check(b['shift'] and b['shift'][0] > 10, 'det     first-cycle dP/sig is now a multiple of sigma (> 10)')


def item5():
    print('Item 05: skip repeated centre run (c4 -> c5)')
    for name in OLD_CASES:
        a, b = run_refine(4, name), run_refine(5, name)
        same = ok_run(b) and results_key(a) == results_key(b)
        fewer = b['runs'] is not None and a['runs'] is not None and b['runs'] < a['runs']
        check(same and fewer, '%-7s identical results, macro runs %d -> %d' % (name, a['runs'] or -1, b['runs'] or -1))
    a, b = run_refine(4, 'plot'), run_refine(5, 'plot')
    check(a['runs'] - b['runs'] == 1 and a['plots'] == b['plots'],
          'plot    with plot: only the final repeat is saved (%d -> %d), plots %s -> %s'
          % (a['runs'], b['runs'], a['plots'], b['plots']))


def item6():
    print('Item 06: derivatives only for accepted steps (c5 -> c6)')
    for name in OLD_CASES:
        a, b = run_refine(5, name), run_refine(6, name)
        same = ok_run(b) and results_key(a) == results_key(b)
        rejected = sum(1 for k in range(1, len(a['shift'])) if a['shift'][k] == 0.0)
        if rejected:
            ok = same and b['runs'] < a['runs']
        else:
            ok = same and b['runs'] <= a['runs']
        check(ok, '%-7s identical results, macro runs %d -> %d' % (name, a['runs'] or -1, b['runs'] or -1))


def item7():
    print('Item 07: points:2 one-sided derivatives (c6 -> c7)')
    for name in OLD_CASES:
        a, b = run_refine(6, name), run_refine(7, name)
        check(ok_run(b) and results_key(a) == results_key(b) and a['runs'] == b['runs'],
              '%-7s points:3/5 fit unchanged' % name)
    for name in P2_CASES:
        a = run_refine(6, name)
        check(bool(a['errors']), '%-7s before: points:2 is rejected with an error' % name)
        b, ref = run_refine(7, name), run_refine(7, P2_REF[name])
        ok = ok_run(b)
        dmax = max(abs(b['pars'][k] - ref['pars'][k]) / max(ref['esd'][k], 1e-12) for k in ref['pars']) if ok else None
        check(ok and dmax < 0.2 and b['runs'] < ref['runs'],
              '%-7s converges to the points:3 result within %.3f esd, macro runs %d (points:3: %d)'
              % (name, dmax if dmax is not None else -1, b['runs'] or -1, ref['runs'] or -1))


def summary():
    print('Macro runs per fit, base (c0) and every commit:')
    cols = list(range(0, 8))
    print('%-8s' % 'case' + ''.join('%7s' % ('c%d' % i) for i in cols))
    for name in OLD_CASES + P2_CASES:
        row = []
        for i in cols:
            r = run_refine(i, name)
            row.append('err' if not ok_run(r) else str(r['runs']))
        print('%-8s' % name + ''.join('%7s' % x for x in row))


def sha(i):
    return open(os.path.join(BINS, 'c%d.sha' % i)).read().strip()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    arg = sys.argv[1]
    if arg == 'all':
        summary()
        sys.exit(0)
    {'1': item1, '2': item2, '3': item3, '4': item4, '5': item5, '6': item6, '7': item7}[arg.lstrip('0')]()
    print('RESULT: %s' % ('PASS' if not failures else 'FAIL (%d)' % len(failures)))
