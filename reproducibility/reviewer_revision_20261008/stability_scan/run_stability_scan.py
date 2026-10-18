"""Forced-response boundaries of the frozen magnetic-beam FE model.

All outputs stay in this folder. No original evidence is modified.
Physical seconds, metres and tesla are used in the integration. The nominal
plant, not a learned surrogate, is scanned. Lyapunov exponents condition on the
prescribed forcing phase; the trivial autonomous clock exponent is excluded.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import numpy as np
from numba import njit, prange, set_num_threads
from scipy.integrate._ivp import dop853_coefficients
from threadpoolctl import threadpool_limits
threadpool_limits(1)

HERE = Path(__file__).resolve().parent
MODEL = np.load(HERE / 'inputs/reference_model_matrices.npz')
SPEC = json.loads((HERE / 'inputs/reference_model_spec.json').read_text(encoding='utf-8'))
P = SPEC['parameters_SI']
W0 = float(MODEL['omega0_rad_s'])
ALPHA = P['drag_mass_factor_s_inv']
B0 = P['bias_field_T']
E0 = float(MODEL['energy_scale_J'])
RK_A = np.ascontiguousarray(dop853_coefficients.A[:12, :12])
RK_B = np.ascontiguousarray(dop853_coefficients.B)
RK_C = np.ascontiguousarray(dop853_coefficients.C[:12])
AMPLITUDES = np.array([1., 2., 4., 6., 10., 15., 20.])
RATIOS = np.unique(np.round(np.r_[np.arange(.7, 1.301, .02),
                                 np.arange(.94, 1.061, .005)], 6))
INITIAL_LOW = np.array([[0., 0., 0., 0.], [.18, .04, 0., 0.],
                        [-.18, -.04, 0., 0.]])


def dump(name, value):
    text = json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n'
    (HERE / name).write_bytes(text.replace('\r\n', '\n').replace('\n', '\r\n').encode('utf-8'))


def geometry(n):
    omega = np.asarray(MODEL['all_omega'][:n])
    scale = np.sqrt(E0 / (W0 * omega))
    shape = np.asarray(MODEL['all_phi'][:, :n]) * scale
    nodal = np.vstack([np.zeros((2, n)), shape])
    element_length = P['length_m'] / 32
    wb, sb = [], []
    for e in range(32):
        w0, t0, w1, t1 = nodal[2*e:2*e+4]
        wb.extend([w0, w0+element_length*t0/3,
                   w1-element_length*t1/3, w1])
        sb.extend([t0, 3*(w1-w0)/element_length-t0-t1, t1])
    return omega, shape[-1].copy(), shape[-2].copy(), np.array(wb), np.array(sb)


@njit(cache=True)
def rhs(z, out, t, omega, d, amp, drive, ramp):
    n = len(omega)
    angle = 0.
    for j in range(n):
        angle += d[j]*z[j, 0]
    u = amp/B0 * math.sin(drive*t)
    if ramp:
        u *= 1.-math.exp(-(t/.25)**2)
    nonlinear = angle-math.sin(angle)+u*math.cos(angle)
    derivative = 1.-math.cos(angle)-u*math.sin(angle)
    for j in range(n):
        out[j, 0] = omega[j]*z[n+j, 0]
        out[n+j, 0] = -omega[j]*z[j, 0]-ALPHA*z[n+j, 0]+W0*d[j]*nonlinear
    for k in range(1, z.shape[1]):
        delta_angle = 0.
        for j in range(n):
            delta_angle += d[j]*z[j, k]
        for j in range(n):
            out[j, k] = omega[j]*z[n+j, k]
            out[n+j, k] = (-omega[j]*z[j, k]-ALPHA*z[n+j, k]
                           +W0*d[j]*derivative*delta_angle)


@njit(cache=True)
def step(z, work, trial, out, t, h, omega, d, amp, drive, ramp):
    nr, nc = z.shape
    for stage in range(12):
        for row in range(nr):
            for col in range(nc):
                v = z[row, col]
                for j in range(stage):
                    v += h*RK_A[stage, j]*work[j, row, col]
                trial[row, col] = v
        rhs(trial, out, t+RK_C[stage]*h, omega, d, amp, drive, ramp)
        for row in range(nr):
            for col in range(nc):
                work[stage, row, col] = out[row, col]
    for row in range(nr):
        for col in range(nc):
            v = z[row, col]
            for stage in range(12):
                v += h*RK_B[stage]*work[stage, row, col]
            z[row, col] = v


@njit(cache=True)
def spatial_bounds(x, wb, sb):
    displacement, slope = 0., 0.
    for i in range(len(wb)):
        v = 0.
        for j in range(len(x)):
            v += wb[i, j]*x[j]
        displacement = max(displacement, abs(v))
    for i in range(len(sb)):
        v = 0.
        for j in range(len(x)):
            v += sb[i, j]*x[j]
        slope = max(slope, abs(v))
    return displacement, slope


@njit(cache=True)
def orbit(x0, amp, ratio, omega, d, tip, wb, sb, burn, keep, spp, ramp, phase_samples=128):
    n = len(omega)
    period = 2*math.pi/(W0*ratio)
    h = period/spp
    z = np.zeros((2*n, 2))
    z[:, 0] = x0
    for j in range(2*n):
        z[j, 1] = 1./math.sqrt(2*n)
    work = np.zeros((12, 2*n, 2)); trial = z.copy(); out = z.copy()
    periods = burn+keep
    strobes = np.empty((periods, 2*n))
    peaks = np.zeros((periods, 3))
    exponents = np.zeros(periods)
    trace = np.zeros((spp+1, 2*n))
    # Spatial checks use 128 phases/cycle; the final cycle retains every step.
    skip = max(1, spp//phase_samples)
    for p in range(periods):
        for s in range(spp):
            if p == periods-1:
                trace[s] = z[:, 0]
            if s % skip == 0:
                w, slope = spatial_bounds(z[:n, 0], wb, sb)
                ti = 0.
                for j in range(n): ti += tip[j]*z[j, 0]
                peaks[p, 0] = max(peaks[p, 0], abs(ti))
                peaks[p, 1] = max(peaks[p, 1], w)
                peaks[p, 2] = max(peaks[p, 2], slope)
            t = (p*spp+s)*h
            step(z, work, trial, out, t, h, omega, d, amp, W0*ratio, ramp)
        norm = math.sqrt(np.sum(z[:, 1]*z[:, 1]))
        exponents[p] = math.log(norm)/period
        z[:, 1] /= norm
        strobes[p] = z[:, 0]
    trace[-1] = z[:, 0]
    return strobes, peaks, exponents, trace


@njit(cache=True, parallel=True)
def independent_batch(params, initial, omega, d, tip, wb, sb, burn, keep, spp):
    count, nd = initial.shape
    strobes = np.empty((count, burn+keep, nd))
    peaks = np.empty((count, burn+keep, 3))
    exponents = np.empty((count, burn+keep))
    traces = np.empty((count, spp+1, nd))
    for i in prange(count):
        a, b, c, e = orbit(initial[i], params[i, 0]*.001, params[i, 1],
                            omega, d, tip, wb, sb, burn, keep, spp, True)
        strobes[i] = a; peaks[i] = b; exponents[i] = c; traces[i] = e
    return strobes, peaks, exponents, traces


@njit(cache=True, parallel=True)
def continuation_batch(amplitudes, ratios, omega, d, tip, wb, sb, burn, keep, spp):
    count = len(amplitudes)*2
    nd = 2*len(omega); nr = len(ratios)
    strobes = np.empty((count, nr, burn+keep, nd))
    peaks = np.empty((count, nr, burn+keep, 3))
    exponents = np.empty((count, nr, burn+keep))
    traces = np.empty((count, nr, spp+1, nd))
    for row in prange(count):
        x = np.zeros(nd)
        for k in range(nr):
            ri = k if row % 2 == 0 else nr-1-k
            a, b, c, e = orbit(x, amplitudes[row//2]*.001, ratios[ri],
                              omega, d, tip, wb, sb, burn, keep, spp, k == 0)
            strobes[row, ri] = a; peaks[row, ri] = b
            exponents[row, ri] = c; traces[row, ri] = e
            x = a[-1].copy()
    return strobes, peaks, exponents, traces


@njit(cache=True)
def monodromy(x, amp, ratio, omega, d, spp):
    nd = len(x)
    z = np.zeros((nd, nd+1)); z[:, 0] = x
    for j in range(nd): z[j, j+1] = 1.
    work = np.zeros((12, nd, nd+1)); trial = z.copy(); out = z.copy()
    period = 2*math.pi/(W0*ratio); h = period/spp
    for k in range(spp):
        step(z, work, trial, out, k*h, h, omega, d, amp, W0*ratio, False)
    return z[:, 0], z[:, 1:]


def summarize(params, data, burn, keep):
    strobes, peaks, exponents, traces = data
    rows = []
    for i, (amp, ratio, initial) in enumerate(params):
        st = strobes[i, burn:]
        period = 2*np.pi/(W0*ratio)
        rows.append(dict(amplitude_mT=float(amp), frequency_ratio=float(ratio),
            initial_condition=int(initial), duration_s=float((burn+keep)*period),
            steady_tip_peak_mm=float(peaks[i, burn:, 0].max()*1e3),
            steady_whole_bound_mm=float(peaks[i, burn:, 1].max()*1e3),
            all_time_sampled_whole_bound_mm=float(peaks[i, :, 1].max()*1e3),
            steady_slope_bound_rad=float(peaks[i, burn:, 2].max()),
            conditional_FTLE_per_s=float(exponents[i, burn:].mean()),
            first_half_FTLE_per_s=float(exponents[i, burn:burn+keep//2].mean()),
            second_half_FTLE_per_s=float(exponents[i, burn+keep//2:].mean()),
            last_strobe_step_max=float(np.max(np.abs(np.diff(st[-16:], axis=0)))),
            strobe_diameter_last_16=float(np.max(np.ptp(st[-16:], axis=0))),
            last_two_cycle_peak_change_mm=float(abs(peaks[i, -1, 1]-peaks[i, -2, 1])*1e3)))
    return rows


def save(name, params, initial, data, burn, keep, spp, modes):
    np.savez_compressed(HERE/(name+'_raw.npz'), parameters=params, initial_states=initial,
        strobes=data[0], period_peaks=data[1], local_lyapunov_per_s=data[2],
        last_cycle_states=data[3], burn_periods=burn, keep_periods=keep,
        steps_per_period=spp, modes=modes)
    rows = summarize(params, data, burn, keep)
    dump(name+'_results.json', rows)
    with (HERE/(name+'_results.csv')).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['pilot', 'grid', 'continuation'])
    parser.add_argument('--threads', type=int, default=12)
    args = parser.parse_args(); set_num_threads(args.threads)
    start = time.perf_counter(); geom = geometry(2)
    if args.mode == 'pilot':
        params = np.array([[a, r, 0] for a in [2., 20.] for r in [.8, 1., 1.2]])
        initial = np.zeros((len(params), 4)); burn, keep, spp = 128, 64, 256
        data = independent_batch(params, initial, *geom, burn, keep, spp)
        rows = save('pilot', params, initial, data, burn, keep, spp, 2)
    elif args.mode == 'grid':
        params = np.array([[a, r, i] for a in AMPLITUDES for r in RATIOS for i in range(3)])
        initial = np.array([INITIAL_LOW[int(p[2])] for p in params]); burn, keep, spp = 128, 64, 256
        data = independent_batch(params, initial, *geom, burn, keep, spp)
        rows = save('grid', params, initial, data, burn, keep, spp, 2)
    else:
        burn, keep, spp = 128, 64, 256
        data = continuation_batch(AMPLITUDES, RATIOS, *geom, burn, keep, spp)
        params = np.array([[a, r, -1 if direction == 0 else -2]
                          for a in AMPLITUDES for direction in range(2) for r in RATIOS])
        flat = tuple(v.reshape((-1,)+v.shape[2:]) for v in data)
        rows = save('continuation', params, np.empty((0, 4)), flat, burn, keep, spp, 2)
    elapsed = time.perf_counter()-start
    info = dict(mode=args.mode, cases=len(rows), elapsed_s=elapsed,
        max_FTLE=max(v['conditional_FTLE_per_s'] for v in rows),
        max_strobe_step=max(v['last_strobe_step_max'] for v in rows),
        max_whole_mm=max(v['steady_whole_bound_mm'] for v in rows),
        max_slope=max(v['steady_slope_bound_rad'] for v in rows))
    dump(args.mode+'_run.json', info); print(json.dumps(info), flush=True)


if __name__ == '__main__':
    main()
