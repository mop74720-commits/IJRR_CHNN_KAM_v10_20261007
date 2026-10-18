"""Independent solver checks, Floquet multipliers and full-FE response boundaries."""
from pathlib import Path
import argparse, json, math, time
import numpy as np
from numba import njit, prange, set_num_threads
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from run_stability_scan import (HERE, MODEL, W0, ALPHA, B0, E0, P, geometry,
    rhs, orbit, independent_batch, monodromy, save, dump)


@njit(cache=True, parallel=True)
def floquet_batch(parameters, states, omega, d, spp):
    count, nd = states.shape
    matrices = np.empty((count, nd, nd)); periodic = np.empty(count)
    for k in prange(count):
        endpoint, matrix = monodromy(states[k], parameters[k, 0]*.001,
                                     parameters[k, 1], omega, d, spp)
        matrices[k] = matrix
        periodic[k] = np.max(np.abs(endpoint-states[k]))
    return matrices, periodic


def floquet():
    records = []
    for name in ['grid', 'continuation']:
        raw = np.load(HERE/(name+'_raw.npz'))
        parameters = raw['parameters']; states = raw['strobes'][:, -1]
        matrices, residual = floquet_batch(parameters, states, *geometry(2)[:2], 512)
        eigenvalues = np.linalg.eigvals(matrices)
        radii = np.max(abs(eigenvalues), axis=1)
        periods = 2*np.pi/(W0*parameters[:, 1])
        expected_det = np.exp(-2*ALPHA*periods)
        det_error = abs(np.linalg.det(matrices)-expected_det)
        np.savez_compressed(HERE/(name+'_floquet.npz'), parameters=parameters,
            starting_states=states, monodromy=matrices, multipliers=eigenvalues,
            spectral_radius=radii, periodic_residual=residual,
            determinant_error=det_error, steps_per_period=512)
        records.append(dict(source=name, cases=len(parameters),
            spectral_radius_min=float(radii.min()), spectral_radius_max=float(radii.max()),
            max_periodic_residual=float(residual.max()),
            max_determinant_identity_error=float(det_error.max()),
            max_floquet_exponent_per_s=float(np.max(np.log(radii)/periods))))
    dump('floquet_summary.json', records); print(json.dumps(records), flush=True)


def reference_checks():
    rng = np.random.default_rng(614)
    omega, d, tip, wb, sb = geometry(2)
    worst = 0.
    for k in range(20):
        x = rng.uniform(-.3, .3, 4); t = rng.uniform(0., 3.)
        z = np.column_stack([x, np.eye(4)]); out = np.zeros_like(z)
        rhs(z, out, t, omega, d, .01, W0*.95, True)
        for j in range(4):
            eps = 1e-6
            zp = z[:, :1].copy(); zm = zp.copy(); zp[j, 0] += eps; zm[j, 0] -= eps
            op = np.zeros_like(zp); om = np.zeros_like(zm)
            rhs(zp, op, t, omega, d, .01, W0*.95, True)
            rhs(zm, om, t, omega, d, .01, W0*.95, True)
            worst = max(worst, float(abs((op[:, 0]-om[:, 0])/(2*eps)-out[:, j+1]).max()))
    checks = []
    for amp, ratio in [(.006, .82725), (.002, 1.), (.020, .8)]:
        period = 2*np.pi/(W0*ratio)
        dat = orbit(np.zeros(4), amp, ratio, *geometry(2), 8, 4, 256, True)
        # Separately written vector field, adaptive SciPy DOP853.
        def fun(t, x):
            th = d@x[:2]
            u = amp/B0*(1-np.exp(-(t/.25)**2))*np.sin(W0*ratio*t)
            return np.r_[omega*x[2:], -omega*x[:2]-ALPHA*x[2:]
                         +W0*d*(th-np.sin(th)+u*np.cos(th))]
        ts = np.arange(1, 13)*period
        sol = solve_ivp(fun, (0., ts[-1]), np.zeros(4), t_eval=ts,
                        method='DOP853', rtol=2e-12, atol=2e-14)
        assert sol.success
        discrepancy = dat[0]-sol.y.T
        checks.append(dict(amplitude_mT=amp*1e3, frequency_ratio=ratio,
                           maximum_state_difference=float(abs(discrepancy).max()),
                           maximum_tip_difference_m=float(abs(discrepancy[:, :2]@tip).max())))
    assert worst < 3e-8
    assert max(q['maximum_state_difference'] for q in checks) < 1e-9
    dump('independent_solver_checks.json', dict(analytic_jacobian_max_difference=worst,
        finite_difference_step=1e-6, adaptive_DOP853=checks,
        note='Independent RHS spelling; not a rigorous floating-point certificate.'))


def refinements():
    parameters = np.array([[a, r, 0] for a in [2., 6., 20.]
                           for r in [.8, .95, 1., 1.05, 1.2]])
    for name, n, burn, keep, spp in [('time_refined', 2, 256, 128, 512),
                                    ('four_mode', 4, 128, 64, 1024),
                                    ('eight_mode', 8, 128, 64, 4096)]:
        start = time.perf_counter()
        data = independent_batch(parameters, np.zeros((len(parameters), 2*n)),
                                 *geometry(n), burn, keep, spp)
        save(name, parameters, np.zeros((len(parameters), 2*n)), data, burn, keep, spp, n)
        print(json.dumps(dict(stage=name, seconds=time.perf_counter()-start)), flush=True)


def periodic_response(amplitude_mT, ratio, n=64, samples=1024, tolerance=2e-13):
    """FFT convolution with exact damped-oscillator transfer functions.

    This is a numerical periodic solution calculation. Attraction/stability is
    tested separately, and no contraction theorem is inferred from iteration.
    """
    omega, d, tip, wb, sb = geometry(n)
    drive = W0*ratio
    phase = 2*np.pi*np.arange(samples)/samples
    u = amplitude_mT*.001/B0*np.sin(phase)
    harmonics = np.fft.rfftfreq(samples, 1/samples)*drive
    transfer = (W0*omega[:, None]*d[:, None] /
        (omega[:, None]**2-harmonics[None, :]**2+1j*ALPHA*harmonics[None, :]))
    angle_transfer = d@transfer
    theta = np.fft.irfft(angle_transfer*np.fft.rfft(u), n=samples)
    error = np.inf
    for iteration in range(1000):
        force = theta-np.sin(theta)+u*np.cos(theta)
        updated = np.fft.irfft(angle_transfer*np.fft.rfft(force), n=samples)
        error = np.max(abs(updated-theta))
        theta = updated
        if error < tolerance:
            break
    else:
        raise RuntimeError(f'Periodic iteration did not converge: {amplitude_mT}, {ratio}, {error}')
    force = theta-np.sin(theta)+u*np.cos(theta)
    xf = transfer*np.fft.rfft(force)[None, :]
    x = np.fft.irfft(xf, n=samples, axis=1)
    y = np.fft.irfft(1j*harmonics[None, :]/omega[:, None]*xf, n=samples, axis=1)
    dx = np.fft.irfft(1j*harmonics[None, :]*xf, n=samples, axis=1)
    ddy = np.fft.irfft(-harmonics[None, :]**2/omega[:, None]*xf, n=samples, axis=1)
    reconstructed_theta = d@x
    nonlinear_rhs = (reconstructed_theta-np.sin(reconstructed_theta)
                     +u*np.cos(reconstructed_theta))
    force_rhs = -omega[:, None]*x-ALPHA*y+W0*d[:, None]*nonlinear_rhs
    residual = max(np.max(abs(ddy-force_rhs)), np.max(abs(dx-omega[:, None]*y)))
    whole = float(np.max(abs(wb@x)))
    slope = float(np.max(abs(sb@x)))
    return dict(amplitude_mT=float(amplitude_mT), frequency_ratio=float(ratio),
        modes=n, phase_samples=samples, iterations=iteration+1,
        iteration_residual_rad=float(error), ODE_residual_per_s=float(residual),
        whole_bound_mm=whole*1e3, tip_peak_mm=float(np.max(abs(tip@x)))*1e3,
        slope_bound_rad=slope, state=np.vstack([x, y]).T, phase=phase)


def roots():
    # More resolved near the narrow first-mode resonance; outside [0,20] mT
    # report no in-range crossing rather than extrapolate a threshold.
    ratios = np.unique(np.round(np.r_[np.arange(.7, 1.301, .025),
                                      np.arange(.94, 1.061, .01)], 6))
    records = []; waveforms = []; root_values = []
    for ratio in ratios:
        at20 = periodic_response(20., ratio)
        if at20['whole_bound_mm'] < 4.:
            records.append(dict(frequency_ratio=float(ratio), crossing_below_20mT=False,
                whole_at_20mT_mm=at20['whole_bound_mm']))
            continue
        function = lambda a: periodic_response(a, ratio)['whole_bound_mm']-4.
        root = brentq(function, 0., 20., xtol=1e-9)
        ans = periodic_response(root, ratio)
        waveforms.append(ans.pop('state')); ans.pop('phase')
        ans['crossing_below_20mT'] = True
        records.append(ans); root_values.append((root, ratio))
    np.savez_compressed(HERE/'full_FE_boundary_waveforms.npz',
        amplitude_and_ratio=np.asarray(root_values), states=np.asarray(waveforms),
        phase=2*np.pi*np.arange(1024)/1024)
    dump('full_FE_steady_clearance_boundary.json', records)
    # Modal and temporal convergence on representative clearance crossings.
    rows = []; traces = {}
    for ratio in [.8, .9, 1., 1.1, 1.2]:
        baseline = periodic_response(20., ratio)
        if baseline['whole_bound_mm'] < 4.: continue
        root = brentq(lambda a: periodic_response(a, ratio)['whole_bound_mm']-4.,
                      0., 20., xtol=1e-10)
        entry = dict(frequency_ratio=ratio, full_FE_root_mT=root, modal_results=[], phase_results=[])
        for n in [2, 4, 8, 16, 32, 64]:
            value = periodic_response(root, ratio, n)
            value.pop('state'); value.pop('phase'); entry['modal_results'].append(value)
        for samples in [512, 1024, 2048, 4096]:
            value = periodic_response(root, ratio, 64, samples)
            state = value.pop('state'); value.pop('phase'); entry['phase_results'].append(value)
            if samples == 1024: traces[str(ratio)] = (root, state)
        rows.append(entry)
    dump('boundary_modal_phase_refinement.json', rows)
    print(json.dumps(dict(stage='full_FE_boundaries', count=len(records), crossings=len(root_values))), flush=True)
    # Independent implicit full 64-mode integration from the Fourier orbit.
    # One initial state on this numerical periodic solution; no startup claim.
    checks = []
    omega, d, tip, wb, sb = geometry(64)
    for ratio in [.8, 1., 1.2]:
        if str(ratio) not in traces: continue
        root, states = traces[str(ratio)]
        period = 2*np.pi/(W0*ratio)
        ts = np.arange(1025)*period/1024
        reference = np.vstack([states, states[0]])
        def fun(t, z):
            theta = d@z[:64]; u = root*.001/B0*np.sin(W0*ratio*t)
            return np.r_[omega*z[64:], -omega*z[:64]-ALPHA*z[64:]
                +W0*d*(theta-np.sin(theta)+u*np.cos(theta))]
        def jac(t, z):
            theta = d@z[:64]; u = root*.001/B0*np.sin(W0*ratio*t)
            J = np.zeros((128, 128)); J[:64, 64:] = np.diag(omega)
            J[64:, :64] = -np.diag(omega)+W0*(1-np.cos(theta)-u*np.sin(theta))*np.outer(d, d)
            J[64:, 64:] = -ALPHA*np.eye(64)
            return J
        sol = solve_ivp(fun, (0., period), states[0], method='Radau', jac=jac,
                        rtol=2e-10, atol=2e-12, t_eval=ts, max_step=period/128)
        assert sol.success
        error = sol.y.T-reference
        record = dict(frequency_ratio=ratio, amplitude_mT=root,
            maximum_tip_difference_m=float(np.max(abs(error[:, :64]@tip))),
            maximum_state_difference=float(abs(error).max()), nfev=sol.nfev,
            periodic_endpoint_state_residual=float(abs(sol.y[:, -1]-sol.y[:, 0]).max()))
        np.savez_compressed(HERE/f'full_FE_Radau_r{ratio:g}.npz', time_s=ts,
            states=sol.y.T, Fourier_reference=reference, amplitude_mT=root, frequency_ratio=ratio)
        checks.append(record)
    dump('full_FE_Radau_crosscheck.json', checks)


def startup():
    records = []; saved = {}
    for ratio in [.8, .9, 1., 1.1, 1.2]:
        entry = dict(frequency_ratio=ratio, modes=[])
        for n, spp in [(2, 512), (8, 4096)]:
            memo = {}
            def evaluate(a):
                data = orbit(np.zeros(2*n), a*.001, ratio, *geometry(n),
                             128, 64, spp, True, 512)
                peak = float(data[1][:, 1].max()*1e3)
                memo[a] = (peak, data)
                return peak-4.
            root = brentq(evaluate, 0., 20., xtol=2e-7)
            peak, data = memo[root] if root in memo else (None, None)
            if data is None:
                evaluate(root); peak, data = memo[root]
            name = f'startup_{n}mode_r{ratio:g}'
            np.savez_compressed(HERE/(name+'.npz'), amplitude_mT=root,
                frequency_ratio=ratio, strobes=data[0], period_peaks=data[1],
                local_lyapunov_per_s=data[2], last_cycle_states=data[3],
                steps_per_period=spp, samples_per_period=512, burn_periods=128,
                keep_periods=64, modes=n, initial_state=np.zeros(2*n))
            entry['modes'].append(dict(modes=n, startup_crossing_mT=float(root),
                all_time_sampled_whole_mm=peak,
                steady_whole_mm=float(data[1][128:, 1].max()*1e3),
                peak_period=int(np.argmax(data[1][:, 1]))+1,
                max_slope_rad=float(data[1][:, 2].max()), evaluations=len(memo)))
            if n == 8: saved[ratio] = root
        records.append(entry)
        print(json.dumps(entry), flush=True)
    dump('startup_clearance_boundary.json', records)
    # Refine the temporal observation grid and solve the full FE transient at
    # representative 8-mode crossings. All starts are the identical zero state.
    checks = []
    for ratio in [.8, 1., 1.2]:
        root = saved[ratio]
        data = orbit(np.zeros(16), root*.001, ratio, *geometry(8),
                     256, 128, 8192, True, 1024)
        checks.append(dict(frequency_ratio=ratio, amplitude_mT=root,
            method='8-mode DOP853: doubled time and steps, doubled observations',
            maximum_whole_mm=float(data[1][:, 1].max()*1e3),
            steady_whole_mm=float(data[1][256:, 1].max()*1e3),
            maximum_slope_rad=float(data[1][:, 2].max())))
        # For off-resonant beating, the startup maximum is early. Thirty-two
        # periods cover that peak; the independently solved steady orbit is
        # also compared, without claiming an analytic tail bound.
        n = 64; omega, d, tip, wb, sb = geometry(n)
        period = 2*np.pi/(W0*ratio); ts = np.arange(32*512+1)*period/512
        def fun(t, z):
            theta = d@z[:n]
            u = root*.001/B0*(1-np.exp(-(t/.25)**2))*np.sin(W0*ratio*t)
            return np.r_[omega*z[n:], -omega*z[:n]-ALPHA*z[n:]
                +W0*d*(theta-np.sin(theta)+u*np.cos(theta))]
        def jac(t, z):
            theta = d@z[:n]
            u = root*.001/B0*(1-np.exp(-(t/.25)**2))*np.sin(W0*ratio*t)
            J = np.zeros((2*n, 2*n)); J[:n, n:] = np.diag(omega)
            J[n:, :n] = -np.diag(omega)+W0*(1-np.cos(theta)-u*np.sin(theta))*np.outer(d, d)
            J[n:, n:] = -ALPHA*np.eye(n)
            return J
        sol = solve_ivp(fun, (0., ts[-1]), np.zeros(2*n), method='Radau', jac=jac,
                        rtol=2e-9, atol=2e-11, t_eval=ts, max_step=period/128)
        assert sol.success
        whole = np.max(abs(wb@sol.y[:n]), axis=0)
        slope = np.max(abs(sb@sol.y[:n]), axis=0)
        steady = periodic_response(root, ratio, 64, 2048)
        np.savez_compressed(HERE/f'full_FE_startup_r{ratio:g}.npz', time_s=ts,
            states=sol.y.T, whole_Bernstein_bound_m=whole, slope_Bernstein_bound_rad=slope,
            amplitude_mT=root, frequency_ratio=ratio)
        checks.append(dict(frequency_ratio=ratio, amplitude_mT=root,
            method='Full-64-mode Radau from rest, first 32 periods plus periodic solution',
            maximum_whole_mm=float(whole.max()*1e3),
            steady_whole_mm=steady['whole_bound_mm'], duration_s=float(ts[-1]),
            peak_time_s=float(ts[np.argmax(whole)]), nfev=sol.nfev))
    dump('startup_boundary_refinement.json', checks)


def resonance_minimum():
    def critical_amplitude(r):
        return brentq(lambda a: periodic_response(a, r, 64, 2048)['whole_bound_mm']-4.,
                      0., 4., xtol=1e-10)
    optimum = minimize_scalar(critical_amplitude, bounds=(.98, 1.02), method='bounded',
                             options=dict(xatol=1e-9))
    value = periodic_response(optimum.fun, optimum.x, 64, 4096)
    value.pop('state'); value.pop('phase')
    value['purpose'] = 'Minimum of numerical steady 4mm amplitude contour near first resonance; not a chaos threshold or rigorous safe input.'
    value['search_interval_ratio'] = [.98, 1.02]
    value['physical_frequency_rad_s'] = float(W0*optimum.x)
    dump('minimum_steady_clearance_amplitude.json', value)


def main():
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['floquet', 'refine', 'roots', 'checks', 'startup', 'minimum'])
    p.add_argument('--threads', type=int, default=12); args = p.parse_args(); set_num_threads(args.threads)
    start = time.perf_counter()
    if args.mode == 'floquet': floquet()
    elif args.mode == 'refine': refinements()
    elif args.mode == 'roots': roots()
    elif args.mode == 'startup': startup()
    elif args.mode == 'minimum': resonance_minimum()
    else: reference_checks()
    print(json.dumps(dict(mode=args.mode, elapsed_s=time.perf_counter()-start)), flush=True)


if __name__ == '__main__': main()
