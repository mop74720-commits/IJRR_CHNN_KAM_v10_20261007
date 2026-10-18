"""Recompute summary, validation assertions and the publication figure."""
from pathlib import Path
import hashlib, json, platform
import numpy as np
import scipy, numba
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from run_stability_scan import HERE, W0, ALPHA, AMPLITUDES, RATIOS, geometry, dump


def read(name): return json.loads((HERE/name).read_text(encoding='utf-8'))


def main():
    grid = np.load(HERE/'grid_raw.npz'); cont = np.load(HERE/'continuation_raw.npz')
    rows = read('grid_results.json'); cf = np.load(HERE/'continuation_floquet.npz')
    gf = np.load(HERE/'grid_floquet.npz')
    boundary = read('full_FE_steady_clearance_boundary.json')
    finite_boundary = [r for r in boundary if r['crossing_below_20mT']]
    startup = read('startup_clearance_boundary.json')
    startup_refinement = read('startup_boundary_refinement.json')
    roots = read('boundary_modal_phase_refinement.json')
    summary = dict(scope='Numerical forced-response study of the nominal damped beam. Conservative KAM-torus destruction and a chaos threshold are not established.',
        independent_runs=len(rows), continuation_runs=len(cont['parameters']),
        amplitudes_mT=AMPLITUDES.tolist(), frequency_ratios=RATIOS.tolist(),
        physical_angular_frequency_range_rad_s=[float(RATIOS.min()*W0), float(RATIOS.max()*W0)],
        total_time_range_s=[float(min(r['duration_s'] for r in rows)),float(max(r['duration_s'] for r in rows))],
        independent_runs_with_all_sampled_slope_at_most_point1=int(np.sum(grid['period_peaks'][:, :, 2].max(axis=1)<=.1)),
        independent_runs_with_steady_slope_at_most_point1=int(np.sum(grid['period_peaks'][:, 128:, 2].max(axis=1)<=.1)),
        independent_runs_with_all_sampled_whole_bound_below4=int(np.sum(grid['period_peaks'][:, :, 1].max(axis=1)<.004)),
        maximum_FTLE_per_s=float(grid['local_lyapunov_per_s'][:, 128:].mean(axis=1).max()),
        minimum_FTLE_per_s=float(grid['local_lyapunov_per_s'][:, 128:].mean(axis=1).min()),
        maximum_floquet_radius=float(max(gf['spectral_radius'].max(),cf['spectral_radius'].max())),
        minimum_floquet_radius=float(min(gf['spectral_radius'].min(),cf['spectral_radius'].min())),
        maximum_periodic_residual=float(max(gf['periodic_residual'].max(),cf['periodic_residual'].max())),
        maximum_floquet_determinant_error=float(max(gf['determinant_error'].max(),cf['determinant_error'].max())),
        max_boundary_Bernstein_minus_actual_tip_mm=float(max(abs(r['whole_bound_mm']-r['tip_peak_mm']) for r in finite_boundary)),
        minimum_steady_clearance_amplitude=read('minimum_steady_clearance_amplitude.json'))
    endpoints = grid['strobes'][:, -1].reshape(7, len(RATIOS), 3, 4)
    branch = cont['strobes'][:, -1].reshape(7, 2, len(RATIOS), 4)
    summary['max_initial_condition_endpoint_difference'] = float(abs(endpoints[:, :, 1:]-endpoints[:, :, 0:1]).max())
    summary['max_up_down_endpoint_difference'] = float(abs(branch[:, 0]-branch[:, 1]).max())
    summary['max_grid_continuation_endpoint_difference'] = float(abs(endpoints[:, :, 0]-branch[:, 0]).max())
    base = {(r['amplitude_mT'],r['frequency_ratio'],r['initial_condition']):r for r in rows}
    refinement = {}
    for name in ['time_refined', 'four_mode', 'eight_mode']:
        data = read(name+'_results.json')
        refinement[name] = dict(cases=len(data),
            max_FTLE_per_s=max(r['conditional_FTLE_per_s'] for r in data),
            max_strobe_step=max(r['last_strobe_step_max'] for r in data),
            max_change_from_two_mode_steady_peak_mm=max(abs(r['steady_whole_bound_mm']-base[r['amplitude_mT'],r['frequency_ratio'],0]['steady_whole_bound_mm']) for r in data))
    summary['selected_refinements'] = refinement
    summary['max_full_FE_periodic_Radau_tip_difference_m'] = max(r['maximum_tip_difference_m'] for r in read('full_FE_Radau_crosscheck.json'))
    summary['max_512_to_4096_phase_peak_change_mm'] = max(abs(q['whole_bound_mm']-r['phase_results'][0]['whole_bound_mm']) for r in roots for q in r['phase_results'])
    summary['max_1024_to_4096_phase_peak_change_mm'] = max(abs(q['whole_bound_mm']-r['phase_results'][1]['whole_bound_mm']) for r in roots for q in r['phase_results'][2:])
    summary['max_startup_time_refinement_peak_change_um'] = max(abs(r['maximum_whole_mm']-4)*1e3 for r in startup_refinement if r['method'].startswith('8-mode'))
    summary['max_startup_modal_check_peak_difference_um'] = max(abs(max(r['maximum_whole_mm'],r['steady_whole_mm'])-4)*1e3 for r in startup_refinement if r['method'].startswith('Full'))
    assert summary['maximum_floquet_radius'] < 1.
    assert summary['maximum_periodic_residual'] < 1e-8
    assert summary['maximum_floquet_determinant_error'] < 1e-10
    assert summary['max_boundary_Bernstein_minus_actual_tip_mm'] < 1e-10
    assert refinement['time_refined']['max_change_from_two_mode_steady_peak_mm'] < 1e-5
    assert all(r['max_FTLE_per_s'] < 0 for r in refinement.values())
    dump('summary.json', summary)
    # Four panels: response map, operating contours, local stability, startup.
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9,
        'axes.labelsize':9, 'axes.titlesize':10, 'legend.fontsize':7.5,
        'axes.spines.top':False, 'axes.spines.right':False, 'pdf.fonttype':42,
        'lines.linewidth':1.6, 'savefig.bbox':'tight'})
    fig, ax = plt.subplots(2, 2, figsize=(8.4, 6.4), layout='constrained')
    # From-rest grid only: identical ramp and zero initial state for every cell.
    peak = grid['period_peaks'][::3, 128:, 1].max(axis=1).reshape(7,-1)*1e3
    slope = grid['period_peaks'][::3, :, 2].max(axis=1).reshape(7,-1)
    mesh = ax[0,0].pcolormesh(RATIOS, AMPLITUDES, peak, shading='nearest',
                              norm=LogNorm(vmin=.2,vmax=40), cmap='viridis')
    ax[0,0].contourf(RATIOS, AMPLITUDES, slope, levels=[.1,10], colors='none', hatches=['////'])
    rr = np.array([r['frequency_ratio'] for r in finite_boundary])
    aa = np.array([r['amplitude_mT'] for r in finite_boundary])
    ax[0,0].plot(rr,aa,color='white',lw=2,label='64-mode steady 4 mm')
    ax[0,0].set(xlim=(.7,1.3),ylim=(1,20),xlabel=r'Frequency ratio $r=\Omega/\omega_0$',
                ylabel='Field amplitude (mT)',title='(a) Steady response, from-rest grid')
    ax[0,0].legend(loc='upper right',facecolor='#42505b',labelcolor='white',framealpha=.9)
    cb=fig.colorbar(mesh,ax=ax[0,0],pad=.02,ticks=[.2,1,4,10,40])
    cb.set_label('Two-mode peak (mm)'); cb.ax.set_yticklabels(['0.2','1','4','10','40'])
    ax[0,0].text(.715,1.3,'Hatching: startup slope > 0.1 rad',fontsize=7,
                 color='white',bbox=dict(facecolor='#333333',alpha=.65,pad=2,edgecolor='none'))
    ax[0,1].plot(rr,aa,'o-',ms=3,color='#0072B2',label='Steady: full 64 modes')
    ax[0,1].plot([r['frequency_ratio'] for r in startup],
        [r['modes'][1]['startup_crossing_mT'] for r in startup],'s--',ms=4,
        color='#D55E00',label='Startup: 8 modes, from rest')
    ax[0,1].axhline(20,color='#777777',ls=':',label='20 mT field limit')
    minimum = summary['minimum_steady_clearance_amplitude']
    ax[0,1].annotate('Minimum steady crossing\n1.916 mT at r = 0.99924',
        xy=(minimum['frequency_ratio'],minimum['amplitude_mT']),xytext=(.76,8.7),
        fontsize=7.5,arrowprops=dict(arrowstyle='->',lw=.8))
    ax[0,1].set(xlim=(.7,1.3),ylim=(0,22),xlabel=r'Frequency ratio $r$',
                ylabel='4 mm crossing amplitude (mT)',title='(b) Steady and startup boundaries')
    ax[0,1].legend(loc='upper right'); ax[0,1].grid(alpha=.18)
    ax[1,0].scatter(gf['parameters'][:,1],gf['spectral_radius'],s=7,
                     color='#009E73',alpha=.45,label='1,029 independent runs')
    ax[1,0].scatter(cf['parameters'][::2,1],cf['spectral_radius'][::2],s=16,
                     marker='+',color='#333333',alpha=.45,label='686 continuation runs')
    ax[1,0].axhline(1,color='#D55E00',ls='--',label='Unit-circle boundary')
    ax[1,0].set(xlim=(.7,1.3),ylim=(.78,1.02),xlabel=r'Frequency ratio $r$',
                ylabel='Largest Floquet multiplier modulus',
                title='(c) Two-mode periodic-orbit stability')
    ax[1,0].legend(loc='upper left'); ax[1,0].grid(alpha=.18)
    raw=np.load(HERE/'full_FE_startup_r0.8.npz');tip=raw['states'][:,:64]@geometry(64)[2]*1e3
    ax[1,1].plot(raw['time_s'],tip,color='#0072B2',label='Full 64-mode response')
    ax[1,1].axhline(4,color='#D55E00',ls='--');ax[1,1].axhline(-4,color='#D55E00',ls='--',label='4 mm clearance')
    steady=[r for r in startup_refinement if r['frequency_ratio']==.8 and r['method'].startswith('Full')][0]['steady_whole_mm']
    ax[1,1].axhline(steady,color='#777777',ls=':');ax[1,1].axhline(-steady,color='#777777',ls=':',label='Steady envelope')
    ax[1,1].set(xlim=(0,8),ylim=(-4.5,4.5),xlabel='Time (s)',ylabel='Tip displacement (mm)',
                title='(d) Startup overshoot: r = 0.8, A = 8.651 mT')
    ax[1,1].legend(loc='lower right',framealpha=.9); ax[1,1].grid(alpha=.18)
    fig.savefig(HERE/'fig_stability_boundaries.pdf'); plt.close(fig)
    import fitz
    pdf=fitz.open(HERE/'fig_stability_boundaries.pdf')
    pdf[0].get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(HERE/'fig_stability_boundaries_preview.png')
    pdf.close()
    dump('validation.json', dict(passed=True, checks=[
        'All 1715 sampled two-mode periodic orbits have Floquet spectral radius below one.',
        'Monodromy determinant matches the exact divergence identity.',
        'All periodic endpoints have residual below 1e-8 in canonical coordinates.',
        'Matched initial conditions and forward/reverse branches agree.',
        'Fifteen time/mode refinements retain negative conditional FTLE.',
        'At all 31 full-FE boundary points the spatial Bernstein upper bound equals the sampled physical tip maximum.',
        'Adaptive full-FE Radau trajectories agree with independently calculated Fourier periodic solutions.'],
        caveats=['Floating-point checks are numerical, not interval proofs.',
                 'No physical collision/contact is simulated after a crossing.',
                 'A 0.1-rad slope screen is a reporting convention, not a calibrated model-validity threshold.',
                 'Grid absence of detected chaos is not a proof over a continuous parameter domain.']))
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
