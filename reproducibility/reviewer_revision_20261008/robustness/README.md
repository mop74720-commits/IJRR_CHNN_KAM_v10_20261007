# Finite initialization and port-error tolerances

This self-contained supplement reproduces the finite-horizon initialization,
sensing and execution-error budget in `sections/operational_tolerances.tex`
(label `sec:operational-tolerances`). It uses the full 64-DOF frozen nodal
beam matrices. The archived CHNN only supplies the requested commands.
The robust work governor can change that request; its clearance guarantee
does not assert exact invariant-torus motion after activation.

## Reproduction

Use Python 3.12 and the versions in `requirements.txt`. From this directory:

```text
python run_operational_tolerances.py --stage all
```

The script needs no historical archive or absolute source path. It reads
only `inputs/`. Stages `certificate`, `primary`, `refined`, and `compare`
can also be run separately. Every stage regenerates the analytic certificate
and the environment record. `compare` requires the two trajectory files.
A complete run took approximately 202 seconds on the development machine;
timing is not part of the scientific result. Numba may create the local
`numba_cache/` directory, which is not a scientific input and can be omitted
from the distribution.

The original local validation used python-flint 0.9.0 from the historical
archive's vendor directory via a process-local `PYTHONPATH`. That directory
is not referenced by the delivered code. A normal installation of the listed
dependencies supplies the same import. `environment.json` records the exact
versions actually used. No global package or system setting was changed.

## Inputs and scientific object

- `model_matrices.npz`: original frozen nodal M, K and all 64 modal columns.
- `integral_fit.npz`: original COMSOL integral-fit learned coefficients.
- `feature_model.py`: original learned-energy derivative implementation,
  copied unchanged, with its original file hash and line endings preserved.
- `finite_fourier_initial_state.npz`: only the phase-(0,0) state and
  frequency data extracted from the finite Fourier candidate.
- `reference_nodal_state.npz`: the frozen dyadic nodal position and velocity
  reconstruction used directly in every initialization bound.
- `manifest.json`: original archive-relative locations, SHA-256 hashes,
  and definitions of the two small derived inputs.

Hashes of every delivered input are also recorded in `environment.json`.
No COMSOL native project or full trajectory source archive is required.
The three original inputs are copied byte-for-byte. The reference is a finite
state, not a certified exact-torus member. The full nodal mechanical plant
has small diagonalization differences from the separately frozen diagonal
modal plant in the KAM theorem; the simulation retains the full transformed
mass and stiffness matrices instead of dropping those differences.

## Analytic guarantee

The proof is for the exact exported nodal matrices, positive Rayleigh damping
alpha = 0.55558524 /s, and the stipulated deterministic error envelopes.
Arb uses 256 bits. It verifies positive LDL pivots for M and K, exact
Hermite Bernstein rows, matrix inverse compliances, and all reported scalar
bounds. Float exports are moved outwards with `nextafter`; internal ball
calculations are retained. The actual initialized error is checked with
Arb products Phi*a0 - q_reference and Phi*v0 - v_reference, thereby including
the coordinate-solve and reconstruction rounding of the exported initial state.

The initial full-state error ellipsoid is Q(error) <= 5 pJ; the nominal
quadratic energy is below 21.45855219 nJ. This implies initial physical
energy below 22.11866366 nJ. Its tip-error projection is 56.006 micrometres;
this is not an independent coordinate box, nor a sufficient tip-only test.

For 4 seconds, the specified envelopes are angle 1 mrad, angular rate
0.01 rad/s, tip velocity 0.1 mm/s, gain error 1% on both ports, additive
force 0.5 nN, and additive field 0.05 microtesla. The additive-work reserve
is below 0.09677917 nJ. The 22.66 nJ cap guarantees centreline displacement
below 3.898405 mm against a 4 mm net corridor allowance. Command limits
reserve both relative and additive errors before clipping:

```text
F_command_max = (20e-6 N - 0.5e-9 N) / 1.01
B_command_max = (0.02 T - 0.05e-6 T) / 1.01
```

Thus the stated actual limits remain 20 micronewtons and 20 mT under the
specified envelopes. The force limit corresponds to 1 T/m for mu=2e-5 A m^2.
With constant nonzero additive errors this claim has the stated 4 s horizon;
it is not an indefinite robustness claim. The section also gives sufficient
sample-age and delay/slew conditions. They are conditional formulas, not a
numerical hardware-latency validation. The demonstrated sensors and errors
are deterministic test inputs, not measured hardware specifications.

## Numerical check and outputs

The nominal KAM-realizing feedback uses the original fitted derivatives and
a velocity-to-momentum solve. The request receives a bounded low-state
perturbation. Distinct sinusoidal sensor, gain and additive-port errors act
at their specified amplitudes. The deterministic initialization has verified
error energy below 4.95000001 pJ. All 64 modes and the energy/work ledger are
integrated by Radau. `primary` and `refined` reduce both relative and absolute
tolerances by a factor of ten; 2001 common times are saved.

- `analytic_certificate.json`: validated energy, compliance and error-set bounds.
- `primary_results.json`, `refined_results.json`: implementation checks and extrema.
- `primary_trajectory.npz`, `refined_trajectory.npz`: states, ports, work,
  dissipation, Bernstein bounds and ledger data.
- `refinement_comparison.json`: differences at the saved common times.
- `environment.json`: actual runtime and dependency versions plus input hashes.

The refined maximum energy is 21.7448979265 nJ; the saved-time centreline
bound is 1.9441442281 mm. Actual field and force maxima are 0.2167243896 mT
and 0.1815663330 micronewtons. The governor is active at 97.8011% of saved
times. The refined energy-ledger residual is 2.4803e-16 J; refinement changes
the tip by 9.6951e-9 m and energy by 3.1596e-15 J. These sampled numerical
results are not a validated ODE enclosure. Continuous-time clearance follows
from the analytic work theorem under its hypotheses, not interpolation of
saved times. The high activation fraction records substantial filtering of
the nominal KAM request and is not evidence of retained exact KAM motion.
