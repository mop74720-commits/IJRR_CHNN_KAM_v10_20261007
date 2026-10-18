# Computational attachments and result-to-source map

Revision date: 8 October 2026. The article uses two separately identified local supplementary archives. Deliver both archives with the manuscript when sharing the complete computational record.

## Archive S1: complete original scientific record

- Supplied filename: `IJRR_CHNN_KAM_导师意见凝练稿_20261002.zip`.
- Size: 808,408,246 bytes.
- SHA-256: `f335542f8eec42067bda474306ebdcb23dcca0450f132fcba7723621eab34170`.
- ZIP root: `IJRR_CHNN_KAM_导师意见凝练稿_20261002/`.
- The archive contains native COMSOL models, their construction scripts, original training data, fitted weights, interval certificates, and the original learning, integration and control experiments. Its internal historical manuscript is v9; the revised article and Supplement are supplied in S2.
- The working extraction of this archive is named `IJRR_CHNN_KAM_主线补齐稿_20261002/`. The extraction-directory name is different; the scientific inputs are identified by the recorded file hashes.

Paths in the next table are relative to S1's `reproducibility/` directory.

| Result | Source directory or entry | Object and interpretation |
|---|---|---|
| Matched real/complex and prior comparisons | `unified_model/learning/README_中文.md`, `run_learning.py`, `raw/` | Three paired seeds; noiseless states with specified derivative noise; 35/45 effective real parameters |
| Independent FE construction and response | `unified_model/comsol/README_复现说明.md` | Native COMSOL beam, matrices, modal map and independent time integration |
| State-only learning-to-KAM certificate | `mainline_completion/integral_certificate/` | Two-mode COMSOL integral fit; this is the fitted Hamiltonian used for the main error-to-KAM proof |
| Full-FE data-driven predictor | `unified_model/learning/raw/comsol_identification_refined/` | Separate full-FE position/velocity integral fit, used for prediction and S-path control |
| Finite-amplitude torus certificates | `unified_model/theory/README_认证说明.md`, `verify_torus.py` | Frozen frequencies, magnetic coefficients and complete quartic readouts; Fourier candidates and analytic-tail validation |
| Central4 convergence | `unified_model/integrator/` | GC17 energy; same fitted function across the compared integrators |
| Full-model magnetic realization | `mainline_completion/resonance/` | Certified two-mode motion lifted to the controlled 64-mode model; online COMSOL refinement |
| Physical work and centreline clearance | `unified_model/safety/README.md` | Full physical energy and its actual magnetic work port |
| S-shaped endpoint task | `mainline_completion/s_path/` | Full-FE-trained controller, axial base feed, measured tip feedback and independent field replay |
| Local measure and finite-time classification | `mainline_completion/measure_completion/README.md` | Exact local KAM tube and separately defined finite-time classifier; their probabilities are different objects |

## Archive S2: revised article and additional evidence

- Supplied filename: `S2_citation_revision_20261008.zip`.
- The companion `S2_citation_revision_20261008.zip.sha256` records the archive checksum. Inside the ZIP, `verification/citation_artifact_manifest_20261008.json` lists every archived payload file's size and SHA-256.
- This citation revision retains the earlier scientific computations. The previous `S2_reviewer_revision_20261008.zip` remains the frozen pre-citation attachment.
- Extract S2 into a new directory. Its root contains `main.tex`, `supplement.tex`, the corresponding PDFs, the unmodified SAGE class, `sections/`, `figures/`, and the additional experiments below.
- The new experiments carry their minimal frozen inputs locally. They do not require absolute paths to the original author's workstation or a new COMSOL solve. S1 remains necessary to rerun the original study in full.

Paths in this table are relative to S2's `reproducibility/reviewer_revision_20261008/` directory.

| Revision issue | Evidence location | What can be checked |
|---|---|---|
| Forced-response and operational boundary | `stability_scan/` | Frozen protocol, input hashes, independent initial states, bidirectional continuation, periodic-response stability diagnostics and numerical refinement |
| Nonzero initialization and port errors | `robustness/` | Physical-energy tolerance calculation, uncertainty-aware work bounds, controller implementation and recorded full-model responses |
| Closest related work | `literature/related_work_patch.json` | Verified formal publications, exact inserted text, bibliography and claim-to-primary-source mapping |
| Citation applicability corrections | `verification/citation_revision_20261008.json` and `核对说明/引文适用性修订说明_20261008.md` in the archive root | Current source, mathematical corrections, compilation and unresolved original-text verification limits |
| Source, layout and result verification | `verification/four_issue_revision_20261008.json` in the archive root | Historical four-issue compilation and scientific result checks; current document hashes are in the citation-revision record |
| Revision correspondence | `核对说明/四项审阅问题处理说明_20261008.md` in the archive root | Four issues, concrete changes, evidence and remaining scientific scope |

See each experiment's README for its commands and recorded environment. Rebuild the article after the experiments are available with:

```text
python reproducibility/compile_documents.py
```

This compiles `supplement.tex` first and then `main.tex`, resolving the external supplementary references. Keep the two PDFs together. The class and bibliography styles match the supplied original IJRR/SAGE V1.20 template.

## Scope of access

These are supplied supplementary files, identified by filenames and hashes; no public repository DOI is claimed. The full scientific package comprises S1 and S2. The original v10 baseline is retained locally under `verification/backups/`, separately from the clean revised attachment.
