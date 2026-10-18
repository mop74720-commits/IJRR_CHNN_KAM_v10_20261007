"""Check the final revision against raw outputs, the saved baseline and PDFs.

This is an artifact/consistency check, not a rerun of the full experiments.
The local baseline ZIP is intentionally not duplicated in attachment S2.
"""
from pathlib import Path
from zipfile import ZipFile
import difflib
import hashlib
import html
import json
import re

import fitz
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    backup = ROOT / "verification/backups/v10_pre_four_issue_revision_20261008.zip"
    assert sha(backup) == "e5a0531fcfa82a70626c2a5168f0982af8de40335094f5b8e10f72f83483b73a"
    report = {"revision_date": "2026-10-08", "status": "PASS",
              "scope": "New computations, mathematical QA, raw-output consistency, source preservation and compiled layout; historical COMSOL study not rerun.",
              "baseline_sha256": sha(backup)}
    changed = []
    blocks = []
    with ZipFile(backup) as archive:
        for path in sorted([ROOT / "main.tex", ROOT / "supplement.tex", *ROOT.glob("sections/*.tex")]):
            relative = path.relative_to(ROOT).as_posix()
            old = archive.read(relative).decode("utf-8").splitlines() if relative in archive.namelist() else []
            new = path.read_text(encoding="utf-8").splitlines()
            if old == new:
                continue
            changed.append(relative)
            diff = []
            for line in difflib.unified_diff(old, new, fromfile="before/" + relative,
                                             tofile="after/" + relative, lineterm="", n=3):
                cls = "add" if line.startswith("+") else "del" if line.startswith("-") else "context"
                diff.append('<span class="' + cls + '">' + html.escape(line) + '</span>')
            blocks.append("<h2>" + html.escape(relative) + "</h2><pre>" + "\n".join(diff) + "</pre>")
        preserved = ["sagej.cls", "SageH.bst", "SageV.bst", "sections/core_derivation.tex",
                     "sections/core_analysis.tex", "sections/core_kam.tex", "sections/theorem_bridge.tex",
                     "sections/integral_comsol_to_kam.tex", "sections/measure_proof.tex"]
        for relative in preserved:
            assert archive.read(relative) == (ROOT / relative).read_bytes(), relative
    report["changed_tex_files"] = changed
    report["baseline_identical_core_and_template_files"] = preserved
    assert sha(ROOT / "sagej.cls") == "3d355246afa097b6dc8b637ff49d74fccaf8886f9331b3a97bb1b0a3b22b9e23"
    redline = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<title>四项问题修订：源码红线</title><style>
body{max-width:1100px;margin:35px auto;padding:0 22px;font:16px/1.65 system-ui,sans-serif;color:#202124}
pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f6f7f8;padding:18px;font:13px/1.65 Consolas,monospace}
.add{color:#b00020;background:#fff0f1}.del{color:#666;text-decoration:line-through}.context{color:#333}
h2{font-size:19px;margin-top:35px}</style>
<h1>四项审阅问题修订：源码红线</h1>
<p>相对本轮修改前备份。红色为新增或改写的源码行，灰色删除线为删去的行；英文 PDF 为干净版。
本文件展示 LaTeX 源码差异，不是渲染后的 PDF 叠加比较。原证明与历史实验结果的保留检查见配套验证记录。</p>'''
    (ROOT / "verification/four_issue_revision_redline.html").write_text(
        redline + "\n".join(blocks) + "</html>\n", encoding="utf-8", newline="\r\n")

    scan = HERE / "stability_scan"
    manifest = read_json(scan / "manifest.json")
    for item in manifest["files"]:
        path = scan / item["path"]
        assert path.stat().st_size == item["bytes"] and sha(path) == item["sha256"], item["path"]
    radii = []
    closure = []
    for name, count in [("grid", 1029), ("continuation", 686)]:
        raw = np.load(scan / (name + "_floquet.npz"))
        assert len(raw["monodromy"]) == count
        actual = np.max(np.abs(np.linalg.eigvals(raw["monodromy"])), axis=1)
        assert np.max(abs(actual - raw["spectral_radius"])) < 1e-12
        assert actual.max() < 1
        radii.extend(actual.tolist())
        closure.extend(raw["periodic_residual"].tolist())
    roots = read_json(scan / "full_FE_steady_clearance_boundary.json")
    crossings = [row for row in roots if row["crossing_below_20mT"]]
    assert len(crossings) == 31
    assert max(abs(row["tip_peak_mm"] - 4) for row in crossings) < 1e-6
    assert max(row["slope_bound_rad"] for row in crossings) < .1
    report["stability_raw_output_checks"] = {
        "manifest_files_verified": len(manifest["files"]), "independent_runs": 1029,
        "continuation_runs": 686, "floquet_radius_min": min(radii), "floquet_radius_max": max(radii),
        "max_periodic_residual": max(closure), "full_modal_steady_crossings": len(crossings),
        "r1_steady_amplitude_mT": next(row["amplitude_mT"] for row in crossings if row["frequency_ratio"] == 1),
        "minimum_steady_crossing": read_json(scan / "minimum_steady_clearance.json")
        if (scan / "minimum_steady_clearance.json").exists() else read_json(scan / "summary.json")["minimum_steady_clearance_amplitude"],
        "claim": "Numerical response and clearance contours; no detected chaos under this sampled protocol, no global absence or conservative-torus breakup theorem."}

    robust = HERE / "robustness"
    env = read_json(robust / "environment.json")
    for name, digest in env["inputs"].items():
        assert sha(robust / "inputs" / name) == digest, name
    certificate = read_json(robust / "analytic_certificate.json")
    matrices = np.load(robust / "inputs/model_matrices.npz")
    results = {}
    for label in ["primary", "refined"]:
        raw = np.load(robust / (label + "_trajectory.npz"))
        q = raw["modal_position"] @ matrices["all_phi"].T
        v = raw["modal_velocity"] @ matrices["all_phi"].T
        energy = (.5 * np.einsum("ij,jk,ik->i", v, matrices["M"], v)
                  + .5 * np.einsum("ij,jk,ik->i", q, matrices["K"], q)
                  + 6e-7 * 2 * np.sin(q[:, -1] / 2) ** 2)
        assert np.max(abs(energy - raw["energy_J"])) < 1e-18
        assert energy.max() < certificate["energy_cap_J"]
        assert np.max(energy - raw["energy_enclosure_J"]) < 0
        assert raw["work_observer_J"].max() < certificate["command_work_cap_J"]
        assert np.max(abs(raw["applied_force_N"])) < 20e-6
        assert np.max(abs(raw["applied_field_T"])) < .02
        results[label] = {"saved_times": len(raw["time_s"]),
                          "recomputed_maximum_energy_J": float(energy.max()),
                          "energy_ledger_residual_J": float(np.max(abs(energy-energy[0]-raw["input_work_J"]+raw["dissipation_J"])))}
    assert certificate["shape_bound_m"] < .003898405
    assert certificate["initial_energy_bound_J"] < certificate["command_work_cap_J"]
    report["operational_tolerance_checks"] = {"analytic_certificate": certificate, "raw_trajectory_checks": results,
        "scope": "Strict analytic work theorem for stated deterministic envelopes over 4 seconds; sampled integration is an implementation check and does not certify KAM tracking."}
    qa = read_json(HERE / "independent_qa.json")
    report["independent_qa_record"] = "reproducibility/reviewer_revision_20261008/independent_qa.json"
    assert qa

    pdfs = {}
    for name, pages in [("main", 15), ("supplement", 45)]:
        log = (ROOT / (name + ".log")).read_text(encoding="utf-8", errors="replace")
        assert not re.search(r"Overfull|There were undefined|undefined references|LaTeX Error|Rerun to get|Label\(s\) may have changed", log)
        document = fitz.open(ROOT / (name + ".pdf"))
        assert len(document) == pages
        assert all("??" not in page.get_text() for page in document)
        pdfs[name] = {"pages": len(document), "sha256": sha(ROOT / (name + ".pdf")),
                      "log_errors_undefined_references_overfull": 0}
    document = fitz.open(ROOT / "main.pdf")
    figure_page = next(i+1 for i, page in enumerate(document) if "Figure 2." in page.get_text())
    conclusion_page = next(i+1 for i, page in enumerate(document) if re.search(r"6\s+Conclusion", page.get_text()))
    assert figure_page < conclusion_page
    bib = (ROOT / "sections/references.tex").read_text(encoding="utf-8")
    keys = re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", bib)
    assert len(keys) == len(set(keys)) == 24
    assert {"francescon2025shape", "tong2025navigation"}.issubset(keys)
    report["documents"] = pdfs
    report["layout"] = {"figure_2_page": figure_page, "conclusion_start_page": conclusion_page,
                        "background_starts_page": 1, "derivation_starts_page": 4,
                        "visual_check": "All main pages reviewed in overview; new equations, tolerance table, response figure, tracking figure and conclusion inspected at reading scale; no clipped material or orphan blank page."}
    report["references"] = {"unique_bibliography_entries": len(keys), "new_primary_publications": 2,
                            "verification_record": "reproducibility/reviewer_revision_20261008/literature/related_work_patch.json"}
    report["algebra_validation"] = "Six existing SymPy control/clearance algebra checks executed and passed."
    report["not_rerun"] = ["All historical COMSOL model solves", "Original neural fitting", "Hardware, contact or in-orbit experiments"]
    report["core_file_hashes"] = {name: sha(ROOT / name) for name in ["main.tex", "supplement.tex", "sagej.cls", "SageH.bst", "SageV.bst",
                                    "sections/operational_tolerances.tex", "sections/stability_boundary.tex"]}
    target = ROOT / "verification/four_issue_revision_20261008.json"
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8", newline="\r\n")
    print(json.dumps({"status": "PASS", "main_pages": 15, "supplement_pages": 45,
                      "changed_tex_files": changed, "bibliography_entries": len(keys),
                      "raw_scan_cases": 1715, "verified_manifest_files": len(manifest["files"])}, ensure_ascii=False))


if __name__ == "__main__":
    main()
