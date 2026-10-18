"""Verify the citation-only revision without rerunning scientific experiments.

Requires SymPy/PyMuPDF, the preserved local backup and the old S2 archive.
Run from any working directory. Writes current QA and a textual revision diff.
"""
from pathlib import Path
from zipfile import ZipFile
import difflib
import hashlib
import json
import re
import sys

import fitz
import sympy as sp

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent.parent
BACKUP = ROOT / "verification/backups/v10_pre_citation_revision_20261008.zip"
OLD_S2 = ROOT / "S2_reviewer_revision_20261008.zip"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    assert sha256(BACKUP) == "7102fc942222037f3c4a212cfb95bb7bdcb6e4b666a8bba249069fbcc3388f5d"
    assert sha256(OLD_S2) == "5fbb87da73c68b729b5d452aea36bad57a24e25c72721da9307e1972fa9f93a4"
    m = sp.symbols("m0:3", real=True)
    b = sp.symbols("b0:3", real=True)
    g = sp.Matrix(3, 3, sp.symbols("g0:9", real=True))
    force = []
    torque = []
    for i in range(3):
        contracted = sum(sp.LeviCivita(i, j, k) * sp.LeviCivita(l, j, c) * m[c] * g[k, l]
                         for j in range(3) for k in range(3) for l in range(3) for c in range(3))
        expected = sum(m[c] * g[c, i] for c in range(3)) - m[i] * sp.trace(g)
        force.append(sp.simplify(contracted - expected))
        contracted_t = sum(sp.LeviCivita(l, i, c) * m[c] * b[l] for l in range(3) for c in range(3))
        torque.append(sp.simplify(contracted_t - sp.Matrix(m).cross(sp.Matrix(b))[i]))
    aq, ap, q, p, bias = sp.symbols("A_Q A_P Q P bias", real=True)
    actual = sp.sqrt(2) * (aq - sp.I * ap) * (q + sp.I * p) / sp.sqrt(2) + bias
    expected = aq*q + ap*p + bias + sp.I*(aq*p-ap*q)
    complex_difference = sp.simplify(actual - expected)
    assert all(value == 0 for value in force + torque + [complex_difference])

    templates = {}
    changes = []
    diff = []
    with ZipFile(BACKUP) as archive:
        assert archive.testzip() is None
        for name in archive.namelist():
            if name.startswith("workspace_indexes/"):
                basename = Path(name).name
                path = ROOT.parents[1] / basename if basename.startswith("论文") else ROOT.parent / basename
            else:
                path = ROOT / name
            before = archive.read(name)
            after = path.read_bytes()
            if name in {"sagej.cls", "SageH.bst", "SageV.bst"}:
                assert before == after, name
                templates[name] = sha256(path)
            if before != after and path.suffix in {".tex", ".md"}:
                changes.append(name)
                diff.extend(difflib.unified_diff(before.decode("utf-8-sig").splitlines(True),
                                                after.decode("utf-8-sig").splitlines(True),
                                                fromfile="before/"+name, tofile="after/"+name))
    (ROOT / "verification/citation_revision_20261008.diff").write_text("".join(diff), encoding="utf-8", newline="")

    scientific_files = []
    with ZipFile(OLD_S2) as archive:
        for name in archive.namelist():
            if not name.startswith("reproducibility/reviewer_revision_20261008/") or name.endswith("/"):
                continue
            if Path(name).name in {"artifact_manifest.json", "package_build.json"}:
                continue
            assert archive.read(name) == (ROOT/name).read_bytes(), name
            scientific_files.append(name)

    documents = {}
    all_tex = set()
    for base in ("main", "supplement"):
        pending = [ROOT / (base+".tex")]
        seen = set()
        while pending:
            path = pending.pop()
            if path in seen:
                continue
            assert path.is_file(), path
            seen.add(path)
            source = re.sub(r"(?<!\\)%[^\n]*", "", path.read_text(encoding="utf-8"))
            for include in re.findall(r"\\(?:input|include)\{([^}]+)\}", source):
                child = ROOT / include
                pending.append(child if child.suffix else child.with_suffix(".tex"))
        all_tex.update(seen)
        log = (ROOT/(base+".log")).read_text(encoding="utf-8", errors="replace")
        failures = re.findall(r"^.*(?:Overfull|undefined|LaTeX Error|Rerun to get|Label\(s\) may have changed).*", log, re.M)
        assert not failures, failures
        doc = fitz.open(ROOT/(base+".pdf"))
        expected_pages = 15 if base == "main" else 45
        assert len(doc) == expected_pages
        outside = []
        for i, page in enumerate(doc):
            for word in page.get_text("words"):
                if not page.rect.contains(fitz.Rect(word[:4])):
                    outside.append({"page": i+1, "word": word[4]})
        assert not outside, outside
        documents[base] = {"pages": len(doc), "pdf_sha256": sha256(ROOT/(base+".pdf")),
                           "input_files": sorted(p.relative_to(ROOT).as_posix() for p in seen),
                           "unresolved_reference_or_overfull_or_error_messages": failures,
                           "words_outside_page": len(outside), "log_sha256": sha256(ROOT/(base+".log"))}
        doc.close()

    bibliography = (ROOT/"sections/references.tex").read_text(encoding="utf-8")
    keys = re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", bibliography)
    assert len(keys) == len(set(keys)) == 24
    cited = set()
    for path in all_tex:
        source = re.sub(r"(?<!\\)%[^\n]*", "", path.read_text(encoding="utf-8"))
        for group in re.findall(r"\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}", source):
            cited.update(key.strip() for key in group.split(","))
    assert cited == set(keys), {"missing": sorted(cited-set(keys)), "uncited": sorted(set(keys)-cited)}
    report = {
        "revision_date": "2026-10-08", "status": "PASS",
        "scope": "Implemented citation corrections, document consistency and packaging inputs; not a blanket full-text audit or new experimental replication.",
        "backup": {"file": BACKUP.relative_to(ROOT).as_posix(), "sha256": sha256(BACKUP), "crc_verified": True},
        "preserved_old_S2": {"file": OLD_S2.name, "sha256": sha256(OLD_S2)},
        "changed_existing_text_files": changes,
        "source_sha256": {p.relative_to(ROOT).as_posix(): sha256(p) for p in sorted(all_tex)},
        "template_files_unchanged_sha256": templates,
        "documents": documents,
        "citations": {"bibliography_count": len(keys), "cited_keys": sorted(cited), "missing": [], "uncited": []},
        "symbolic_checks": {"scalar_identities_passed": 7, "force_component_differences": list(map(str, force)),
                            "torque_component_differences": list(map(str, torque)), "complex_preactivation_difference": str(complex_difference),
                            "comparison": "simplify(lhs-rhs)==0; replaces a prior structural-equality false alarm"},
        "scientific_files_unchanged_vs_old_S2": {"count": len(scientific_files), "files": scientific_files},
        "visual_review": {"supplement_pages": [1,2,6,11,40,45], "scope": "Edited passages, plus automated page-bound checks on all pages; not a fresh visual review of every page."},
        "independent_review": ["Complex embedding, theorem scope, version locators and least-squares distinction: PASS",
                               "Magnetic reduction, signs, remainder dimensions and Haro Case III constant: PASS"],
        "limits": ["Antman and Simo original full texts not fully obtained/verified.",
                   "Hairer checked against author lectures, not the entire monograph.",
                   "No retraining, COMSOL rerun, or rerun of all historical 512-bit interval certificates."]
    }
    target = ROOT/"verification/citation_revision_20261008.json"
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "pages": [documents[x]["pages"] for x in ("main","supplement")],
                      "citations":len(keys), "symbolic_identities":7, "unchanged_scientific_files":len(scientific_files),
                      "changed_text_files":changes}, ensure_ascii=False))


if __name__ == "__main__":
    main()
