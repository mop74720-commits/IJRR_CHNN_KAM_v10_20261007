"""Build the local revision attachment, with exact input/output checksums.

Run after compile_documents.py and the revision validation have passed.
The original 808 MB scientific archive is identified, not duplicated.
"""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
NAME = "S2_reviewer_revision_20261008.zip"


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def selected_files():
    names = [
        "main.tex", "main.pdf", "main.aux", "supplement.tex",
        "supplement.pdf", "supplement.aux", "sagej.cls", "SageH.bst",
        "SageV.bst", "README.md", "reproducibility/compile_documents.py",
        "reproducibility/ARTIFACT_INDEX.md",
        "verification/four_issue_revision_20261008.json",
        "verification/verify_control_algebra.py",
        "verification/core_control_algebra_validation.json",
        "verification/four_issue_revision_redline.html",
        "核对说明/四项审阅问题处理说明_20261008.md",
    ]
    paths = [ROOT / name for name in names]
    for folder in (ROOT / "sections", ROOT / "figures", HERE):
        for path in folder.rglob("*"):
            if not path.is_file():
                continue
            relative = path.relative_to(ROOT)
            if any(part in {"__pycache__", "numba_cache"} for part in relative.parts) or path.suffix in {".pyc", ".nbc", ".nbi"}:
                continue
            if path.name in {"artifact_manifest.json", "package_build.json"}:
                continue
            if path.name.startswith(("tmp_", "pilot")) or path.suffix == ".tmp":
                continue
            paths.append(path)
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required completed artifacts missing: " + repr(missing))
    return sorted(set(paths), key=lambda path: path.relative_to(ROOT).as_posix())


def main():
    paths = selected_files()
    files = [{"path": p.relative_to(ROOT).as_posix(), "bytes": p.stat().st_size,
              "sha256": sha256(p)} for p in paths]
    manifest = {"revision_date": "2026-10-08", "archive": "S2", "files": files,
                "scope": "Revised manuscript and new experiments; S1 retains the original complete scientific record."}
    encoded = (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    (HERE / "artifact_manifest.json").write_bytes(encoded)
    target = ROOT / NAME
    temporary = target.with_suffix(".zip.tmp")
    with ZipFile(temporary, "w", ZIP_DEFLATED, compresslevel=6) as archive:
        for path in paths:
            archive.write(path, path.relative_to(ROOT).as_posix())
        archive.writestr("reproducibility/reviewer_revision_20261008/artifact_manifest.json", encoded)
    with ZipFile(temporary) as archive:
        corrupt = archive.testzip()
        if corrupt:
            raise RuntimeError("Archive CRC failure: " + corrupt)
        for item in files:
            actual = hashlib.sha256(archive.read(item["path"])).hexdigest()
            if actual != item["sha256"]:
                raise RuntimeError("Archived content changed: " + item["path"])
    temporary.replace(target)
    digest = sha256(target)
    target.with_suffix(".zip.sha256").write_text(digest + "  " + NAME + "\n", encoding="utf-8")
    report = {"archive": NAME, "bytes": target.stat().st_size, "sha256": digest,
              "file_count": len(files) + 1, "crc_and_hashes_verified": True}
    (HERE / "package_build.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    main()
