"""Build a new citation-revision attachment without replacing earlier archives.

Run the existing compile_documents.py and verify_citation_revision_20261008.py
first. Reuses the original payload selection without running its package builder.
"""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent.parent
NAME = "S2_citation_revision_20261008.zip"


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    target = ROOT/NAME
    temporary = ROOT/(NAME+".tmp")
    if target.exists() or temporary.exists():
        raise FileExistsError("Refusing to overwrite an existing archive or pending build: " + str(target))
    qa = json.loads((ROOT/"verification/citation_revision_20261008.json").read_text(encoding="utf-8"))
    assert qa["status"] == "PASS"
    for name, expected in qa["source_sha256"].items():
        assert sha256(ROOT/name) == expected, name
    for name, record in qa["documents"].items():
        assert sha256(ROOT/(name+".pdf")) == record["pdf_sha256"], name
    old_builder = ROOT/"reproducibility/reviewer_revision_20261008/build_revision_package.py"
    spec = importlib.util.spec_from_file_location("prior_package_selection", old_builder)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    paths = set(module.selected_files())
    for name in ["核对说明/引文适用性修订说明_20261008.md", "verification/citation_revision_20261008.json",
                 "verification/citation_revision_20261008.diff", "verification/verify_citation_revision_20261008.py",
                 "verification/build_citation_package_20261008.py", "main.log", "supplement.log"]:
        paths.add(ROOT/name)
    paths = sorted(paths, key=lambda p:p.relative_to(ROOT).as_posix())
    assert all(p.is_file() and p.resolve().is_relative_to(ROOT.resolve()) for p in paths)
    files = [{"path":p.relative_to(ROOT).as_posix(), "bytes":p.stat().st_size, "sha256":sha256(p)} for p in paths]
    manifest = {"revision_date":"2026-10-08", "archive":NAME,
                "scope":"Current citation revision and preserved earlier scientific payload. S1 retains the original full scientific record.",
                "files":files, "manifest_excludes_itself":True}
    encoded = (json.dumps(manifest,ensure_ascii=False,indent=2)+"\n").encode("utf-8")
    manifest_name = "verification/citation_artifact_manifest_20261008.json"
    (ROOT/manifest_name).write_bytes(encoded)
    with ZipFile(temporary,"x",ZIP_DEFLATED,compresslevel=6) as archive:
        for path in paths:
            archive.write(path,path.relative_to(ROOT).as_posix())
        archive.writestr(manifest_name,encoded)
    with ZipFile(temporary) as archive:
        assert archive.testzip() is None, "Archive CRC failure"
        assert len(archive.namelist()) == len(files)+1
        for item in files:
            payload=archive.read(item["path"])
            assert len(payload)==item["bytes"] and hashlib.sha256(payload).hexdigest()==item["sha256"],item["path"]
        assert archive.read(manifest_name)==encoded
    # On Windows rename also refuses to replace an existing destination.
    assert not target.exists()
    temporary.rename(target)
    digest = sha256(target)
    with target.with_suffix(".zip.sha256").open("x",encoding="utf-8",newline="\n") as handle:
        handle.write(digest+"  "+NAME+"\n")
    report={"archive":NAME,"bytes":target.stat().st_size,"sha256":digest,
            "file_count":len(files)+1,"crc_and_hashes_verified":True}
    (ROOT/"verification/citation_package_build_20261008.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False))


if __name__ == "__main__":
    main()
