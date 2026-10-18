"""Export byte-preserving subsets of the frozen S1 archive for GitHub Releases."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


S1_SHA256 = "f335542f8eec42067bda474306ebdcb23dcca0450f132fcba7723621eab34170"
SUBSETS = {
    "S1_comsol_20261002.zip": "unified_model/comsol/",
    "S1_integral_certificate_20261002.zip": "mainline_completion/integral_certificate/",
    "S1_comsol_identification_refined_20261002.zip": "unified_model/learning/raw/comsol_identification_refined/",
    "S1_theory_20261002.zip": "unified_model/theory/",
    "S1_s_path_20261002.zip": "mainline_completion/s_path/",
}


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path, help="New output directory; must not exist")
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    assert sha256(source) == S1_SHA256, "Source is not the frozen S1 archive"
    output = args.output.resolve()
    output.mkdir(parents=False, exist_ok=False)
    manifest = {
        "source_filename": source.name,
        "source_bytes": source.stat().st_size,
        "source_sha256": S1_SHA256,
        "scope": "Original S1 and exact directory subsets; no scientific results regenerated",
        "archives": [],
    }
    complete = output / "S1_complete_original_20261002.zip"
    with source.open("rb") as src, complete.open("xb") as dst:
        shutil.copyfileobj(src, dst)
    assert sha256(complete) == S1_SHA256
    with ZipFile(source) as original:
        entries = [entry for entry in original.infolist() if not entry.is_dir()]
        assert len({entry.filename for entry in entries}) == len(entries)
        roots = {entry.filename.split("/")[0] for entry in entries}
        assert len(roots) == 1
        prefix = next(iter(roots)) + "/"
        manifest["original_archive_root"] = prefix
        assert original.testzip() is None, "Original ZIP integrity check failed"
        manifest["archives"].append({
            "filename": complete.name,
            "bytes": complete.stat().st_size,
            "sha256": S1_SHA256,
            "file_count": len(entries),
            "description": "Byte-identical copy of the complete original S1 ZIP",
        })
        for filename, relative in SUBSETS.items():
            selected = [e for e in entries if e.filename.startswith(prefix + "reproducibility/" + relative)]
            assert selected, relative
            payload = []
            path = output / filename
            with ZipFile(path, "x", compression=ZIP_DEFLATED, compresslevel=6) as subset:
                for entry in selected:
                    name = entry.filename[len(prefix):]
                    assert ".." not in PurePosixPath(name).parts
                    data = original.read(entry)
                    info = ZipInfo(name, date_time=entry.date_time)
                    info.compress_type = ZIP_DEFLATED
                    info.external_attr = entry.external_attr
                    subset.writestr(info, data)
                    payload.append({
                        "path": name,
                        "source_entry": entry.filename,
                        "bytes": len(data),
                        "sha256": hashlib.sha256(data).hexdigest(),
                    })
            with ZipFile(path) as subset:
                assert subset.namelist() == [item["path"] for item in payload]
                for item in payload:
                    assert hashlib.sha256(subset.read(item["path"])).hexdigest() == item["sha256"]
            record = {
                "filename": filename,
                "source_directory": "reproducibility/" + relative,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "file_count": len(payload),
                "files": payload,
            }
            manifest["archives"].append(record)
            print(json.dumps({key: value for key, value in record.items() if key != "files"}), flush=True)
    checksums = []
    for record in manifest["archives"]:
        checksum = record["sha256"] + "  " + record["filename"] + "\n"
        (output / (record["filename"] + ".sha256")).write_text(checksum, encoding="utf-8", newline="\n")
        checksums.append(checksum)
    (output / "S1_SHA256SUMS.txt").write_text("".join(checksums), encoding="utf-8", newline="\n")
    (output / "S1_export_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    readme = """# S1 原始材料与目录分包

`S1_complete_original_20261002.zip` 是原始 S1 的逐字节副本，仅使用便于下载的英文文件名。原始名称与 SHA-256 见 `S1_export_manifest.json`。其中的历史论文为 v9；v10 正文和补充请使用同一 Release 的 S2。

五个目录分包逐文件取自同一 S1，未修改源码、工程、数据或历史记录，ZIP 内从 `reproducibility/` 开始保留路径。将所需分包解压到同一个工作目录，可合并原有目录结构；完整 S1 保留原始外层目录。

这些分包是原始目录切片，未额外补入跨目录依赖，不能据此声称每个分包可独立复算。需要完整依赖时使用完整 S1，并按包内 README 准备环境。重新求解原 COMSOL 工程需要相应许可、COMSOL 及 LiveLink for MATLAB。此次仅完成打包与完整性核验，未重跑 COMSOL、训练或数学认证。

每个 ZIP 均附 `.sha256`，汇总见 `S1_SHA256SUMS.txt`。`S1_export_manifest.json` 列出每个分包的原始路径、文件大小及逐文件 SHA-256。

"""
    for name, relative in SUBSETS.items():
        readme += f"- `{name}` → `reproducibility/{relative}`\n"
    (output / "S1_DOWNLOAD_README.md").write_text(readme, encoding="utf-8", newline="\n")
    print(json.dumps({"output": str(output), "archives": len(manifest["archives"]), "verified": True}), flush=True)


if __name__ == "__main__":
    main()
