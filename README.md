# v10 引文适用性修订稿

初版日期：2026-10-07；四项问题与引文适用性修订：2026-10-08。原 v9 和修订前 v10 均保留。

本次引文修订澄清复网络导数逼近的适用范围、补充磁应力到偶极模型的标准化约、说明 Haro 常数来源并锁定定理版本。详见 `核对说明/引文适用性修订说明_20261008.md` 和 `verification/citation_revision_20261008.json`。四项审阅问题的既有计算结果继续随包提供；其原验证记录保留为历史证据。

- `main.pdf`：英文主稿，15 页，含参考文献。
- `supplement.pdf`：英文补充，45 页。两份 PDF 放在同一目录，保留跨文件链接。
- `main.tex`、`supplement.tex`、`sections/`、`figures/`：可编辑稿件。新增边界图位于下述实验目录，也包含在 S2 中。
- `核对说明/四项审阅问题处理说明_20261008.md`：四项问题、具体处理、证据和适用范围。
- `verification/four_issue_revision_20261008.json`：此前四项问题修订的验证；`four_issue_revision_redline.html` 为相应源码红线差异。
- `reproducibility/ARTIFACT_INDEX.md`：S1/S2 身份、校验和结果到源文件映射。

此前四项问题修订实际补充了 1,715 组受迫响应扫描、完整 FE 周期清隙边界及启动复核，以及非零初始化/传感/执行误差下的清隙证明与完整节点模型仿真。新输入、代码和原始结果自包含于 `reproducibility/reviewer_revision_20261008/stability_scan/` 与 `robustness/`，各目录 README 给出复算命令及依赖。机器人相关文献来源记录在 `literature/`。

使用未修改的原 IJRR/SAGE 模板。已有 Python 与 pdflatex 环境可执行：

```text
python reproducibility/compile_documents.py
```

脚本先编译补充，再编译正文并检查未定义引用和越栏。无需重新运行实验即可编译交付的源文件。

本次交付包为 `S2_citation_revision_20261008.zip`，附 `.zip.sha256`。完整历史科学材料 S1 仍为同级父目录的 `IJRR_CHNN_KAM_导师意见凝练稿_20261002.zip`，SHA-256：`f335542f8eec42067bda474306ebdcb23dcca0450f132fcba7723621eab34170`。S1 含原 COMSOL 工程及全部历史计算记录，内部论文是 v9；阅读本次稿件用 S2。分享完整科学材料时同时提供 S1 和 S2。

本次引文修订未修改数据、重新训练拟合模型或重跑历史 COMSOL 工程。此前四项问题修订中的扫描未发现规定采样协议内的混沌转变；它提供数值共振/清隙边界。既有严格误差容限针对所列条件下 4 秒清隙，并不声称噪声下保持精确 KAM 环面。

四项问题修订前备份：`verification/backups/v10_pre_four_issue_revision_20261008.zip`，SHA-256：`e5a0531fcfa82a70626c2a5168f0982af8de40335094f5b8e10f72f83483b73a`。此前 `v10_*.md/json` 属于相应历史修订记录，其旧页数与旧哈希不代表本次结果。

本次引文修订前备份：`verification/backups/v10_pre_citation_revision_20261008.zip`，SHA-256：`7102fc942222037f3c4a212cfb95bb7bdcb6e4b666a8bba249069fbcc3388f5d`。原 `S2_reviewer_revision_20261008.zip` 保留为本次引文修订前的完整交付包。最新包内文件清单见 `verification/citation_artifact_manifest_20261008.json`；打包脚本为 `verification/build_citation_package_20261008.py`。

## GitHub 下载

稿件源码、PDF、复现数据与核对记录可直接从本仓库获取。两个原始 S2 交付包超过 GitHub 普通文件大小限制，保留原始内容及校验和，作为 [v10-20261008 Release 附件](https://github.com/mop74720-commits/IJRR_CHNN_KAM_v10_20261007/releases/tag/v10-20261008) 提供：

- [最新 S2 引文修订包](https://github.com/mop74720-commits/IJRR_CHNN_KAM_v10_20261007/releases/download/v10-20261008/S2_citation_revision_20261008.zip)
- [此前 S2 审阅修订包](https://github.com/mop74720-commits/IJRR_CHNN_KAM_v10_20261007/releases/download/v10-20261008/S2_reviewer_revision_20261008.zip)

下载后将 ZIP 放在仓库根目录，可供需要原始交付包的验证脚本使用。对应 `.zip.sha256` 文件同时保留在仓库和 Release 附件中。编译临时文件和 Python/Numba 缓存不纳入版本控制。

本仓库发布范围为 `IJRR_CHNN_KAM_v10_20261007` 目录；上述位于父目录的 S1 完整历史材料包未随本仓库发布。完整复现原始研究仍需 S1，详见 `reproducibility/ARTIFACT_INDEX.md`。
