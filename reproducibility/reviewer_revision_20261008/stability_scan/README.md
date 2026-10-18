# 同一磁驱动梁的受迫响应、共振及清隙边界

本目录新增于 2026-10-08。只复用冻结的原模型矩阵，没有修改旧数据或拟合权重。
研究对象为现稿同一 32 单元磁驱动梁，静偏场 30 mT、质量阻尼
0.5555852368222173 s^-1、完整磁正弦/余弦。扫描没有 governor、磁梯度或接触力。
全部新实验真正运行；不是将既有图改名。英文补充为稿件的 `sections/stability_boundary.tex`。

## 结论与边界

- 7 个幅值、49 个频率比、3 个初态给出 1029 次独立积分；升降频延拓再给出 686 次。
  在规定网格均收敛到同驱动周期响应，Floquet 最大模 0.803809--0.889055；未检测到混沌或多稳态。
- 这是受迫阻尼模型的有限网格局部稳定性证据，不是连续参数域的无混沌定理，亦不是保守 KAM 环面破裂阈值。
- 主网格以两模态进行，15 个时间加密/四模态/八模态检查保留周期性和负条件有限时间 Lyapunov 指数。
  完整 64 模态周期响应另用精确线性模态 Fourier 传递及非线性固定点求解，直接检验重构后的非线性 ODE 残差。
- 全 64 模态稳态 4 mm 轮廓在近共振处的最小场幅值约 1.915911 mT，频率比约 0.999245。
  这条曲线所有根处的空间 Bernstein 上界等于实际 tip 采样峰值；一般 Bernstein 上界本身只是空间充分界。
- 启动轮廓依赖零初态及 0.25 s 指数建立的正弦场，不能当作任意初态的安全阈值。
  八模态启动根经过加密及完整 FE 抽核；后者峰值差不超过 2.44 µm。
- 0.1 rad 是预先声明的小斜率报告筛选值，不是经标定的力学有效性定理。域外计算保留但不作物理安全结论。
  原输入功 governor 的严格保证与这里的数值响应轮廓分别成立。

## 输入和可移植性

`inputs/reference_model_matrices.npz` 是 v9 完整科学资料中
`reproducibility/unified_model/learning/raw/model_matrices.npz` 的精确副本，SHA-256：
`ae1b221b5f67ff455c478c17552bc26ce4d8c155169b9e2682f53b992cca4288`。
`inputs/reference_model_spec.json` 为同目录模型说明的副本。
代码只从自身目录读取这些输入，不依赖旧 E 盘路径、COMSOL、MATLAB或历史 ZIP。
说明文件中的旧来源路径仅是历史元数据。

实测环境：Python 3.12.7、NumPy 1.26.4、SciPy 1.13.1、Numba 0.60.0、Matplotlib、
threadpoolctl；绘图预览另外使用 PyMuPDF。完整版本和输入/输出哈希见 `manifest.json`。
脚本的数值库线程设置仅对本进程生效。新文本使用 UTF-8、CRLF。

## 复算顺序

在本目录执行：

```text
python run_stability_scan.py grid --threads 16
python run_stability_scan.py continuation --threads 14
python verify_and_boundaries.py checks
python verify_and_boundaries.py floquet --threads 12
python verify_and_boundaries.py refine --threads 12
python verify_and_boundaries.py roots
python verify_and_boundaries.py minimum
python verify_and_boundaries.py startup
python summarize_scan.py
python build_manifest.py
```

这些命令只会重写本目录自身的相应实验输出，请在新复制的目录中复算以保留交付结果。
主网格与延拓在记录机器上各约 33/26 秒，模态和时间细化约 15 秒；完整启动检查约 153 秒。
首次 Numba 编译增加耗时。最后几位可能随环境变化，记录值不是区间严格误差界。
`pilot*` 为六个早期试运行，只用于性能/范围预检，不进入任何正式表格或结论。

## 证据入口

- `protocol.json`：在主网格前冻结的参数、初态、时间和判据；`protocol_amendment.json`：后续全 FE 边界求解与复核方法。
- `grid_raw.npz`、`continuation_raw.npz`：参数、每周期状态、位移/斜率峰值、局部指数、末周期完整波形。
- `*_floquet.npz`：每个周期解的完整 monodromy、复乘子、周期余差和行列式检验。
- `full_FE_boundary_waveforms.npz`：31 个完整 FE 稳态根的 1024 相位、128 维真实状态。
- `boundary_modal_phase_refinement.json`、`full_FE_Radau_r*.npz`：模态/相位和独立积分复核。
- `startup_*mode_r*.npz`、`full_FE_startup_r*.npz`：固定 ramp/零初态的启动根和完整 FE 响应。
- `summary.json`、`validation.json`：可由原始结果重新生成的汇总和断言；不替代原始数据。
- `fig_stability_boundaries.pdf`：矢量结果图；同名 `preview.png` 仅用于布局检查。

最大 Lyapunov 指数是固定外部驱动相位的条件有限时间指数，未把外加时钟的零指数计入。
Floquet 分析目前覆盖完整两模态网格；四/八模态做所列 15 例积分检查；完整 64 模态做稳态
非线性周期解、3 例周期 Radau 对照以及 3 例启动 Radau 抽核。不得把其中一种范围写成另一种。
