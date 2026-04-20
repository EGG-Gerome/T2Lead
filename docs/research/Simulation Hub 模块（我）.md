# Simulation Hub 模块

### 要求上游提交的

**清理好的靶点 PDB (`protein.pdb`)**

- **不要全套：** 如果是多聚体，只保留药物结合的那个**单体/结构域**。
- **必须干净：** 删掉所有无关的结晶水、原本的配体和硫酸根等杂质。

**物理合规的配体 SDF (`ligand.sdf`)**

- **必须 3D：** 必须是已经对接（Docking）到口袋中的三维坐标，绝对不要 `SMILES`。
- **自带氢原子：** 必须在 pH 7.4 环境下加好氢原子。
- **无奇异元素：** 只允许 C, H, O, N, S, P, F, Cl, Br, I，禁止任何过渡金属。

**配置文件 JSON (`meta.json`)**

- 记录这个分子的 AI 预测打分、靶点名称等，方便你跑完 MD 后把数据合并。

---

1. **数据流向图：** BioDriver 如何把分子传给 Verification Hub？用 JSON 还是 Protocol Buffers？ 2.  **错误处理机制：** 如果 AM1-BCC 像你最近遇到的那样失败了，整个流水线是直接 Crash 还是抛出异常并记录 Log？ 3.  **开发标准：** 规定团队必须用 OpenMM 的哪个版本，必须写什么样的单元测试。？4.向华大基因要数据，把产品给广科医？

2. 我要不要用 GitHub Organization管理（单参考对吗）然后我负责**Simulation Hub** (物理仿真/验证中台)因为本来的T2Lead的核心就是Stage4，Dock+MD+RMSD的计算和验证？

    Stage 4 核心是处理复杂依赖、调用 OpenMM 跑隐式溶剂，以及极其关键的 **AM1-BCC 电荷计算**。

    - **目前的接收状态：** 它接收的就是上述的蛋白和配体坐标，然后分配力场（如 AMBER 力场给蛋白，GAFF/OpenFF 给配体）。
    - **你的重头戏（升级计划）：**
        1. **引入显式溶剂 (Explicit Solvent)：** 隐式溶剂算得快但不准（尤其是氢键网络）。你要在代码里加入添加水盒子（Water Box）和抗衡离子（Ions）的逻辑。
        2. **多 GPU 并行提速（先用单GPU验证流程，多GPU初步我已经做了，以后再优化）：** 这是工业级的难点。你要利用 OpenMM 的 CUDA 平台配置，让动辄百万原子的体系在多卡上飞起来。
        3. **解决 AM1-BCC 痛点：** 把你之前复盘总结的报错 fix 彻底自动化，让体系不仅能跑，而且不报错。
        4. 我现在Stage4是“OpenFF + GAFFTemplateGenerator 作为接口层，底层调用 AmberTools (antechamber/sqm) 计算 AM1-BCC”。

### 问题一

现在的是**`OpenFF + Modeller + PDBFixer` 是目前“开源+自动化 Python 流水线”里的准一线配置，但它绝不是最完美的。** 这套组合最大的“阿喀琉斯之踵”在于 **PDBFixer**。PDBFixer 足够轻量、方便自动化，但它在处理**蛋白质质子化状态（Protonation States）**时非常粗糙，而质子化状态直接决定了结合口袋里的氢键网络，甚至能导致 MD 模拟完全跑偏。

所以我是不是依然可以用 PDBFixer 补全缺失残基，但加氢这一步，**必须交由 PDB2PQR 根据 PROPKA 的结果来加**。

---

`sqm`（半经验量子化学）非常慢，而且极度脆弱。只要上游传过来的配体结构哪怕有一点点扭曲，`sqm` 计算 AM1-BCC 时就会无限不收敛或直接 Segmentation Fault。

然后我们**引入 `Espaloma` (基于图神经网络的 ML 电荷)**

如果 `sqm` 崩溃了，触发 `Espaloma` 兜底计算，**但是（最关键的一步）**：系统必须在这个分子的输出报告上打上一个红色的 `[WARNING: ML_CHARGE]` 标签。是不是这样做更好？

对于标准的水盒子（水 + NaCl 离子），OpenMM 的 Modeller.addSolvent 已经做得极其优秀且高度优化，不需要替换。除非你以后要做复杂的细胞膜跨膜蛋白模拟（需要换用 PACKMOL），我们现在先不用那些写上以后未来优化计划那里。

---

为什么先有个项目基础，我们是不是先**第一阶段（锁定现状）：** 坚决执行当前的 `OpenFF + PDBFixer + Modeller` 路线。在这个阶段，把你作为 Lead 的精力放在**严格执行交付契约**上，强制要求上游必须传标准的 3D SDF 和清理好的 PDB 单体。

**第二阶段（灰度升级）：** 等你的团队跑通了第一个完整的端到端（End-to-End）测试，并且在日志里收集到了那些让 `sqm` 崩溃的“坏分子”。拿着这些真实的报错案例，你再去引入 Espaloma 做 Fallback，或者用 PROPKA 替换 PDBFixer，这时的重构才是真正创造价值的。



