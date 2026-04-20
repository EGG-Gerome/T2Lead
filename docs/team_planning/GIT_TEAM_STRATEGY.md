# Git 团队协作策略

## 结论（先说）

**推荐：GitHub Organization + 多仓库（每模块一个仓库 + 一个团队管理仓库）。**

---

## 1) 单仓 vs 多仓的选择

| 维度 | 单仓库（Monorepo） | 多仓库（Multi-repo） | 本项目建议 |
|------|-------------------|---------------------|-----------|
| 权限隔离 | 难（GitHub 不支持目录级权限） | 易（按仓库分配） | **多仓胜** |
| 跨模块耦合开发 | 方便 | 需要发版对接 | — |
| 依赖版本管理 | 复杂 | 清晰（每模块独立 pin） | **多仓胜** |
| CI 成本 | 高（动一点全跑） | 按仓库独立跑 | **多仓胜** |
| 学生/新人上手 | 一进来看到所有代码 | 只看自己的 | **多仓胜** |
| 模块联调 | 直接 import | 用文件交接契约 | 本项目首选文件交接，多仓胜 |

**关键判断**：本项目 6 个模块技术栈差异大（Python + Scanpy / OpenMM / PyTorch / RNA 工具 / 蛋白设计），耦合度低（靠文件交接），人员分工清晰。**这是经典的 Multi-repo 场景**。

---

## 2) 推荐的 Organization + 仓库结构

### 2.1 创建 GitHub Organization
- 建议名字：`T2Lead-Team` 或 `<你定的代号>`
- 你是 **Owner**，拥有所有权限

### 2.2 仓库列表

| 仓库 | 负责人 | 写权限 | 其他成员权限 |
|------|--------|-------|-------------|
| `team-planning` | Lead（你） | Lead 专属 | **Read-only**，所有人只能看 |
| `biodriver` | BioDriver 负责人 | 该负责人 + Lead | **Read-only** + PR |
| `drugreflector-scrna` | 姜可盈 | 她 + Lead | **Read-only** + PR |
| `rnai-therapy` | RNAi 负责人 | 该负责人 + Lead | **Read-only** + PR |
| `precision-delivery-adc` | ADC 负责人 | 该负责人 + Lead | **Read-only** + PR |
| `immunogen-mrna` | ImmunoGen 负责人 | 该负责人 + Lead | **Read-only** + PR |
| `simulation-hub` | Lead（你） | Lead 专属 | **Read-only** + PR |
| `shared-contracts`（可选） | Lead | Lead 专属 | Read-only | 存公共 JSON schema / 契约定义 |

### 2.3 权限实现方式（GitHub Teams）

- 创建 **teams**：
  - `leads`（含你）→ 所有仓库 Admin
  - `biodriver-team` / `drugreflector-team` / ... 每人加入自己模块的 team
  - `all-members` → 所有人，默认 Read
- 每个模块仓库：对应 team 给 **Maintain / Write**，`all-members` 给 **Read**，`leads` 给 **Admin**
- `team-planning` 和 `simulation-hub`：只有 `leads` 有 Write，其他给 Read

> **结果**：每人可以改自己的模块；其他模块只能 Read + 提 PR；只有你能改 team-planning 和 simulation-hub。

---

## 3) 分支与 PR 规则（建议最小集）

- 主分支 `main` 保护：
  - 至少 1 个 review 才能 merge
  - CI 必须通过
  - 禁止直接 push（走 PR）
- Feature 分支命名：`feat/<short-desc>`、`fix/<...>`
- 每 PR 至少包含：
  - 一句话总结（Summary）
  - 测试计划（Test plan）
- Lead 作为 owner 可以跨仓 review 和 merge

---

## 4) 跨模块协作（文件交接，不直接 import）

- 模块间不强行互相 import 代码。
- 交付走"文件 + 约定 schema"（见 `ARCHITECTURE_AND_DATAFLOW.md`）。
- 如果需要共享契约定义（JSON schema / CSV 列名），放到 `shared-contracts` 仓库，所有人 Read。

---

## 5) 共享存储（交付产物）

代码仓库不放大文件。交付产物用共享存储：

| 层 | 位置 | 用途 |
|---|------|-----|
| 代码 | GitHub Organization 各仓库 | 工程代码 + 测试 |
| 交付产物 | autodl-fs / S3 / 公司 NAS | `deliveries/<run_id>/...`（矩阵、结构、轨迹等） |
| 文档 | GitHub（各仓库 `docs/`） + `team-planning` | 设计、接口、报告 |

---

## 6) 联通性测试（第一阶段必做）

**目标**：BioDriver → 4 下游 → Simulation Hub 端到端跑通一次。

**建议做法**：
1. 在 `team-planning` 仓库（或单独一个 `e2e-test` 仓库）写 1 个 orchestration 脚本。
2. 脚本按序调用各模块的 `smoke_test.sh`。
3. 每次联通测试产出一份 `e2e_report_<timestamp>.md` 归档。

---

## 7) 常见问题

**Q: 一个人能不能参与多个模块？**  
A: 可以，加入多个 team 即可。权限累加。

**Q: 新人入职怎么办？**  
A: Lead 把他加入对应 team，给对应仓库 Write。其他仓库默认 Read。

**Q: 如果我不想用 GitHub Organization，个人账号下行不行？**  
A: 也可以，但权限管理会更粗糙（GitHub 个人账号下的协作者权限不能像 Org 那样按 team 管理）。强烈建议开 Organization，免费的就够用。

**Q: 如果第一阶段还没联通就发现接口要改怎么办？**  
A: 在 `team-planning` 仓库改 `ARCHITECTURE_AND_DATAFLOW.md` + 在 `shared-contracts` 改 schema；发 PR 通知所有负责人。

---

## 8) 对 Lead 的操作清单（第一周可以做的事）

1. 创建 GitHub Organization。
2. 创建 7-8 个仓库 + 对应 teams。
3. 每个仓库放初始 `README.md`（可以直接把 `docs/team_planning/modules/<module>.md` 作为初始内容）。
4. `main` 分支保护规则打开。
5. 把 `docs/team_planning/` 整目录推到 `team-planning` 仓库（这是你的"总指挥部"）。
6. 拉所有成员入职 + 分配 team。
