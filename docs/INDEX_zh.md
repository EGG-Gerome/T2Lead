# 文档导航（中文）

## 当前建议结构

- `docs/guide/`：使用指南与配置说明
- `docs/reproducibility/`：复现实验步骤
- `docs/progress/`：阶段进展记录（周报/里程碑）
- `docs/postmortems/`：问题复盘与根因分析
- `docs/experiments/`：每次运行的固定归档（推荐主入口）

## 快速入口

- 归档规范：`docs/experiments/README_zh.md`
- 运行归档目录：`docs/experiments/runs/`
- 使用说明：`docs/guide/`
- Skill/Rule 使用：`docs/guide/skills_and_rules_usage_zh.md`
- 单/多 GPU 指南：`docs/guide/gpu_execution_zh.md`
- `variant_runs` 结构图解：`docs/guide/variant_runs_structure_zh.md`
- 团队执行规则：`docs/guide/team_execution_rules_zh.md`
- 模块改进路线图：`docs/guide/module_improvement_roadmap_zh.md`
- 复盘记录：`docs/postmortems/`

## 整理建议（先不破坏历史）

1. 先加导航与命名规范，不立即移动老文件。
2. 新运行统一落在 `docs/experiments/runs/...`，避免散落。
3. `postmortems` 仅保留“问题复盘”，不存普通运行报告。
4. 是否做历史归档（如 `docs/archive/`）按项目决策执行，不做自动迁移。
