# Skills 与 Rules 使用说明（团队版）

## 放在哪里

- 项目级 Skill：`.cursor/skills/<skill-name>/SKILL.md`
- 项目级 Rule：`.cursor/rules/<rule-name>.mdc`
- Skill 架构总览：`.cursor/skills/README_zh.md`

这些路径在仓库内，提交后团队成员 `git pull` 就能获得。
当前仓库 `.gitignore` 未忽略 `.cursor/`，因此可正常纳入版本控制。

## 二者区别

- **Skill（建议型）**：按任务场景触发，可选择使用，也可自行修改。
- **Rule（约束型）**：用于统一行为规范。若 `alwaysApply: true` 则默认全局生效。

## 本项目建议策略

- 团队规范写在 `docs/guide/team_execution_rules_zh.md`（硬性制度）。
- Cursor Skill/Rule 作为执行辅助（建议或半自动约束）。
- 若希望成员“可选使用”，将 Rule 设为 `alwaysApply: false`。
- 当前 `t2lead-execution-safety.mdc` 已设置为 `alwaysApply: false`（建议型）。

## 个人强制与团队共享如何兼容

- 团队仓库里的 Rule 保持建议型（`alwaysApply: false`）。
- 你个人若想强制生效，不要改仓库文件反复提交；改用个人本地规则（`~/.cursor/rules/`）设置 `alwaysApply: true`。
- 这样可以避免每次 `git push` 都产生无意义配置差异。

## 新增一个 Skill 的最小步骤

1. 新建目录：`.cursor/skills/<name>/`
2. 新建 `SKILL.md`（包含 frontmatter 的 `name` 和 `description`）
3. 在描述里写清触发场景关键词
4. 提交到仓库并在团队文档登记

## 触发建议

- 在 `description` 中写“做什么 + 何时用（触发词）”。
- 关键词尽量覆盖团队常用表述（如“周报”“复盘”“测试”“修复”）。

## 当前可用团队 Skill

- `code-task-checklist`
- `architecture-note`
- `bugfix-workflow`
- `testing-workflow`
- `settings-change-zh`
- `experiment-reporting`
- `t2lead-run-guardrails`

说明：`weekly-reporting-zh` 为个人私有 Skill（不随仓库共享）。
