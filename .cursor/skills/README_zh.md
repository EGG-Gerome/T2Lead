# `.cursor/skills` 技能架构指引（中文）

本目录用于存放团队可共享的项目级 Skills。  
每个 Skill 必须是一个独立目录，包含 `SKILL.md`。

## 目录约定

```text
.cursor/skills/
  <skill-name>/
    SKILL.md
    (可选) reference.md / examples.md / scripts/
```

## 使用原则

- Skill 属于“建议型流程模板”，团队成员可选用、可改造。
- 硬性制度写在 `docs/guide/team_execution_rules_zh.md`，不要全部塞进 Skill。

## 新增 Skill 的强制要求

1. 新增 Skill 时，**必须更新本文件**（新增一条目录说明）。
2. `SKILL.md` frontmatter 必须包含：
   - `name`
   - `description`（写清“做什么 + 何时触发”）。
3. 名称使用小写短横线，避免歧义。

## 当前 Skill 列表

- `experiment-reporting`：实验运行报告与归档产出
- `t2lead-run-guardrails`：运行守则（环境/测试/证据/长任务）
- `code-task-checklist`：小步实现 + 即时测试 + git 范围检查
- `architecture-note`：架构决策与取舍记录
- `bugfix-workflow`：缺陷定位、修复、回归验证
- `testing-workflow`：分层测试流程（单测/冒烟/运行检查）
- `settings-change-zh`：配置修改流程与验证

> `weekly-reporting-zh` 已改为个人私有，不再作为仓库共享 Skill。

## 与 Rule 的关系

- Skill：建议型、场景触发。
- Rule（`.cursor/rules/*.mdc`）：行为约束，可设置 `alwaysApply`。
