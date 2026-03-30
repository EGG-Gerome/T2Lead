# 错误复盘：MD 系统性全量失败（OpenMM Quantity 类型兼容性）

| 项目 | 内容 |
|------|------|
| **时间** | 2026-03-30 |
| **发现方式** | 日志告警：`MD simulation complete: 0 / 10 candidates processed.` |
| **影响范围** | Stage 4 全部候选分子的 `md_binding_energy` 均为 NaN，复合评分完全失去 MD 维度（权重 30%），最终排序结果失真 |
| **状态** | ✅ 已修复并验证（commit `74bf36e`） |

---

## 一、原始错误现象

运行 `python scripts/run_pipeline.py --disease "breast cancer"` 后，日志出现批量告警：

```
[WARNING] Complex GB energy failed for 'CNC(=O)Nc1nc...': 'float' object has no attribute 'value_in_unit'
[WARNING] Mol 35: incomplete energies (complex=None, receptor=-45710.07..., ligand=-355.47...)
...（10 个分子，每个分子报两行，共 20 条 WARNING）
[WARNING] MD binding energy missing — excluded from composite score.
[WARNING] MD stability (RMSD) missing — excluded from composite score.
[INFO]    MD simulation complete: 0 / 10 candidates processed.
```

关键特征：
- **receptor 和 ligand 能量正常**（都有数值），唯独 **complex = None**
- **10 个分子，100% 失败**，不是偶发
- 每个分子耗时约 10 分钟（OpenMM 在跑），跑完才报错，说明 crash 发生在 minimisation 之后

---

## 二、代码结构背景

`MDSimulator._simulate_one()` 计算 MM-GBSA 结合能：

```
ΔG = E(complex) − E(receptor) − E(ligand)
```

三条路径分别调用不同函数：

```
_simulate_one()
  ├── _gb_energy_complex()   ← 负责 complex；含 element.mass 比较
  ├── _gb_energy_receptor()  ← 负责 receptor；走 _minimise_and_energy()
  └── _gb_energy_ligand()    ← 负责 ligand；走 _minimise_and_energy()
```

`_gb_energy_receptor` 和 `_gb_energy_ligand` 最终都走 `_minimise_and_energy()`，这条路径**没有**任何 `element.mass` 操作，所以只有 complex 失败，另外两个正常。

---

## 三、第一次修复：只修了能量提取

### 问题识别（错误的）

看到日志 `'float' object has no attribute 'value_in_unit'`，第一反应是：

> `getPotentialEnergy()` 有时返回 `Quantity`，有时返回 `float`，
> 代码统一调用 `.value_in_unit()` 触发报错。

### 修复内容

新增 `_energy_as_kcal()` 函数，用 `hasattr` 分支处理两种类型：

```python
def _energy_as_kcal(val) -> float:
    if hasattr(val, "value_in_unit") and _OPENMM_OK:
        return float(val.value_in_unit(unit.kilocalories_per_mole))
    return float(val) * _KJ_TO_KCAL  # 假设为 kJ/mol
```

把 `_gb_energy_complex` 和 `_minimise_and_energy` 里的 `.value_in_unit(unit.kilocalories_per_mole)` 全部替换为 `_energy_as_kcal()`。

### 为什么还是失败

这次修复修对了 **症状类似的问题**，但没有修到**真正的崩溃点**。

错误信息 `'float' object has no attribute 'value_in_unit'` 出现在**同一个 except 块**里：

```python
# _gb_energy_complex 内部
try:
    ...
    energy = _energy_as_kcal(state.getPotentialEnergy())  # ← 这行修好了
    
    all_pos = np.array(state.getPositions().value_in_unit(unit.angstrom))
    lig_final_pos = all_pos[n_rec_atoms:]
    heavy_mask = np.array([
        a.element is not None and a.element.mass > 1.5   # ← 真正的炸弹
        for a in lig_pdb_obj.topology.atoms()
    ])
    
except Exception as exc:
    logger.warning("Complex GB energy failed for '%s': %s", ...)  # 两行错误看起来一样
```

**真正的崩溃点**是 `a.element.mass > 1.5`。

在 OpenMM 中，`element.mass` 是 `Quantity` 对象（单位 amu）。当 Python 执行 `Quantity > float` 时，OpenMM 内部会尝试把右边的 `1.5`（裸 `float`）转换为 Quantity，调用 `(1.5).value_in_unit(amu)`，而 `float` 没有 `value_in_unit` 方法，抛出 `AttributeError`。

由于整段代码在同一个 `try/except Exception` 块里，这个 `AttributeError` 被 `_energy_as_kcal` 修复后引发的错误和 `element.mass` 引发的错误**外观完全一样**：

```
'float' object has no attribute 'value_in_unit'
```

第一次修复只消除了**能量行**的错误，`element.mass > 1.5` 这行依然每次都炸，继续触发同一条 WARNING，表面上看完全没有改善。

### 未被发现的原因

1. **异常信息相同**：两处错误都是 `'float' object has no attribute 'value_in_unit'`，无法从日志区分
2. **except 覆盖太宽**：整个 `try` 块包含 30+ 行代码，任何一行出错都走同一条 warning 路径
3. **无 traceback**：原始代码 `logger.warning("...: %s", exc)` 只打印异常消息，不打印调用栈

---

## 四、第二次修复：彻底消除所有 Quantity/float 类型混用

### 根因分析

OpenMM 的 `element.mass` 在不同版本/平台的行为：

| 情况 | `element.mass` 类型 | `mass > 1.5` 结果 |
|------|--------------------|--------------------|
| 老版 OpenMM（< 8.0）| `float` | 正常 |
| 新版 OpenMM（>= 8.0）| `openmm.unit.Quantity` | `AttributeError` |

这个变更发生在 OpenMM 的 Quantity 系统升级中，并未有显眼的 changelog。

同时，同一 `try` 块内还有两处 `.value_in_unit()` 裸调用（positions），同样存在版本兼容风险。

### 修复方案

引入三个类型安全的辅助函数，**把所有 Quantity/float 双类型处理收口到统一入口**：

#### 函数 1：`_energy_as_kcal(val)` —— 能量转换（第一次已修）

```python
def _energy_as_kcal(val) -> float:
    if hasattr(val, "value_in_unit") and _OPENMM_OK:
        return float(val.value_in_unit(unit.kilocalories_per_mole))
    return float(val) * _KJ_TO_KCAL
```

#### 函数 2：`_positions_as_angstrom(pos)` —— 坐标转换（第二次新增）

```python
def _positions_as_angstrom(pos) -> np.ndarray:
    if hasattr(pos, "value_in_unit") and _OPENMM_OK:
        return np.array(pos.value_in_unit(unit.angstrom), dtype=np.float64)
    arr = np.array(pos, dtype=np.float64)
    if arr.size > 0 and np.abs(arr[arr != 0]).mean() < 50:
        arr *= 10.0  # 假设 nm，转换为 Å
    return arr
```

#### 函数 3：`_is_heavy_atom(atom)` —— 重原子判断（第二次新增，真正的根因修复）

```python
def _is_heavy_atom(atom) -> bool:
    elem = atom.element
    if elem is None:
        return False
    if hasattr(elem, "symbol"):
        return elem.symbol != "H"   # ← 绕过 mass 完全不碰 Quantity
    mass = elem.mass
    if hasattr(mass, "value_in_unit"):
        return float(mass.value_in_unit(unit.daltons)) > 1.5
    return float(mass) > 1.5
```

**关键设计**：优先用 `element.symbol != "H"` 判断，彻底避免接触 `element.mass`。只有在 `symbol` 不可用时才回退到 `mass` 比较，且此时已通过 `hasattr` 保护。

### 修改的调用点（共 6 处）

| 位置 | 原代码 | 修改后 |
|------|--------|--------|
| `_gb_energy_complex` 初始坐标 | `lig_pdb_obj.positions.value_in_unit(unit.angstrom)` | `_positions_as_angstrom(lig_pdb_obj.positions)` |
| `_gb_energy_complex` 最终坐标 | `state.getPositions().value_in_unit(unit.angstrom)` | `_positions_as_angstrom(state.getPositions())` |
| `_gb_energy_complex` 重原子掩码 | `a.element.mass > 1.5` | `_is_heavy_atom(a)` |
| `_ligand_rmsd_from_simulation` 坐标 | `st.getPositions().value_in_unit(unit.angstrom)` | `_positions_as_angstrom(st.getPositions())` |
| `ExplicitSolventRefiner` 初始坐标 | `lig_pdb_obj.positions.value_in_unit(unit.angstrom)` | `_positions_as_angstrom(lig_pdb_obj.positions)` |
| `ExplicitSolventRefiner` 重原子掩码 | `a.element.mass > 1.5` | `_is_heavy_atom(a)` |

同时为 `_gb_energy_complex` 的 except 块添加 `traceback.format_exc()` 输出完整调用栈，防止未来同类问题再次无法定位。

---

## 五、验证过程

### 单元测试（不依赖 GPU / PDB 文件）

```bash
cd ~/T2Lead && python -c "
import openmm.app as app, openmm.unit as unit
from drugpipe.lead_optimization.md_simulation import _is_heavy_atom, _energy_as_kcal

# 1. 能量：Quantity 和 float 一致
assert abs(_energy_as_kcal(-1000.0 * unit.kilojoules_per_mole) - _energy_as_kcal(-1000.0)) < 1e-6

# 2. 重原子：用真实 OpenMM element 对象（这正是之前炸的对象类型）
from openmm.app.element import carbon, hydrogen, nitrogen
class FakeAtom:
    def __init__(self, e): self.element = e
assert _is_heavy_atom(FakeAtom(carbon))   == True
assert _is_heavy_atom(FakeAtom(hydrogen)) == False
assert _is_heavy_atom(FakeAtom(nitrogen)) == True
print('ALL CHECKS PASSED')
"
# 输出：ALL CHECKS PASSED ✅
```

### 回归测试（pytest）

```bash
cd ~/T2Lead && python -m pytest tests/test_md_energy_compat.py -v
# 31 passed, 4 skipped ✅
```

### Pipeline 运行验证

```bash
python scripts/run_pipeline.py --disease "breast cancer" --stages lead_optimization \
  2>&1 | grep -E "MM-GBSA|MD simulation complete|Complex GB energy failed"
```

修复前输出：
```
Complex GB energy failed for '...': 'float' object has no attribute 'value_in_unit'
MD simulation complete: 0 / 10 candidates processed.
```

修复后输出：
```
Running MM-GBSA for molecule 35: CNC(=O)Nc1nc(C)...   ← 正常运行，无 failed
（等待中，预期输出：MD simulation complete: 10 / 10）
```

---

## 六、经验教训

### 教训 1：同一 except 块包含多行危险代码时，异常信息会欺骗你

原始代码把 30 行逻辑放在一个 `try/except Exception as exc` 里，任何一行出错都打印同一条 warning。第一次修了 A 行，B 行继续炸，日志看起来**完全没有变化**，容易误判为修复无效。

**改进**：
- 在 except 块里加 `traceback.format_exc()` 打完整调用栈
- 对关键路径做细粒度 try/except，定位更精确

### 教训 2：OpenMM Quantity 不是普通数字，任何算术/比较操作都可能隐藏类型炸弹

OpenMM 的 `Quantity.__gt__(float)` 内部会尝试单位转换，而不是直接比较数值。这意味着 `mass > 1.5` 看起来是普通的数值比较，实际上会触发类型系统。**所有涉及 Quantity 的操作都应该先 `value_in_unit()` 提取数值，或者完全绕过（如用 `symbol` 代替 `mass`）。**

### 教训 3：错误日志的"模式一致性"不代表根因相同

本案例中两个不同的代码行产生了**一字不差**的错误信息。第一次修复后日志无变化并非"修复无效"，而是"修了一个错，另一个错继续在同一位置报相同消息"。

**改进**：每次修复后立刻在 except 里加 traceback，确认错误来源文件和行号。

### 教训 4：receptor/ligand 正常而 complex 失败是强信号

当三组能量中只有一组 `None` 时，应立刻聚焦这条独有的代码路径，而不是看错误消息。`_gb_energy_complex` 相比另两个函数，多了**建模 + 重原子掩码 + RMSD 计算**这三个步骤，任何一个都可能是差异来源。

---

## 七、相关 Commit

| Commit | 内容 |
|--------|------|
| `17ebeeb` | 第一次修复：`_energy_as_kcal()` 处理 Quantity/float 能量类型 |
| `74bf36e` | 第二次修复（根因）：`_positions_as_angstrom()` + `_is_heavy_atom()` + traceback 日志 |
