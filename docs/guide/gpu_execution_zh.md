# 单 GPU / 多 GPU 运行指南

## 目标

说明如何在 T2Lead 中切换单卡与多卡，覆盖：

- 配置文件方式（YAML）
- 环境变量方式（`DP_*`）
- CLI 方式（`run_mainpath.py`）

## 1) 配置文件方式（推荐）

文件：`configs/default_config.yaml`

关键项：

- `pipeline.device`: `auto | cuda | cpu | mps`
- `lead_optimization.md_simulation.parallel_workers`
- `lead_optimization.md_simulation.schedule_tasks_across_gpus`
- `lead_optimization.md_simulation.gpu_device_indices`

### 单 GPU 示例

```yaml
pipeline:
  device: cuda

lead_optimization:
  md_simulation:
    parallel_workers: 1
    schedule_tasks_across_gpus: false
    gpu_device_indices: [0]
```

### 多 GPU 示例（4 卡）

```yaml
pipeline:
  device: cuda

lead_optimization:
  md_simulation:
    parallel_workers: 4
    schedule_tasks_across_gpus: true
    gpu_device_indices: [0, 1, 2, 3]
```

## 2) 环境变量方式（`DP_*`）

```bash
# 单 GPU
export DP_PIPELINE__DEVICE=cuda
export DP_LEAD_OPTIMIZATION__MD_SIMULATION__PARALLEL_WORKERS=1
export DP_LEAD_OPTIMIZATION__MD_SIMULATION__SCHEDULE_TASKS_ACROSS_GPUS=false
export DP_LEAD_OPTIMIZATION__MD_SIMULATION__GPU_DEVICE_INDICES='[0]'

# 多 GPU
export DP_PIPELINE__DEVICE=cuda
export DP_LEAD_OPTIMIZATION__MD_SIMULATION__PARALLEL_WORKERS=4
export DP_LEAD_OPTIMIZATION__MD_SIMULATION__SCHEDULE_TASKS_ACROSS_GPUS=true
export DP_LEAD_OPTIMIZATION__MD_SIMULATION__GPU_DEVICE_INDICES='[0,1,2,3]'
```

> 注意：配置系统原生识别 `DP_*`。  
> `PARALLEL_WORKERS` 仅在 `scripts/run_mainpath.py` 中做了兼容读取，不是通用配置键。

## 3) CLI 方式（`run_mainpath.py`）

### 单 GPU

```bash
conda run -n t2lead python scripts/run_mainpath.py \
  --disease "breast cancer" \
  --vcf-path /path/to/sample.ann.vcf.gz \
  --sample-id patient_001 \
  --device cuda \
  --parallel-workers 1 \
  --no-schedule-tasks-across-gpus \
  --gpu-device-indices 0 \
  -v
```

### 多 GPU

```bash
conda run -n t2lead python scripts/run_mainpath.py \
  --disease "breast cancer" \
  --vcf-path /path/to/sample.ann.vcf.gz \
  --sample-id patient_001 \
  --device cuda \
  --parallel-workers 4 \
  --schedule-tasks-across-gpus \
  --gpu-device-indices 0,1,2,3 \
  -v
```

## 4) 建议

- 单卡机器：`parallel_workers=1` 更稳。
- 多卡机器：`parallel_workers` 不建议盲目高于可用 GPU 数。
- 长任务请配合 `tmux/screen/nohup`，避免 SSH 断连中断任务。
