#!/usr/bin/env python3
"""Convenience script to run the full DrugPipe pipeline."""
# Module overview and key intent for maintainers.
# 模块总览与关键设计意图，便于后续维护。

# 一键运行完整 DrugPipe 流水线的便捷脚本。

import os
import sys
from pathlib import Path

conda_lib = Path(sys.executable).resolve().parents[1] / "lib"
if conda_lib.exists():
    os.environ["LD_LIBRARY_PATH"] = str(conda_lib) + ":" + os.environ.get("LD_LIBRARY_PATH", "")

# Ensure src/ is on the import path when running from project root
# 从项目根目录运行时确保 src/ 在导入路径中
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from drugpipe.pipeline import main

if __name__ == "__main__":
    main()
