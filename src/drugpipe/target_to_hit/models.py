"""ML / DL models for pIC50 regression: RandomForest + optional Torch MLP."""
# pIC50 回归的 ML/DL 模型：随机森林 + 可选 Torch MLP。

from __future__ import annotations

import logging
from typing import Any, Dict, Optional, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

logger = logging.getLogger(__name__)

# Torch is optional / Torch 为可选依赖
try:
    import torch
    import torch.nn as nn
    import torch.optim as optim

    _TORCH_OK = True
except ImportError:
    _TORCH_OK = False


# ======================================================================
# Torch MLP
# Torch 多层感知机
# ======================================================================

if _TORCH_OK:
    class _MLPRegressor(nn.Module):
        def __init__(self, in_dim: int):
            super().__init__()
            self.net = nn.Sequential(
                nn.Linear(in_dim, 1024),
                nn.ReLU(),
                nn.Dropout(0.25),
                nn.Linear(1024, 512),
                nn.ReLU(),
                nn.Dropout(0.25),
                nn.Linear(512, 256),
                nn.ReLU(),
                nn.Dropout(0.15),
                nn.Linear(256, 1),
            )

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            return self.net(x).squeeze(-1)


# ======================================================================
# Unified trainer
# 统一训练器
# ======================================================================

class ModelTrainer:
    """Train RF (always) + Torch MLP (if available), expose predict/evaluate."""
    # 训练 RF（必选）+ Torch MLP（可选），提供 predict/evaluate 接口。

    def __init__(self, cfg: Dict[str, Any]):
        mcfg = cfg.get("target_to_hit", {}).get("model", {})
        self.seed = int(cfg.get("pipeline", {}).get("seed", 42))
        self.test_size = float(mcfg.get("test_size", 0.2))

        # RF params / 随机森林参数
        self.rf_n_est = int(mcfg.get("rf_n_estimators", 600))
        self.rf_min_leaf = int(mcfg.get("rf_min_samples_leaf", 2))

        # MLP params / MLP 参数
        self.mlp_epochs = int(mcfg.get("mlp_epochs", 25))
        self.mlp_bs = int(mcfg.get("mlp_batch_size", 256))
        self.mlp_lr = float(mcfg.get("mlp_lr", 2e-3))
        self.mlp_wd = float(mcfg.get("mlp_weight_decay", 1e-4))

        self.rf: Optional[RandomForestRegressor] = None
        self.mlp: Optional[Any] = None  # torch model or None / Torch 模型或 None
        self._device: Optional[Any] = None

    # ------------------------------------------------------------------
    def train(self, X: np.ndarray, y: np.ndarray) -> Dict[str, float]:
        """Train both models, return test-set RMSE metrics."""
        # 训练两个模型，返回测试集 RMSE 指标。
        X_tr, X_te, y_tr, y_te = train_test_split(
            X, y, test_size=self.test_size, random_state=self.seed,
        )

        # Random Forest
        logger.info("Training RandomForest (n_estimators=%d) ...", self.rf_n_est)
        self.rf = RandomForestRegressor(
            n_estimators=self.rf_n_est,
            max_depth=None,
            n_jobs=-1,
            random_state=self.seed,
            min_samples_leaf=self.rf_min_leaf,
        )
        self.rf.fit(X_tr, y_tr)
        rmse_rf = _rmse(y_te, self.rf.predict(X_te))
        logger.info("  RF test RMSE = %.4f", rmse_rf)

        metrics = {"rmse_rf": rmse_rf}

        # Torch MLP
        if _TORCH_OK:
            logger.info("Training Torch MLP (epochs=%d) ...", self.mlp_epochs)
            self.mlp, rmse_mlp = self._train_mlp(X_tr, y_tr, X_te, y_te)
            logger.info("  MLP test RMSE = %.4f", rmse_mlp)
            metrics["rmse_mlp"] = rmse_mlp
        else:
            logger.info("Torch not available — skipping MLP.")

        return metrics

    # ------------------------------------------------------------------
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Ensemble prediction: average of RF and MLP (if available)."""
        # 集成预测：RF 与 MLP（若有）的平均。
        pred_rf = self.rf.predict(X)
        if self.mlp is not None:
            pred_mlp = self._predict_mlp(X)
            return (pred_rf + pred_mlp) / 2.0
        return pred_rf

    def predict_rf(self, X: np.ndarray) -> np.ndarray:
        return self.rf.predict(X)

    def predict_mlp(self, X: np.ndarray) -> Optional[np.ndarray]:
        if self.mlp is None:
            return None
        return self._predict_mlp(X)

    # ------------------------------------------------------------------
    # Internal MLP helpers / 内部 MLP 辅助
    # ------------------------------------------------------------------
    def _train_mlp(
        self, X_tr: np.ndarray, y_tr: np.ndarray, X_te: np.ndarray, y_te: np.ndarray,
    ) -> Tuple[Any, float]:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._device = device

        ds = torch.utils.data.TensorDataset(
            torch.tensor(X_tr, dtype=torch.float32),
            torch.tensor(y_tr, dtype=torch.float32),
        )
        dl = torch.utils.data.DataLoader(ds, batch_size=self.mlp_bs, shuffle=True)

        model = _MLPRegressor(X_tr.shape[1]).to(device)
        opt = optim.AdamW(model.parameters(), lr=self.mlp_lr, weight_decay=self.mlp_wd)
        loss_fn = nn.MSELoss()

        model.train()
        for ep in range(1, self.mlp_epochs + 1):
            losses = []
            for xb, yb in dl:
                xb, yb = xb.to(device), yb.to(device)
                opt.zero_grad(set_to_none=True)
                loss = loss_fn(model(xb), yb)
                loss.backward()
                opt.step()
                losses.append(loss.item())
            if ep % 5 == 0 or ep == 1 or ep == self.mlp_epochs:
                logger.info("    epoch %d/%d  loss=%.4f", ep, self.mlp_epochs, float(np.mean(losses)))

        pred = self._predict_torch(model, device, X_te)
        return model, _rmse(y_te, pred)

    def _predict_mlp(self, X: np.ndarray) -> np.ndarray:
        return self._predict_torch(self.mlp, self._device, X)

    @staticmethod
    def _predict_torch(model: Any, device: Any, X: np.ndarray, batch_size: int = 4096) -> np.ndarray:
        model.eval()
        preds = []
        with torch.no_grad():
            for i in range(0, len(X), batch_size):
                t = torch.tensor(X[i : i + batch_size], dtype=torch.float32, device=device)
                preds.append(model(t).detach().cpu().numpy())
        return np.concatenate(preds)


def _rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))
