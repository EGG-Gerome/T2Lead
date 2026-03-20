"""ML / DL models for pIC50 regression: RandomForest + optional Torch MLP."""
# EN: Module overview and key intent for maintainers.
# 中文：模块总览与关键设计意图，便于后续维护。

# pIC50 回归的 ML/DL 模型：随机森林 + 可选 Torch MLP。

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import joblib
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, train_test_split

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
    # EN: _MLPRegressor core behavior and intent.
    # 中文：_MLPRegressor 的核心行为与设计意图。
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

# EN: ModelTrainer core behavior and intent.
# 中文：ModelTrainer 的核心行为与设计意图。
class ModelTrainer:
    """Train RF (always) + Torch MLP (if available), expose predict/evaluate."""
    # 训练 RF（必选）+ Torch MLP（可选），提供 predict/evaluate 接口。

    # EN: __init__ core behavior and intent.
    # 中文：__init__ 的核心行为与设计意图。
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

        # Device: auto | cuda | mps | cpu
        self._device_cfg = cfg.get("pipeline", {}).get("device", "auto")

        self.rf: Optional[RandomForestRegressor] = None
        self.mlp: Optional[Any] = None
        self._device: Optional[Any] = None
        self._cache_dir: Optional[Path] = None

    # ------------------------------------------------------------------
    # EN: train core behavior and intent.
    # 中文：train 的核心行为与设计意图。
    def train(
        self, X: np.ndarray, y: np.ndarray, cache_dir: Optional[Path] = None,
    ) -> Dict[str, float]:
        """Train both models with k-fold CV evaluation, return metrics.

        If *cache_dir* is given and a matching model cache exists (same
        data hash + hyperparameters), loads it instead of retraining.
        """
        self._cache_dir = cache_dir
        if cache_dir is not None:
            loaded_metrics = self._try_load(X, y)
            if loaded_metrics is not None:
                return loaded_metrics

        X_tr, X_te, y_tr, y_te = train_test_split(
            X, y, test_size=self.test_size, random_state=self.seed,
        )

        cv_metrics = self._cross_validate(X_tr, y_tr)

        logger.info("Training final RandomForest (n_estimators=%d) ...", self.rf_n_est)
        self.rf = RandomForestRegressor(
            n_estimators=self.rf_n_est,
            max_depth=None,
            n_jobs=-1,
            random_state=self.seed,
            min_samples_leaf=self.rf_min_leaf,
        )
        self.rf.fit(X_tr, y_tr)
        pred_te_rf = self.rf.predict(X_te)
        rmse_rf = _rmse(y_te, pred_te_rf)
        r2_rf = float(r2_score(y_te, pred_te_rf))
        logger.info("  RF held-out RMSE = %.4f, R² = %.4f", rmse_rf, r2_rf)

        metrics: Dict[str, float] = {
            "rmse_rf": rmse_rf,
            "r2_rf": r2_rf,
            **cv_metrics,
        }

        if _TORCH_OK:
            logger.info("Training final Torch MLP (epochs=%d) ...", self.mlp_epochs)
            self.mlp, rmse_mlp = self._train_mlp(X_tr, y_tr, X_te, y_te)
            pred_te_mlp = self._predict_mlp(X_te)
            r2_mlp = float(r2_score(y_te, pred_te_mlp))
            logger.info("  MLP held-out RMSE = %.4f, R² = %.4f", rmse_mlp, r2_mlp)
            metrics["rmse_mlp"] = rmse_mlp
            metrics["r2_mlp"] = r2_mlp
        else:
            logger.info("Torch not available — skipping MLP.")

        if cache_dir is not None:
            self._save(X, y, metrics)

        return metrics

    # ------------------------------------------------------------------
    # Model persistence
    # ------------------------------------------------------------------
    # EN: _data_hash core behavior and intent.
    # 中文：_data_hash 的核心行为与设计意图。
    def _data_hash(self, X: np.ndarray, y: np.ndarray) -> str:
        h = hashlib.sha256()
        h.update(f"n={len(y)},d={X.shape[1]},seed={self.seed}".encode())
        h.update(y[:64].tobytes())
        return h.hexdigest()[:12]

    # EN: _save core behavior and intent.
    # 中文：_save 的核心行为与设计意图。
    def _save(self, X: np.ndarray, y: np.ndarray, metrics: Dict[str, float]) -> None:
        d = self._cache_dir
        d.mkdir(parents=True, exist_ok=True)
        tag = self._data_hash(X, y)
        joblib.dump(self.rf, d / f"rf_{tag}.joblib")
        joblib.dump(metrics, d / f"metrics_{tag}.joblib")
        if self.mlp is not None and _TORCH_OK:
            import torch
            torch.save(self.mlp.state_dict(), d / f"mlp_{tag}.pt")
            joblib.dump(X.shape[1], d / f"mlp_dim_{tag}.joblib")
        logger.info("Saved model cache: %s (tag=%s)", d, tag)

    # EN: _try_load core behavior and intent.
    # 中文：_try_load 的核心行为与设计意图。
    def _try_load(self, X: np.ndarray, y: np.ndarray) -> Optional[Dict[str, float]]:
        d = self._cache_dir
        if d is None or not d.exists():
            return None
        tag = self._data_hash(X, y)
        rf_path = d / f"rf_{tag}.joblib"
        metrics_path = d / f"metrics_{tag}.joblib"
        if not rf_path.exists() or not metrics_path.exists():
            return None
        self.rf = joblib.load(rf_path)
        metrics = joblib.load(metrics_path)
        mlp_path = d / f"mlp_{tag}.pt"
        dim_path = d / f"mlp_dim_{tag}.joblib"
        if _TORCH_OK and mlp_path.exists() and dim_path.exists():
            import torch
            in_dim = joblib.load(dim_path)
            self._device = _resolve_torch_device(self._device_cfg)
            self.mlp = _MLPRegressor(in_dim).to(self._device)
            self.mlp.load_state_dict(torch.load(mlp_path, map_location=self._device, weights_only=True))
            self.mlp.eval()
        logger.info("Loaded cached model: %s (tag=%s)", d, tag)
        logger.info("  Cached metrics: %s", metrics)
        return metrics

    # ------------------------------------------------------------------
    # EN: _cross_validate core behavior and intent.
    # 中文：_cross_validate 的核心行为与设计意图。
    def _cross_validate(
        self, X: np.ndarray, y: np.ndarray, n_folds: int = 5,
    ) -> Dict[str, float]:
        """Run k-fold CV on the training set and report per-fold RMSE/R²."""
        logger.info("Running %d-fold cross-validation ...", n_folds)
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=self.seed)
        fold_rmse: list[float] = []
        fold_r2: list[float] = []

        for fold, (train_idx, val_idx) in enumerate(kf.split(X), 1):
            X_f_tr, X_f_val = X[train_idx], X[val_idx]
            y_f_tr, y_f_val = y[train_idx], y[val_idx]

            rf = RandomForestRegressor(
                n_estimators=self.rf_n_est,
                max_depth=None,
                n_jobs=-1,
                random_state=self.seed,
                min_samples_leaf=self.rf_min_leaf,
            )
            rf.fit(X_f_tr, y_f_tr)
            pred = rf.predict(X_f_val)
            rmse = _rmse(y_f_val, pred)
            r2 = float(r2_score(y_f_val, pred))
            fold_rmse.append(rmse)
            fold_r2.append(r2)
            logger.info("  Fold %d/%d: RMSE=%.4f  R²=%.4f", fold, n_folds, rmse, r2)

        mean_rmse = float(np.mean(fold_rmse))
        std_rmse = float(np.std(fold_rmse))
        mean_r2 = float(np.mean(fold_r2))
        std_r2 = float(np.std(fold_r2))
        logger.info(
            "  CV summary: RMSE=%.4f ± %.4f, R²=%.4f ± %.4f",
            mean_rmse, std_rmse, mean_r2, std_r2,
        )
        return {
            "cv_rmse_mean": mean_rmse,
            "cv_rmse_std": std_rmse,
            "cv_r2_mean": mean_r2,
            "cv_r2_std": std_r2,
        }

    # ------------------------------------------------------------------
    # EN: predict core behavior and intent.
    # 中文：predict 的核心行为与设计意图。
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Ensemble prediction: average of RF and MLP (if available)."""
        # 集成预测：RF 与 MLP（若有）的平均。
        pred_rf = self.rf.predict(X)
        if self.mlp is not None:
            pred_mlp = self._predict_mlp(X)
            return (pred_rf + pred_mlp) / 2.0
        return pred_rf

    # EN: predict_rf core behavior and intent.
    # 中文：predict_rf 的核心行为与设计意图。
    def predict_rf(self, X: np.ndarray) -> np.ndarray:
        return self.rf.predict(X)

    # EN: predict_mlp core behavior and intent.
    # 中文：predict_mlp 的核心行为与设计意图。
    def predict_mlp(self, X: np.ndarray) -> Optional[np.ndarray]:
        if self.mlp is None:
            return None
        return self._predict_mlp(X)

    # ------------------------------------------------------------------
    # Internal MLP helpers / 内部 MLP 辅助
    # ------------------------------------------------------------------
    # EN: _train_mlp core behavior and intent.
    # 中文：_train_mlp 的核心行为与设计意图。
    def _train_mlp(
        self, X_tr: np.ndarray, y_tr: np.ndarray, X_te: np.ndarray, y_te: np.ndarray,
    ) -> Tuple[Any, float]:
        device = _resolve_torch_device(self._device_cfg)
        self._device = device
        logger.info("  MLP device: %s", device)

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

    # EN: _predict_mlp core behavior and intent.
    # 中文：_predict_mlp 的核心行为与设计意图。
    def _predict_mlp(self, X: np.ndarray) -> np.ndarray:
        return self._predict_torch(self.mlp, self._device, X)

    @staticmethod
    # EN: _predict_torch core behavior and intent.
    # 中文：_predict_torch 的核心行为与设计意图。
    def _predict_torch(model: Any, device: Any, X: np.ndarray, batch_size: int = 4096) -> np.ndarray:
        model.eval()
        preds = []
        with torch.no_grad():
            for i in range(0, len(X), batch_size):
                t = torch.tensor(X[i : i + batch_size], dtype=torch.float32, device=device)
                preds.append(model(t).detach().cpu().numpy())
        return np.concatenate(preds)


# EN: _resolve_torch_device core behavior and intent.
# 中文：_resolve_torch_device 的核心行为与设计意图。
def _resolve_torch_device(device_cfg: str) -> "torch.device":
    """Resolve config device string to a torch.device.

    Supported values: ``"auto"``, ``"cuda"``, ``"mps"``, ``"cpu"``.
    """
    d = device_cfg.lower().strip()
    if d == "cuda":
        return torch.device("cuda")
    if d == "mps":
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps")
        logger.warning("MPS requested but not available, falling back to CPU.")
        return torch.device("cpu")
    if d == "cpu":
        return torch.device("cpu")
    # auto
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


# EN: _rmse core behavior and intent.
# 中文：_rmse 的核心行为与设计意图。
def _rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))
