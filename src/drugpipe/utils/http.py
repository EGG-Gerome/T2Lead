"""HTTP client with automatic retry, exponential back-off, and polite sleep."""
# EN: Module overview and key intent for maintainers.
# 中文：模块总览与关键设计意图，便于后续维护。

# 带自动重试、指数退避与礼貌延时的 HTTP 客户端。

from __future__ import annotations

import random
import time
from typing import Any, Dict, Optional

import requests


# EN: HTTPClient core behavior and intent.
# 中文：HTTPClient 的核心行为与设计意图。
class HTTPClient:
    """Thin wrapper around ``requests.get`` with retry / back-off logic."""
    # 对 requests 的薄封装，支持重试与退避。

    # EN: __init__ core behavior and intent.
    # 中文：__init__ 的核心行为与设计意图。
    def __init__(
        self,
        timeout: int = 60,
        retries: int = 6,
        backoff_base: float = 1.6,
        polite_sleep: float = 0.05,
    ):
        self.timeout = timeout
        self.retries = retries
        self.backoff_base = backoff_base
        self.polite_sleep = polite_sleep
        self._session = requests.Session()
        self._session.headers.update({"Accept": "application/json"})

    # ------------------------------------------------------------------
    # EN: get_json core behavior and intent.
    # 中文：get_json 的核心行为与设计意图。
    def get_json(self, url: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """GET URL and return JSON; retry on 429/5xx. / GET 请求并返回 JSON，遇 429/5xx 重试。"""
        last_err: Optional[Exception] = None
        for attempt in range(self.retries):
            try:
                r = self._session.get(url, params=params, timeout=self.timeout)
                if r.status_code == 200:
                    time.sleep(self.polite_sleep)
                    return r.json()
                if r.status_code in (429, 500, 502, 503, 504):
                    last_err = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
                    self._backoff(attempt)
                    continue
                r.raise_for_status()
            except requests.RequestException as exc:
                last_err = exc
                self._backoff(attempt)
        raise RuntimeError(f"GET failed after {self.retries} retries: {url} — {last_err}")

    # EN: post_json core behavior and intent.
    # 中文：post_json 的核心行为与设计意图。
    def post_json(self, url: str, payload: Dict[str, Any]) -> Dict[str, Any]:
        """POST JSON payload and return JSON. / 发送 JSON 并返回响应 JSON。"""
        last_err: Optional[Exception] = None
        for attempt in range(self.retries):
            try:
                r = self._session.post(url, json=payload, timeout=self.timeout)
                if r.status_code == 200:
                    time.sleep(self.polite_sleep)
                    return r.json()
                if r.status_code in (429, 500, 502, 503, 504):
                    last_err = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
                    self._backoff(attempt)
                    continue
                r.raise_for_status()
            except requests.RequestException as exc:
                last_err = exc
                self._backoff(attempt)
        raise RuntimeError(f"POST failed after {self.retries} retries: {url} — {last_err}")

    # ------------------------------------------------------------------
    # EN: _backoff core behavior and intent.
    # 中文：_backoff 的核心行为与设计意图。
    def _backoff(self, attempt: int) -> None:
        """Exponential backoff before next retry. / 指数退避后再次重试。"""
        sleep_s = (self.backoff_base ** attempt) + random.random()
        time.sleep(sleep_s)
