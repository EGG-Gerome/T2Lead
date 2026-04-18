"""HTTP client with automatic retry, exponential back-off, and polite sleep."""
# HTTP client with automatic retry, exponential back-off, and polite sleep.
# 说明模块职责、上下游关系与维护注意事项。

# 带自动重试、指数退避与礼貌延时的 HTTP 客户端。

from __future__ import annotations

import random
import time
from typing import Any, Dict, Optional

import requests


class HTTPClient:
    """Thin wrapper around ``requests.get`` with retry / back-off logic."""
    # 对 requests 的薄封装，支持重试与退避。

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

    def post_json(self, url: str, payload: Dict[str, Any]) -> Dict[str, Any]:
        """POST JSON payload and return JSON. / 发送 JSON 并返回响应 JSON。"""
        last_err: Optional[Exception] = None
        for attempt in range(self.retries):
            try:
                r = self._session.post(url, json=payload, timeout=self.timeout)
                if r.status_code == 204:
                    time.sleep(self.polite_sleep)
                    return {}
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

    def get_text(
        self,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        *,
        headers: Optional[Dict[str, str]] = None,
    ) -> str:
        """GET URL and return decoded text (non-JSON bodies: PDB, plain text)."""
        last_err: Optional[Exception] = None
        hdr = dict(headers or {})
        if "Accept" not in hdr:
            hdr["Accept"] = "*/*"
        for attempt in range(self.retries):
            try:
                r = self._session.get(
                    url, params=params, timeout=self.timeout, headers=hdr,
                )
                if r.status_code == 200:
                    time.sleep(self.polite_sleep)
                    return r.text
                if r.status_code in (429, 500, 502, 503, 504):
                    last_err = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
                    self._backoff(attempt)
                    continue
                r.raise_for_status()
            except requests.RequestException as exc:
                last_err = exc
                self._backoff(attempt)
        raise RuntimeError(f"GET failed after {self.retries} retries: {url} — {last_err}")

    def post_text(
        self,
        url: str,
        data: Any,
        *,
        headers: Optional[Dict[str, str]] = None,
    ) -> str:
        """POST raw body and return text (e.g. ESMFold PDB)."""
        last_err: Optional[Exception] = None
        hdr = dict(headers or {})
        if "Accept" not in hdr:
            hdr["Accept"] = "*/*"
        for attempt in range(self.retries):
            try:
                r = self._session.post(
                    url, data=data, timeout=self.timeout, headers=hdr,
                )
                if r.status_code == 200:
                    time.sleep(self.polite_sleep)
                    return r.text
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
    def _backoff(self, attempt: int) -> None:
        """Exponential backoff before next retry. / 指数退避后再次重试。"""
        sleep_s = (self.backoff_base ** attempt) + random.random()
        time.sleep(sleep_s)
