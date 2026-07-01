from __future__ import annotations

import hashlib
import os
from collections.abc import Mapping
from pathlib import Path
import re
import sys


CACHE_ENV_VAR = "VSBTOOLS_CACHE_DIR"
_WINDOWS_RESERVED_NAMES = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{i}" for i in range(1, 10)),
    *(f"LPT{i}" for i in range(1, 10)),
}
_INVALID_FILENAME_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]+')


def user_cache_dir(
        *,
        environ: Mapping[str, str] | None = None,
        platform: str | None = None,
        home: str | Path | None = None,
) -> Path:
    """Return the per-user vsbtools cache directory without creating it."""
    environ = os.environ if environ is None else environ
    platform = sys.platform if platform is None else platform
    home = Path.home() if home is None else Path(home)

    configured = environ.get(CACHE_ENV_VAR)
    if configured:
        return Path(configured).expanduser()

    if platform == "win32":
        local_app_data = environ.get("LOCALAPPDATA")
        base = Path(local_app_data) if local_app_data else home / "AppData" / "Local"
        return base / "vsbtools" / "Cache"

    if platform == "darwin":
        return home / "Library" / "Caches" / "vsbtools"

    xdg_cache_home = environ.get("XDG_CACHE_HOME")
    base = Path(xdg_cache_home).expanduser() if xdg_cache_home else home / ".cache"
    return base / "vsbtools"


def database_cache_dir(**kwargs) -> Path:
    """Return the database cache directory without creating it."""
    return user_cache_dir(**kwargs) / "DB_caches"


def safe_cache_component(value: object, *, max_length: int = 80) -> str:
    """Return a bounded filename component valid on Windows and POSIX."""
    original = str(value)
    component = _INVALID_FILENAME_CHARS.sub("_", original)
    component = re.sub(r"\s+", "_", component).strip(" ._")
    if not component:
        component = "cache"
    if component.upper() in _WINDOWS_RESERVED_NAMES:
        component = f"_{component}"

    if len(component) > max_length:
        digest = hashlib.sha256(original.encode("utf-8")).hexdigest()[:8]
        component = f"{component[:max_length - len(digest) - 1]}_{digest}"
    return component
