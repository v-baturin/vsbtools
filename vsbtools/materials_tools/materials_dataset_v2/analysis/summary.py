from __future__ import annotations

from typing import Any, Callable, Collection, Dict, Iterable, List, Mapping, Sequence

#--optional pandas import ----------------------------------------------
try:
    import pandas as _pd
except ModuleNotFoundError as _err:          # fail late, not at import time
    _pd = None
    _PANDAS_ERR = _err


# -- helpers -------------------------------------------------------------
def _value_from_native(entry, col: str) -> Any:
    """Return a native attribute or metadata.<key> value (dot syntax)."""
    if col.startswith("metadata."):
        k = col.split(".", 1)[1]
        return entry.metadata.get(k)
    return getattr(entry, col, None)


def _rows_from_dataset(
    ds, native_cols: Sequence[str], callables: Mapping[str, Callable]
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []

    for e in ds:
        row: Dict[str, Any] = {col: _value_from_native(e, col) for col in native_cols}

        for name, fn in callables.items():
            val = fn(e)
            if isinstance(val, dict):
                # flatten dict into dotted sub-columns
                for k, v in val.items():
                    row[f"{name}.{k}"] = v
            else:
                row[name] = val

        rows.append(row)
    return rows


# -- public API ---------------------------------------------------------
def collect_summary_df(
    ds,
    native_columns: Collection[str] | None = None,
    extra_columns: Collection[Dict[str, Sequence[Any]]] | None = None,
    callables: Mapping[str, Callable] | None = None,
):
    """
    Build an overview table for `ds`.

    Parameters
    ----------
    ds : CrystalDataset
    native_columns : list/tuple of str
        Column names drawn directly from CrystalEntry attributes
        (`energy`, `id` ...) or metadata keys using `"metadata.<key>"`.
        If None → `("id", "energy")`.
    extra_columns : iterable of dicts `{column_name: sequence}`
        Pre-computed per-row vectors aligned with `ds` length.
    callables : dict `{column_name: fn(entry)}`
        Per-entry functions producing scalars or dicts. Dict results expand
        into dotted sub-columns `column.key`.

    Returns
    -------
    pandas.DataFrame
    """
    if _pd is None:  # pandas not installed
        raise RuntimeError(
            "collect_summary_df() requires pandas. "
            "Install with:  pip install crystaldata[table]"
        ) from _PANDAS_ERR

    native_cols = tuple(native_columns) if native_columns else ("id", "composition", "energy")
    callables = callables or {}
    extra_columns = extra_columns or ()

    # 1) build row dicts
    rows = _rows_from_dataset(ds, native_cols, callables)

    # 2) create DataFrame
    df = _pd.DataFrame(rows)

    # 3) attach extra columns
    for mapping in extra_columns:
        for name, values in mapping.items():
            if len(values) != len(ds):
                raise ValueError(
                    f"extra_columns[{name!r}] length {len(values)} "
                    f"!= dataset length {len(ds)}"
                )
            df[name] = list(values)

    return df
