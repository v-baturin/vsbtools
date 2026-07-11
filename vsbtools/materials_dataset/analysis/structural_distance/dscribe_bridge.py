from __future__ import annotations

from collections import OrderedDict
from typing import Callable, Iterable, Any

import numpy as np
from pymatgen.core.periodic_table import Element

from ...crystal_entry import CrystalEntry

USPEX_RMAX_DEFAULT = 10.0
USPEX_SIGMA_DEFAULT = 0.03
USPEX_DELTA_DEFAULT = 0.08
USPEX_TOLERANCE_DEFAULT = 0.012


def _species_by_atomic_number(elements: Iterable[str]) -> list[str]:
    return sorted({str(el) for el in elements}, key=lambda symbol: Element(symbol).Z)


def _as_dense_vector(values: Any) -> np.ndarray:
    if hasattr(values, "toarray"):
        values = values.toarray()
    elif hasattr(values, "todense"):
        values = values.todense()
    return np.asarray(values, dtype=float).reshape(-1)


def cosine_distance(vec1: np.ndarray, vec2: np.ndarray) -> float:
    denom = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    if denom == 0:
        return 0.5
    return float(0.5 * (1.0 - np.dot(vec1, vec2) / denom))


def dscribe_uspex_mbtr_kwargs(
        species: Iterable[str],
        Rmax: float = USPEX_RMAX_DEFAULT,
        sigma: float = USPEX_SIGMA_DEFAULT,
        delta: float = USPEX_DELTA_DEFAULT,
        **kwargs
) -> dict[str, Any]:
    """
    DScribe MBTR parameters that reproduce USPEX structure fingerprint channels.

    USPEX bins intervals [k*delta, (k+1)*delta]. DScribe MBTR samples grid
    points, so the equivalent grid uses the USPEX bin centers.
    """
    return {
        "species": list(species),
        "periodic": True,
        "geometry": {"function": "distance"},
        "grid": {
            "min": delta / 2,
            "max": Rmax - delta / 2,
            "sigma": sigma,
            "n": int(round(Rmax / delta)),
        },
        "weighting": {"function": "inverse_square", "r_cut": Rmax},
        "normalization": "valle_oganov",
        "normalize_gaussians": True,
        "sparse": False,
        **kwargs,
    }


class DScribeBridge:
    """
    Similarity bridge backed by DScribe descriptors.

    By default this uses a DScribe MBTR parametrization matching the USPEX
    radial-distribution structure fingerprint. For arbitrary DScribe workflows,
    pass either a descriptor instance or a descriptor class/name with kwargs.
    """

    def __init__(
            self,
            elements: Iterable[str] | None = None,
            *,
            preset: str | None = "uspex",
            descriptor: Any | None = None,
            descriptor_cls: type | str | None = None,
            descriptor_kwargs: dict[str, Any] | None = None,
            create_kwargs: dict[str, Any] | None = None,
            feature_transform: Callable[[np.ndarray, CrystalEntry, Any], np.ndarray] | None = None,
            metric: str | Callable[[np.ndarray, np.ndarray], float] = "cosine",
            tol_FP: float = USPEX_TOLERANCE_DEFAULT,
            fingerprint_cache_size: int | None = None,
            **preset_kwargs,
    ):
        if preset == "uspex" and (descriptor is not None or descriptor_cls is not None):
            preset = None

        self.species = _species_by_atomic_number(elements or [])
        self.preset = preset
        self.create_kwargs = dict(create_kwargs or {})
        self.feature_transform = feature_transform
        self.tol_FP = tol_FP

        if preset == "uspex":
            descriptor = descriptor or self._make_uspex_descriptor(**preset_kwargs)
        elif descriptor is None:
            descriptor = self._make_custom_descriptor(descriptor_cls, descriptor_kwargs)
        self.descriptor = descriptor
        self.metric = metric
        self.fingerprint_cache_size = fingerprint_cache_size
        self._fingerprint_cache: OrderedDict[tuple[int, str], np.ndarray] = OrderedDict()
        self._sig = (
            tuple(self.species),
            preset,
            repr(self.descriptor),
            repr(self.create_kwargs),
            metric if isinstance(metric, str) else repr(metric),
            tol_FP,
            fingerprint_cache_size,
        )

    def __hash__(self):
        return hash(self._sig)

    def __eq__(self, other):
        return isinstance(other, DScribeBridge) and self._sig == other._sig

    @staticmethod
    def _make_custom_descriptor(descriptor_cls, descriptor_kwargs):
        if descriptor_cls is None:
            raise ValueError("Provide descriptor, descriptor_cls, or preset='uspex'")
        descriptor_kwargs = dict(descriptor_kwargs or {})
        if isinstance(descriptor_cls, str):
            from dscribe import descriptors
            descriptor_cls = getattr(descriptors, descriptor_cls)
        return descriptor_cls(**descriptor_kwargs)

    def _make_uspex_descriptor(self, **preset_kwargs):
        if not self.species:
            raise ValueError("elements must be provided for preset='uspex'")

        Rmax = preset_kwargs.pop("Rmax", preset_kwargs.pop("r_max", USPEX_RMAX_DEFAULT))
        sigma = preset_kwargs.pop("sigma", USPEX_SIGMA_DEFAULT)
        delta = preset_kwargs.pop("delta", USPEX_DELTA_DEFAULT)

        from dscribe.descriptors import MBTR
        return MBTR(**dscribe_uspex_mbtr_kwargs(
            self.species,
            Rmax=Rmax,
            sigma=sigma,
            delta=delta,
            **preset_kwargs,
        ))

    def raw_descriptor(self, entry: CrystalEntry) -> np.ndarray:
        return _as_dense_vector(self.descriptor.create(entry.structure.to_ase_atoms(), **self.create_kwargs))

    def uspex_fingerprint_channels(self, entry: CrystalEntry) -> dict[tuple[str, str], np.ndarray]:
        if self.preset != "uspex":
            raise RuntimeError("uspex_fingerprint_channels is only available for preset='uspex'")

        raw = self.raw_descriptor(entry)
        volume = entry.structure.volume
        return {
            (a, b): (raw[self.descriptor.get_location((a, b))] - 1.0) / volume
            for i, a in enumerate(self.species)
            for b in self.species[i:]
        }

    def _uspex_weighted_vector(self, entry: CrystalEntry) -> np.ndarray:
        channels = self.uspex_fingerprint_channels(entry)
        counts = {el: 0 for el in self.species}
        for site in entry.structure:
            symbol = site.specie.symbol
            if symbol in counts:
                counts[symbol] += 1

        n_atoms = len(entry.structure)
        parts = []
        for i, a in enumerate(self.species):
            for b in self.species[i:]:
                if a == b:
                    weight = counts[a] * counts[b] / n_atoms ** 2
                else:
                    weight = 2 * counts[a] * counts[b] / n_atoms ** 2
                parts.append(np.sqrt(weight) * channels[(a, b)])
        return np.concatenate(parts)

    def fingerprint(self, entry: CrystalEntry) -> np.ndarray:
        key = (id(entry), str(entry.id))
        if key in self._fingerprint_cache:
            self._fingerprint_cache.move_to_end(key)
            return self._fingerprint_cache[key]

        if self.preset == "uspex":
            vector = self._uspex_weighted_vector(entry)
            self._cache_fingerprint(key, vector)
            return vector

        vector = self.raw_descriptor(entry)
        if self.feature_transform is not None:
            vector = self.feature_transform(vector, entry, self.descriptor)
        vector = _as_dense_vector(vector)
        self._cache_fingerprint(key, vector)
        return vector

    def _cache_fingerprint(self, key: tuple[int, str], vector: np.ndarray) -> None:
        if self.fingerprint_cache_size == 0:
            return
        self._fingerprint_cache[key] = vector
        self._fingerprint_cache.move_to_end(key)
        if self.fingerprint_cache_size is not None:
            while len(self._fingerprint_cache) > self.fingerprint_cache_size:
                self._fingerprint_cache.popitem(last=False)

    def clear_cache(self) -> None:
        self._fingerprint_cache.clear()

    def dist(self, entry_1: CrystalEntry, entry_2: CrystalEntry) -> float:
        fp1 = self.fingerprint(entry_1)
        fp2 = self.fingerprint(entry_2)

        if callable(self.metric):
            return float(self.metric(fp1, fp2))
        if self.metric == "cosine":
            return cosine_distance(fp1, fp2)
        if self.metric in ("euclidean", "l2"):
            return float(np.linalg.norm(fp1 - fp2))
        raise ValueError(f"Unknown metric: {self.metric!r}")

    def fp_dist(self, entry_1: CrystalEntry, entry_2: CrystalEntry) -> float:
        return self.dist(entry_1, entry_2)
