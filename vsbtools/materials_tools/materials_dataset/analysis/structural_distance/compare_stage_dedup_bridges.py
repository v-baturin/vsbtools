from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Iterable

import numpy as np


FIXTURES_PROCESSED_ROOT = (
    Path(__file__).resolve().parents[2]
    / "unittests_datasets" / "MG_postprocess_pipelines" / "PROCESSED"
)

DEFAULT_STAGE_RELATIVE = Path(
    "B-Fe-Nd"
    / "B-Fe-Nd__guidance_environment_mode_huber_B-Fe_3__diffusion_loss_weight_0.5-0.5-True__algo_0"
    / "0_27fa706f8fc4ce3b"
)
DEFAULT_STAGE = FIXTURES_PROCESSED_ROOT / DEFAULT_STAGE_RELATIVE

DEFAULT_WORK_DIR = Path("/tmp/vsbtools_stage_dedup_bridge_compare")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")


def _ensure_repo_root_on_path() -> None:
    """Allow running this file directly, not only as `python -m ...`."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "vsbtools").is_dir():
            sys.path.insert(0, str(parent))
            return


_ensure_repo_root_on_path()

from vsbtools.materials_tools.materials_dataset.analysis.similarity_tools import SimilarityTools
from vsbtools.materials_tools.materials_dataset.analysis.structural_distance.dscribe_bridge import DScribeBridge
from vsbtools.materials_tools.materials_dataset.crystal_dataset import CrystalDataset
from vsbtools.materials_tools.materials_dataset.io.structures_dataset_io import StructureDatasetIO
from vsbtools.materials_tools.materials_dataset.io.uspex_bridge import USPEXBridge
from vsbtools.materials_tools.materials_dataset.io.yaml_csv_poscars import (
    read as read_dataset_manifest,
    write as write_dataset_manifest,
)
from vsbtools.materials_tools.materials_dataset.io.zip_handling import exploded_zip_tree


def _metadata_without_constructor_fields(ds: CrystalDataset) -> dict:
    metadata = dict(getattr(ds, "metadata", {}) or {})
    metadata.pop("message", None)
    metadata.pop("created_on", None)
    return metadata


def _clone_dataset(ds: CrystalDataset, suffix: str = "") -> CrystalDataset:
    entries = [entry.copy_with(metadata=dict(entry.metadata or {})) for entry in ds]
    dataset_id = f"{ds.dataset_id}{suffix}" if suffix else ds.dataset_id
    cloned = CrystalDataset(
        entries,
        dataset_id=dataset_id,
        parent_ids=list(getattr(ds, "parent_ids", []) or []),
        message=(getattr(ds, "metadata", {}) or {}).get("message", ""),
        supplementary_metadata=_metadata_without_constructor_fields(ds),
        base_path=ds.base_path,
    )
    if getattr(ds, "_elements", None) is not None:
        cloned._elements = set(ds.elements)
    return cloned


def _limit_dataset(ds: CrystalDataset, limit: int | None) -> CrystalDataset:
    if limit is None:
        return ds
    limited = CrystalDataset(
        list(ds)[:limit],
        dataset_id=f"{ds.dataset_id}_first{limit}",
        parent_ids=[ds.dataset_id],
        message=f"First {limit} entries from {ds.dataset_id}",
        supplementary_metadata=_metadata_without_constructor_fields(ds),
        base_path=ds.base_path,
    )
    if getattr(ds, "_elements", None) is not None:
        limited._elements = set(ds.elements)
    return limited


def load_stage(stage_path: Path, limit: int | None = None) -> CrystalDataset:
    manifest = stage_path / "manifest.yaml"
    if manifest.is_file():
        ds = read_dataset_manifest(manifest)
    else:
        poscars_dir = stage_path / "POSCARS"
        root = poscars_dir if poscars_dir.is_dir() else stage_path
        ds = StructureDatasetIO(root, pattern="*POSCAR*").load_from_directory()
        ds.override_base_path(stage_path)
    return _limit_dataset(ds, limit)


def _ids(entries: Iterable) -> list[str]:
    return [str(entry.id) for entry in entries]


def _cluster_ids(ds: CrystalDataset, clusters: list[list[int]]) -> set[tuple[str, ...]]:
    return {tuple(sorted(str(ds[i].id) for i in cluster)) for cluster in clusters}


def _write_lines(path: Path, lines: Iterable[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


@contextlib.contextmanager
def _maybe_silence_stdout(enabled: bool):
    if not enabled:
        yield
        return
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def deduplicate_with_bridge(
    ds: CrystalDataset,
    bridge_name: str,
    bridge,
    *,
    work_dir: Path,
    tol_fp: float,
    reuse_distance_matrix: bool,
    enforce_compositions_separation: bool,
    verbose_dist_progress: bool,
) -> tuple[CrystalDataset, list[list[int]], list[int], Path]:
    bridge_dir = work_dir / bridge_name
    bridge_dir.mkdir(parents=True, exist_ok=True)

    similarity = SimilarityTools(bridge.fp_dist, tol_fp)
    dist_matrix_file = bridge_dir / "dist_matrix.pkl"
    clusters_file = bridge_dir / "clusters.pkl"

    print(f"Running {bridge_name} deduplication for {len(ds)} entries")
    start = time.perf_counter()
    with _maybe_silence_stdout(not verbose_dist_progress):
        dedup_ds, clusters, best_idx = similarity.deduplicate(
            ds,
            check_clusters_file=False,
            check_dist_matrix_file=reuse_distance_matrix,
            clusters_file=clusters_file,
            dist_matrix_file=dist_matrix_file,
            tol_FP=tol_fp,
            enforce_compositions_separation=enforce_compositions_separation,
        )
    elapsed_s = time.perf_counter() - start

    representative_ids = [str(ds[i].id) for i in best_idx]
    _write_lines(bridge_dir / "representatives.txt", representative_ids)
    (bridge_dir / "timing.json").write_text(json.dumps({"elapsed_s": elapsed_s}, indent=2) + "\n")
    print(f"{bridge_name}: {len(representative_ids)} representatives in {elapsed_s:.2f} s")
    return dedup_ds, clusters, best_idx, dist_matrix_file


def compare_distance_matrices(
    uspex_matrix_file: Path,
    dscribe_matrix_file: Path,
    *,
    ids: list[str],
    tol_fp: float,
    max_pairs: int,
) -> dict:
    with open(uspex_matrix_file, "rb") as handle:
        uspex_matrix = pickle.load(handle)
    with open(dscribe_matrix_file, "rb") as handle:
        dscribe_matrix = pickle.load(handle)

    delta = np.abs(uspex_matrix - dscribe_matrix)
    if delta.shape[0] > 1:
        pair_delta = np.triu(delta, k=1)
        max_i, max_j = np.unravel_index(np.argmax(pair_delta), pair_delta.shape)
        max_delta_pair = {
            "i": int(max_i),
            "j": int(max_j),
            "id_i": ids[max_i],
            "id_j": ids[max_j],
            "uspex": float(uspex_matrix[max_i, max_j]),
            "dscribe": float(dscribe_matrix[max_i, max_j]),
        }
        max_abs_delta = float(pair_delta[max_i, max_j])
    else:
        max_delta_pair = None
        max_abs_delta = 0.0

    uspex_adj = uspex_matrix <= tol_fp
    dscribe_adj = dscribe_matrix <= tol_fp
    mismatch_mask = np.triu(uspex_adj != dscribe_adj, k=1)
    mismatch_i, mismatch_j = np.where(mismatch_mask)

    pairs = []
    for i, j in zip(mismatch_i[:max_pairs], mismatch_j[:max_pairs]):
        pairs.append(
            {
                "i": int(i),
                "j": int(j),
                "id_i": ids[i],
                "id_j": ids[j],
                "uspex": float(uspex_matrix[i, j]),
                "dscribe": float(dscribe_matrix[i, j]),
            }
        )

    return {
        "max_abs_distance_delta": max_abs_delta,
        "max_abs_distance_delta_pair": max_delta_pair,
        "n_duplicate_decision_mismatches": int(len(mismatch_i)),
        "duplicate_decision_mismatch_examples": pairs,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Deduplicate a processed stage with USPEXBridge and DScribeBridge, "
            "then assert that representative ID sets are identical."
        )
    )
    parser.add_argument("stage", nargs="?", type=Path, default=DEFAULT_STAGE)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--tol-fp", type=float, default=0.012)
    parser.add_argument("--limit", type=int, default=None, help="Optional first-N subset for smoke testing.")
    parser.add_argument(
        "--reuse-distance-matrix",
        action="store_true",
        help="Reuse existing bridge-specific distance matrices from --work-dir if present.",
    )
    parser.add_argument(
        "--no-composition-separation",
        action="store_true",
        help="Do not split duplicate clusters by reduced composition.",
    )
    parser.add_argument(
        "--verbose-dist-progress",
        action="store_true",
        help="Let the underlying distance-matrix progress prints through.",
    )
    parser.add_argument(
        "--dump-deduplicated",
        action="store_true",
        help="Write both deduplicated datasets under --work-dir after the equality check.",
    )
    parser.add_argument(
        "--fingerprint-cache-size",
        type=int,
        default=None,
        help=(
            "Maximum number of DScribe fingerprints to retain in memory. "
            "Default keeps all lazily computed fingerprints; 0 disables caching."
        ),
    )
    parser.add_argument("--max-mismatch-pairs", type=int, default=20)
    return parser.parse_args()


def _run_for_stage(args, stage_path: Path, work_dir: Path) -> int:
    base_ds = load_stage(stage_path, args.limit)
    elements = sorted(base_ds.elements)
    ids = _ids(base_ds)
    print(f"Loaded {len(base_ds)} entries from {stage_path}")
    print(f"Elements: {', '.join(elements)}")
    print(f"tol_FP = {args.tol_fp}")
    print(f"Work dir: {work_dir}")

    uspex_bridge = USPEXBridge(elements=set(elements), legacy=True, tol_FP=args.tol_fp)
    dscribe_bridge = DScribeBridge(
        elements=elements,
        preset="uspex",
        Rmax=uspex_bridge.rdu.Rmax,
        sigma=uspex_bridge.rdu.sigma,
        delta=uspex_bridge.rdu.delta,
        tol_FP=args.tol_fp,
        fingerprint_cache_size=args.fingerprint_cache_size,
    )

    enforce_compositions = not args.no_composition_separation
    uspex_ds, uspex_clusters, uspex_best_idx, uspex_matrix_file = deduplicate_with_bridge(
        _clone_dataset(base_ds, "_uspex_input"),
        "uspex",
        uspex_bridge,
        work_dir=work_dir,
        tol_fp=args.tol_fp,
        reuse_distance_matrix=args.reuse_distance_matrix,
        enforce_compositions_separation=enforce_compositions,
        verbose_dist_progress=args.verbose_dist_progress,
    )
    dscribe_ds, dscribe_clusters, dscribe_best_idx, dscribe_matrix_file = deduplicate_with_bridge(
        _clone_dataset(base_ds, "_dscribe_input"),
        "dscribe",
        dscribe_bridge,
        work_dir=work_dir,
        tol_fp=args.tol_fp,
        reuse_distance_matrix=args.reuse_distance_matrix,
        enforce_compositions_separation=enforce_compositions,
        verbose_dist_progress=args.verbose_dist_progress,
    )

    uspex_representatives = [ids[i] for i in uspex_best_idx]
    dscribe_representatives = [ids[i] for i in dscribe_best_idx]
    uspex_set = set(uspex_representatives)
    dscribe_set = set(dscribe_representatives)
    only_uspex = sorted(uspex_set - dscribe_set)
    only_dscribe = sorted(dscribe_set - uspex_set)
    same_representative_set = not only_uspex and not only_dscribe
    same_representative_order = uspex_representatives == dscribe_representatives
    same_clusters = _cluster_ids(base_ds, uspex_clusters) == _cluster_ids(base_ds, dscribe_clusters)

    matrix_summary = compare_distance_matrices(
        uspex_matrix_file,
        dscribe_matrix_file,
        ids=ids,
        tol_fp=args.tol_fp,
        max_pairs=args.max_mismatch_pairs,
    )
    uspex_timing = json.loads((work_dir / "uspex" / "timing.json").read_text())
    dscribe_timing = json.loads((work_dir / "dscribe" / "timing.json").read_text())

    if args.dump_deduplicated and same_representative_set:
        write_dataset_manifest(uspex_ds, enforce_base_path=work_dir / "uspex" / "deduplicated")
        write_dataset_manifest(dscribe_ds, enforce_base_path=work_dir / "dscribe" / "deduplicated")

    _write_lines(work_dir / "only_uspex.txt", only_uspex)
    _write_lines(work_dir / "only_dscribe.txt", only_dscribe)

    summary = {
        "stage": stage_path.as_posix(),
        "n_entries": len(base_ds),
        "elements": elements,
        "tol_fp": args.tol_fp,
        "composition_separation": enforce_compositions,
        "uspex_n_representatives": len(uspex_representatives),
        "dscribe_n_representatives": len(dscribe_representatives),
        "uspex_elapsed_s": uspex_timing["elapsed_s"],
        "dscribe_elapsed_s": dscribe_timing["elapsed_s"],
        "dscribe_to_uspex_time_ratio": dscribe_timing["elapsed_s"] / uspex_timing["elapsed_s"],
        "same_representative_set": same_representative_set,
        "same_representative_order": same_representative_order,
        "same_clusters": same_clusters,
        "only_uspex": only_uspex,
        "only_dscribe": only_dscribe,
        **matrix_summary,
    }
    (work_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"same representative set: {same_representative_set}")
    print(f"USPEX elapsed: {summary['uspex_elapsed_s']:.2f} s")
    print(f"DScribe elapsed: {summary['dscribe_elapsed_s']:.2f} s")
    print(f"DScribe/USPEX time ratio: {summary['dscribe_to_uspex_time_ratio']:.3f}")
    print(f"same representative order: {same_representative_order}")
    print(f"same clusters: {same_clusters}")
    print(f"duplicate-decision mismatches: {summary['n_duplicate_decision_mismatches']}")
    print(f"max |USPEX - DScribe| distance delta: {summary['max_abs_distance_delta']:.6g}")
    print(f"Summary written to {work_dir / 'summary.json'}")

    if not same_representative_set:
        print(f"Only in USPEX representatives: {only_uspex[:args.max_mismatch_pairs]}")
        print(f"Only in DScribe representatives: {only_dscribe[:args.max_mismatch_pairs]}")
        return 1
    return 0


def main() -> int:
    args = parse_args()
    work_dir = args.work_dir.expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    if args.stage == DEFAULT_STAGE and not DEFAULT_STAGE.exists():
        with exploded_zip_tree(FIXTURES_PROCESSED_ROOT) as processed_root:
            return _run_for_stage(args, processed_root / DEFAULT_STAGE_RELATIVE, work_dir)

    return _run_for_stage(args, args.stage.expanduser().resolve(), work_dir)


if __name__ == "__main__":
    raise SystemExit(main())
