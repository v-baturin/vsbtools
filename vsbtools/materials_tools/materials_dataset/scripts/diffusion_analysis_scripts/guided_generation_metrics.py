from __future__ import annotations
import ast
import json
import re
import shutil
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence
import pandas as pd
from pymatgen.core import Element
from ...io import read
from ...io.yaml_csv_poscars import load_yaml_recursively
from ...analysis.phase_diagram_tools import PhaseDiagramTools
from ...analysis import summary as summary_tools
from ...analysis import symmetry_tools as symmetry_tools_module
from .guidance_stats import callables_from_ds, get_two_proportion_z_test
MAX_ECH = 0.2
DEFAULT_MAX_DEVIATION = {
    "environment": 0.2,
    "dominant_environment": 0.1,
    "volume_pa": 0.2,
}
DEFAULT_POLL_DB_STAGE_CANDIDATES = ("poll_db", "poll_db_grace")
DEFAULT_DEDUP_STAGE_CANDIDATES = ("deduplicate_all", "dedup_relaxed_grace")
DEFAULT_ADD_REF_STAGE_CANDIDATES = (
    "add_ref_add_ref_deduplicated",
    "add_ref_deduplicated",
    "add_ref_dedup_relax_grace",
)
@dataclass(slots=True)
class GuidanceMetricReport:
    eval_repo_dir: Path
    spec_repo_dir: Path
    parse_raw_stage: str
    dedup_stage: str
    poll_db_stage: str
    add_ref_stage: str
    elements: tuple[str, ...]
    entries_df: pd.DataFrame
    counts: dict[str, int]
    ratios: dict[str, float | None]
    summary: dict[str, Any]
    stage_dir: Path
    stage_summary_path: Path | None = None
    @property
    def good_ids(self) -> set[str]:
        return set(self.entries_df.loc[self.entries_df["in_good"], "id"].astype(str))
@dataclass(slots=True)
class GuidanceComparison:
    guided: GuidanceMetricReport
    nonguided: GuidanceMetricReport
    pvalue_summary: dict[str, Any]
def normalize_elements(elements: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted(str(Element(str(el))) for el in elements))
def entry_has_exact_elements(entry, elements: Iterable[str]) -> bool:
    return {str(el) for el in entry.composition.elements} == set(normalize_elements(elements))
def entry_n_atoms(entry) -> int:
    return int(round(float(entry.composition.num_atoms)))
def parse_duplicates(value) -> set[str]:
    if value in (None, "NA", "None", "set()"):
        return set()
    if isinstance(value, set):
        return {str(v) for v in value}
    if isinstance(value, (list, tuple)):
        return {str(v) for v in value}
    if isinstance(value, str):
        text = value.strip()
        if not text or text == "set()":
            return set()
        try:
            parsed = ast.literal_eval(text.replace("' '", "','"))
        except Exception:
            return set(re.findall(r"[A-Za-z0-9_.:-]+", text))
        return parse_duplicates(parsed)
    return {str(value)}
def stage_manifest_map(repo_dir: str | Path) -> dict[str, Path]:
    repo_dir = Path(repo_dir)
    manifests: dict[str, Path] = {}
    for manifest in repo_dir.rglob("manifest.yaml"):
        info = load_yaml_recursively(manifest)
        stage = info.get("metadata", {}).get("pipeline_stage")
        if stage is not None:
            manifests[str(stage)] = manifest
    return manifests
def resolve_stage_manifest(repo_dir: str | Path, *stage_candidates: str, required: bool = True) -> Path | None:
    manifests = stage_manifest_map(repo_dir)
    for stage in stage_candidates:
        if stage in manifests:
            return manifests[stage]
    if required:
        raise KeyError(f"None of the stage candidates {stage_candidates} found in {repo_dir}")
    return None
def load_stage(repo_dir: str | Path, *stage_candidates: str, required: bool = True):
    manifest = resolve_stage_manifest(repo_dir, *stage_candidates, required=required)
    return None if manifest is None else read(manifest)
def resolve_stage_dir(repo_dir: str | Path, *stage_candidates: str, required: bool = True) -> Path | None:
    manifest = resolve_stage_manifest(repo_dir, *stage_candidates, required=required)
    return None if manifest is None else manifest.parent
def resolve_nonguided_repo_dir(guided_repo_dir: str | Path, search_root: str | Path | None = None) -> Path:
    guided_repo_dir = Path(guided_repo_dir)
    parse_raw_ds = load_stage(guided_repo_dir, "parse_raw")
    batch_metadata = parse_raw_ds.metadata.get("batch_metadata", {})
    candidate_names = [
        batch_metadata.get("nonguided_repo"),
        batch_metadata.get("git_info", {}).get("nonguided_repo"),
    ]
    candidate_names = [name for name in candidate_names if name]
    if search_root is None:
        search_root = guided_repo_dir.parent.parent
    search_root = Path(search_root)

    if candidate_names:
        nonguided_name = candidate_names[0]
        matches = list(search_root.glob(f"*/{nonguided_name}"))
        if not matches:
            direct = search_root / nonguided_name
            if direct.exists():
                return direct
            raise FileNotFoundError(f"Could not resolve nonguided repo dir for {nonguided_name} under {search_root}")
        if len(matches) > 1:
            raise RuntimeError(f"Multiple nonguided repo dirs found for {nonguided_name}: {matches}")
        return matches[0]

    system_prefix = guided_repo_dir.name.split("__guidance_", 1)[0]
    sibling_matches = sorted(guided_repo_dir.parent.glob(f"{system_prefix}__guidance_None__*"))
    if len(sibling_matches) == 1:
        return sibling_matches[0]
    if len(sibling_matches) > 1:
        preferred = [path for path in sibling_matches if "__diffusion_loss_weight_1__" in path.name]
        if len(preferred) == 1:
            return preferred[0]
        raise RuntimeError(
            f"Multiple nonguided sibling repo dirs found for {guided_repo_dir}: {sibling_matches}. "
            "Provide nonguided_repo_dir explicitly."
        )
    raise KeyError(
        f"Could not infer nonguided repo dir from metadata or sibling repo names for {guided_repo_dir}"
    )
def guidance_spec_from_repo(repo_dir: str | Path, max_deviation: dict[str, float] | None = None) -> dict[str, Any]:
    max_deviation = max_deviation or DEFAULT_MAX_DEVIATION
    parse_raw_ds = load_stage(repo_dir, "parse_raw")
    callables, targets, guidance_name = callables_from_ds(parse_raw_ds, include_losses=False)
    if guidance_name not in max_deviation:
        raise KeyError(f"No max deviation configured for guidance type {guidance_name!r}")
    margins = {key: max_deviation[guidance_name] for key in callables}
    return {
        "parse_raw_ds": parse_raw_ds,
        "callables": callables,
        "targets": targets,
        "margins": margins,
        "guidance_name": guidance_name,
    }
def satisfies_guidance(entry, callables, targets, margins) -> bool:
    for key, fn in callables.items():
        if abs(fn(entry) - targets[key]) > margins[key]:
            return False
    return True
def known_duplicate_ids_from_reference(add_ref_ds, elements: Iterable[str]) -> set[str]:
    known_ids: set[str] = set()
    for entry in add_ref_ds:
        if not entry_has_exact_elements(entry, elements):
            continue
        if re.search(r"[a-z]", str(entry.id)) is None:
            continue
        known_ids |= parse_duplicates(entry.metadata.get("duplicates"))
    return known_ids
def _safe_ratio(numerator: int, denominator: int) -> float | None:
    if denominator == 0:
        return None
    return numerator / denominator
def _collect_known_denominator_count(
    poll_db_ds,
    spec: dict[str, Any],
    elements: Sequence[str],
    max_ech: float,
    denominator_max_atoms: int | None,
) -> int:
    pd_tools = PhaseDiagramTools(poll_db_ds)
    count = 0
    for entry in poll_db_ds:
        if not entry_has_exact_elements(entry, elements):
            continue
        if denominator_max_atoms is not None and entry_n_atoms(entry) > denominator_max_atoms:
            continue
        e_ch = pd_tools.height_above_hull_pa(entry)
        if e_ch >= max_ech:
            continue
        if not satisfies_guidance(entry, spec["callables"], spec["targets"], spec["margins"]):
            continue
        count += 1
    return count
def collect_guidance_metric_report(
    eval_repo_dir: str | Path,
    spec_repo_dir: str | Path | None = None,
    *,
    elements: Iterable[str] | None = None,
    max_ech: float = MAX_ECH,
    max_deviation: dict[str, float] | None = None,
    denominator_max_atoms: int | None = 20,
    poll_db_stage_candidates: Sequence[str] = DEFAULT_POLL_DB_STAGE_CANDIDATES,
    dedup_stage_candidates: Sequence[str] = DEFAULT_DEDUP_STAGE_CANDIDATES,
    add_ref_stage_candidates: Sequence[str] = DEFAULT_ADD_REF_STAGE_CANDIDATES,
) -> GuidanceMetricReport:
    eval_repo_dir = Path(eval_repo_dir)
    spec_repo_dir = eval_repo_dir if spec_repo_dir is None else Path(spec_repo_dir)
    spec = guidance_spec_from_repo(spec_repo_dir, max_deviation=max_deviation)
    parse_raw_ds = load_stage(eval_repo_dir, "parse_raw")
    elements = normalize_elements(elements or spec["parse_raw_ds"].elements)
    poll_db_manifest = resolve_stage_manifest(eval_repo_dir, *poll_db_stage_candidates)
    dedup_manifest = resolve_stage_manifest(eval_repo_dir, *dedup_stage_candidates)
    add_ref_manifest = resolve_stage_manifest(eval_repo_dir, *add_ref_stage_candidates)
    poll_db_ds = read(poll_db_manifest)
    dedup_ds = read(dedup_manifest)
    add_ref_ds = read(add_ref_manifest)
    pd_tools = PhaseDiagramTools(poll_db_ds)
    known_ids = known_duplicate_ids_from_reference(add_ref_ds, elements)
    rows: list[dict[str, Any]] = []
    for entry in dedup_ds:
        if not entry_has_exact_elements(entry, elements):
            continue
        e_ch = pd_tools.height_above_hull_pa(entry)
        is_stable = e_ch < max_ech
        is_new = str(entry.id) not in known_ids
        passes_constraint = satisfies_guidance(entry, spec["callables"], spec["targets"], spec["margins"])
        rows.append(
            {
                "id": str(entry.id),
                "composition": entry.composition.reduced_formula,
                "e_ch": e_ch,
                "in_U": True,
                "in_S": is_stable,
                "in_N": is_new,
                "in_C": passes_constraint,
                "in_US": is_stable,
                "in_USN": is_stable and is_new,
                "in_good": is_stable and is_new and passes_constraint,
            }
        )
    entries_df = pd.DataFrame(rows)
    if len(entries_df):
        entries_df = entries_df.sort_values(["in_good", "e_ch", "id"], ascending=[False, True, True]).reset_index(drop=True)
    T = len(parse_raw_ds)
    U = len(entries_df)
    S = int(entries_df["in_S"].sum()) if len(entries_df) else 0
    N = int(entries_df["in_N"].sum()) if len(entries_df) else 0
    C = int(entries_df["in_C"].sum()) if len(entries_df) else 0
    US = int(entries_df["in_US"].sum()) if len(entries_df) else 0
    USN = int(entries_df["in_USN"].sum()) if len(entries_df) else 0
    USNC = int(entries_df["in_good"].sum()) if len(entries_df) else 0
    known_denominator = _collect_known_denominator_count(
        poll_db_ds,
        spec,
        elements,
        max_ech=max_ech,
        denominator_max_atoms=denominator_max_atoms,
    )
    counts = {
        "T": T,
        "U": U,
        "S": S,
        "N": N,
        "C": C,
        "US": US,
        "USN": USN,
        "USNC": USNC,
        "known_denominator": known_denominator,
    }
    ratios = {
        "|U|/|T|": _safe_ratio(U, T),
        "|U∩S|/|U|": _safe_ratio(US, U),
        "|U∩S∩N|/|U∩S|": _safe_ratio(USN, US),
        "|U∩S∩N∩C|/|U∩S∩N|": _safe_ratio(USNC, USN),
        "good_per_raw": _safe_ratio(USNC, T),
        "scout": _safe_ratio(USNC, known_denominator),
    }
    summary = {
        "eval_repo_dir": str(eval_repo_dir),
        "spec_repo_dir": str(spec_repo_dir),
        "elements": list(elements),
        "guidance_name": spec["guidance_name"],
        "targets": spec["targets"],
        "margins": spec["margins"],
        "max_ech": max_ech,
        "denominator_max_atoms": denominator_max_atoms,
        "parse_raw_stage": "parse_raw",
        "dedup_stage": dedup_manifest.parent.name,
        "poll_db_stage": poll_db_manifest.parent.name,
        "add_ref_stage": add_ref_manifest.parent.name,
        **counts,
        **ratios,
    }
    return GuidanceMetricReport(
        eval_repo_dir=eval_repo_dir,
        spec_repo_dir=spec_repo_dir,
        parse_raw_stage="parse_raw",
        dedup_stage=str(load_yaml_recursively(dedup_manifest).get("metadata", {}).get("pipeline_stage")),
        poll_db_stage=str(load_yaml_recursively(poll_db_manifest).get("metadata", {}).get("pipeline_stage")),
        add_ref_stage=str(load_yaml_recursively(add_ref_manifest).get("metadata", {}).get("pipeline_stage")),
        elements=elements,
        entries_df=entries_df,
        counts=counts,
        ratios=ratios,
        summary=summary,
        stage_dir=add_ref_manifest.parent,
    )
def compare_guided_and_nonguided(
    guided_repo_dir: str | Path,
    nonguided_repo_dir: str | Path | None = None,
    *,
    spec_repo_dir: str | Path | None = None,
    stage: str = "symmetrize_raw",
    elements: Iterable[str] | None = None,
    max_ech: float = MAX_ECH,
    max_deviation: dict[str, float] | None = None,
    denominator_max_atoms: int | None = 20,
    poll_db_stage_candidates: Sequence[str] = DEFAULT_POLL_DB_STAGE_CANDIDATES,
    dedup_stage_candidates: Sequence[str] = DEFAULT_DEDUP_STAGE_CANDIDATES,
    add_ref_stage_candidates: Sequence[str] = DEFAULT_ADD_REF_STAGE_CANDIDATES,
) -> GuidanceComparison:
    guided_repo_dir = Path(guided_repo_dir)
    spec_repo_dir = guided_repo_dir if spec_repo_dir is None else Path(spec_repo_dir)
    nonguided_repo_dir = resolve_nonguided_repo_dir(guided_repo_dir) if nonguided_repo_dir is None else Path(nonguided_repo_dir)
    spec = guidance_spec_from_repo(spec_repo_dir, max_deviation=max_deviation)
    elements = normalize_elements(elements or spec["parse_raw_ds"].elements)
    guided_stage_ds = load_stage(guided_repo_dir, stage)
    nonguided_stage_ds = load_stage(nonguided_repo_dir, stage)
    def count_satisfying(ds) -> tuple[int, int]:
        eligible = [entry for entry in ds if entry_has_exact_elements(entry, elements)]
        successes = sum(
            int(satisfies_guidance(entry, spec["callables"], spec["targets"], spec["margins"]))
            for entry in eligible
        )
        return successes, len(eligible)
    guided_successes, guided_total = count_satisfying(guided_stage_ds)
    nonguided_successes, nonguided_total = count_satisfying(nonguided_stage_ds)
    z_stats = get_two_proportion_z_test(
        nonguided_stage_ds,
        spec["callables"],
        spec["targets"],
        spec["margins"],
        ds_guided=guided_stage_ds,
        alternative="greater",
    )
    guided_report = collect_guidance_metric_report(
        guided_repo_dir,
        spec_repo_dir=spec_repo_dir,
        elements=elements,
        max_ech=max_ech,
        max_deviation=max_deviation,
        denominator_max_atoms=denominator_max_atoms,
        poll_db_stage_candidates=poll_db_stage_candidates,
        dedup_stage_candidates=dedup_stage_candidates,
        add_ref_stage_candidates=add_ref_stage_candidates,
    )
    nonguided_report = collect_guidance_metric_report(
        nonguided_repo_dir,
        spec_repo_dir=spec_repo_dir,
        elements=elements,
        max_ech=max_ech,
        max_deviation=max_deviation,
        denominator_max_atoms=denominator_max_atoms,
        poll_db_stage_candidates=poll_db_stage_candidates,
        dedup_stage_candidates=dedup_stage_candidates,
        add_ref_stage_candidates=add_ref_stage_candidates,
    )
    good_per_raw_ratio = None
    if guided_report.ratios["good_per_raw"] is not None and nonguided_report.ratios["good_per_raw"] not in (None, 0):
        good_per_raw_ratio = guided_report.ratios["good_per_raw"] / nonguided_report.ratios["good_per_raw"]
    pvalue_summary = {
        "stage": stage,
        "guided_successes": guided_successes,
        "guided_total": guided_total,
        "nonguided_successes": nonguided_successes,
        "nonguided_total": nonguided_total,
        "z_score": z_stats["z_score"],
        "p_value": z_stats["p_value"],
        "good_per_raw_ratio": good_per_raw_ratio,
    }
    return GuidanceComparison(
        guided=guided_report,
        nonguided=nonguided_report,
        pvalue_summary=pvalue_summary,
    )
def build_stage_summary_dataframe(repo_dir: str | Path, stage: str) -> pd.DataFrame:
    repo_dir = Path(repo_dir)
    parse_raw_ds = load_stage(repo_dir, "parse_raw")
    stage_ds = load_stage(repo_dir, stage)
    poll_db_ds = load_stage(repo_dir, *DEFAULT_POLL_DB_STAGE_CANDIDATES)
    callables, _targets, _guidance_name = callables_from_ds(parse_raw_ds)
    stk = symmetry_tools_module.SymmetryToolkit(a_sym_prec=1e-3)
    callables["symmetry"] = stk.sym_group_symbol
    pd_tools = PhaseDiagramTools(poll_db_ds)
    callables["e_ch"] = pd_tools.height_above_hull_pa
    return summary_tools.collect_summary_df(
        stage_ds,
        native_columns=("id", "composition", "energy"),
        callables=callables,
    )
def load_or_build_stage_summary(repo_dir: str | Path, stage: str) -> tuple[pd.DataFrame, Path, Path | None]:
    stage_dir = resolve_stage_dir(repo_dir, stage)
    summary_csv = stage_dir / "summary.csv"
    if summary_csv.exists():
        summary_df = pd.read_csv(summary_csv)
        unnamed = [col for col in summary_df.columns if str(col).startswith("Unnamed:")]
        if unnamed:
            summary_df = summary_df.drop(columns=unnamed)
        return summary_df, stage_dir, summary_csv
    warnings.warn(f"{summary_csv} not found. Rebuilding stage summary before exporting good structures.")
    summary_df = build_stage_summary_dataframe(repo_dir, stage=stage)
    return summary_df, stage_dir, None
def write_guidance_metric_artifacts(
    report: GuidanceMetricReport,
    *,
    output_dir: str | Path | None = None,
    include_good_poscars: bool = True,
    summary_stage_candidates: Sequence[str] = DEFAULT_ADD_REF_STAGE_CANDIDATES,
) -> dict[str, Path]:
    output_dir = report.stage_dir if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    intersections_csv = output_dir / "intersection_membership.csv"
    summary_json = output_dir / "intersection_summary.json"
    summary_csv = output_dir / "intersection_summary.csv"
    report.entries_df.to_csv(intersections_csv, index=False)
    summary_json.write_text(json.dumps(report.summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    pd.DataFrame([report.summary]).to_csv(summary_csv, index=False)
    written = {
        "intersection_membership_csv": intersections_csv,
        "intersection_summary_json": summary_json,
        "intersection_summary_csv": summary_csv,
    }
    summary_df = None
    stage_dir = None
    source_summary_path = None
    for stage_name in summary_stage_candidates:
        try:
            summary_df, stage_dir, source_summary_path = load_or_build_stage_summary(report.eval_repo_dir, stage=stage_name)
            break
        except KeyError:
            continue
    if summary_df is None or stage_dir is None:
        raise KeyError(f"None of {summary_stage_candidates} stages found in {report.eval_repo_dir}")
    if "id" not in summary_df.columns:
        raise KeyError(f"id column not found in summary dataframe for {stage_dir}")
    good_df = summary_df[summary_df["id"].astype(str).isin(report.good_ids)].copy()
    sort_columns = [col for col in ("e_ch", "e_hull/at", "id") if col in good_df.columns]
    if sort_columns:
        good_df = good_df.sort_values(sort_columns)
    summary_good_csv = output_dir / "summary_good.csv"
    table_good_txt = output_dir / "table_good.txt"
    good_df.to_csv(summary_good_csv, index=False)
    summary_tools.print_pretty_df(good_df, table_good_txt, sort_by=("e_ch" if "e_ch" in good_df.columns else None))
    written["summary_good_csv"] = summary_good_csv
    written["table_good_txt"] = table_good_txt
    if include_good_poscars:
        poscars_good_dir = output_dir / "POSCARS_good"
        poscars_src_dir = stage_dir / "POSCARS"
        poscars_good_dir.mkdir(exist_ok=True)
        for poscar_path in poscars_good_dir.iterdir():
            if poscar_path.is_file():
                poscar_path.unlink()
        for structure_id in good_df["id"].astype(str):
            src = poscars_src_dir / f"{structure_id}POSCAR"
            if not src.exists():
                warnings.warn(f"Missing POSCAR for good structure {structure_id} in {stage_dir}")
                continue
            shutil.copy2(src, poscars_good_dir / src.name)
        written["poscars_good_dir"] = poscars_good_dir
    if source_summary_path is not None:
        written["source_summary_csv"] = source_summary_path
    return written
