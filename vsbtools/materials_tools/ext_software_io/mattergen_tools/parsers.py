from pathlib import Path
from typing import List, Any
import yaml, re

STANDARD_PARAM_DICT_KEYS = ['chemical_system', 'guidance', 'diffusion_loss_weight', 'algo']

class SafeLoaderWithNone(yaml.SafeLoader): pass
SafeLoaderWithNone.add_implicit_resolver(
    "tag:yaml.org,2002:null",
    re.compile(r"^(?:~|null|Null|NULL|none|None)$"),
    list("~nN"),
)

_MULTI_US = re.compile(r"_+")

def _fmt_num(x: Any) -> str:
    if isinstance(x, float) and x.is_integer():
        return str(int(x))
    return str(x)

def _sanitize_token(s: str) -> str:
    # Turn any run of braces/brackets/parentheses/quotes/spaces into "_"
    s = re.sub(r"[{}\[\]\(\)\"'\s]+", "_", s.strip())
    s = _MULTI_US.sub("_", s)
    return s

def serialize_structure(obj: Any) -> str:
    """
    Dicts: key and value concatenated with '_'
    Lists/Tuples: elements joined with '-'
    Scalars: sanitized; hyphens and dots are preserved
    Leading/trailing '_' are stripped at the end.
    """
    def ser(x: Any) -> str:
        if isinstance(x, dict):
            parts: list[str] = []
            for k, v in x.items():  # preserves insertion order
                k_str = _sanitize_token(str(k))
                v_str = ser(v)
                parts.append(k_str if not v_str else f"{k_str}_{v_str}")
            return "_".join(p for p in parts if p)
        if isinstance(x, (list, tuple)):
            items = [ser(e) for e in x]
            items = [i for i in items if i]  # drop empties
            return "-".join(items)
        if x is None:
            return "None"
        if isinstance(x, bool):
            return "True" if x else "False"
        if isinstance(x, (int, float)):
            return _fmt_num(x)
        return _sanitize_token(str(x))

    s = ser(obj)
    s = _MULTI_US.sub("_", s).strip("_")
    return s


def input_parameters_to_dict(inp_par_file: Path | str | None = None, raw: str | None = None):
    if inp_par_file is not None:
        raw = Path(inp_par_file).read_text(encoding="utf-8")
    data = yaml.load(raw, Loader=SafeLoaderWithNone)
    if not isinstance(data, dict):
        raise ValueError("Parsed content is not a mapping/dict")
    return data

def fname_friendly_serialize(d: dict, dict_keys: List|Any = None):

    if dict_keys is None: dict_keys = STANDARD_PARAM_DICT_KEYS
    parts = []
    for k in dict_keys:
        if k == 'chemical_system':
            serialized_val = serialize_structure(sorted(d['properties_to_condition_on'][k].split('-')))
        elif k == 'algo':
            val = d.get(k, None)
            serialized_val = str(int(val)) if val is not None else 'None'
        else:
            serialized_val = serialize_structure(d.get(k, None))
        if k and k != 'chemical_system':
            parts.append(f"{k}_{serialized_val}")
        else:
            parts.append(serialized_val)
    return '__'.join(parts)

def input_paramfile_2_dirname(input_param_file, dict_keys=None):
    if dict_keys is None:
        dict_keys = STANDARD_PARAM_DICT_KEYS
    return fname_friendly_serialize(input_parameters_to_dict(inp_par_file=input_param_file), dict_keys=dict_keys)

if __name__ == '__main__':
    inpar_path_guided = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/INBOX/Si-O_guided_20250929_0_SiO6/Si-O_guided_0/input_parameters.txt")
    inpar_path_nonguided = Path("/home/vsbat/work/generations/Si-O_no_guidance/input_parameters.txt")
    param_dict_guided = input_parameters_to_dict(inpar_path_guided)
    param_dict_nonguided = input_parameters_to_dict(inpar_path_nonguided)
    print(param_dict_guided)
    fname = fname_friendly_serialize(param_dict_guided, ['chemical_system', 'guidance', 'diffusion_loss_weight', 'algo'])
    print(f"filename_guided = {fname}")

    print(param_dict_nonguided)
    fname = fname_friendly_serialize(param_dict_nonguided, ['chemical_system', 'guidance', 'diffusion_loss_weight', 'algo'])
    print(f"filename_nonguided = {fname}")

    from vsbtools.materials_tools.materials_dataset.io.yaml_csv_poscars import read
    ds = read("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED/Cu-Si-P-Ca/"
              "Cu-Si-P-Ca__guidance_environment_mode_huber_Cu-P_4-2.2_Cu-Cu_0-2.9__diffusion_loss_weight_None__algo_None_annotated/"
              "0_11dec62bb404005d/manifest.yaml")
    print(ds.metadata["batch_metadata"])
    param_dict_guided = input_parameters_to_dict(raw=ds.metadata["batch_metadata"])
    fname = fname_friendly_serialize(param_dict_guided, ['chemical_system', 'guidance'])
    print(f"filename_from_ds_metadata = {fname}")