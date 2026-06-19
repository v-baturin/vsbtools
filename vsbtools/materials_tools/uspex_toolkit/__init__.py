from ..external_paths import add_sys_path, import_from_path_validator, resolve_external_path

USPEX_PYTHON_PATH = resolve_external_path(
    name="USPEX Python package",
    config_key="uspex_python_path",
    env_var="USPEX_PYTHON_PATH",
    validator=import_from_path_validator("USPEX.components"),
    prompt=False,
    required=False,
    prompt_text="Enter path to USPEX Python root: ",
)
if USPEX_PYTHON_PATH is not None:
    add_sys_path(USPEX_PYTHON_PATH)
