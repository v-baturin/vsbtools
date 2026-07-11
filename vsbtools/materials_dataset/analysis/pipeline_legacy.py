from enum import IntEnum
from typing import Dict

class LegacyStage(IntEnum):
    parse_raw = 0
    symmetrize = 1
    poll_db = 2
    augment_by_ref = 3
    estimate = 4
    filter_hull = 5
    deduplicate = 6
    postprocess_dft = 7

LEGACY_DICTIONARY: Dict["old_name", "new_name"] = {'augment_raw_by_db': "merge_base_into_ref",
                                                   "augment_by_ref": "merge_base_into_ref",
                                                   'symmetrize_raw': "symmetrize"}
LEGACY_INDEX_TO_NAME = {st.value: st.name for st in LegacyStage}
LEGACY_NAME_TO_INDEX = {st.name: st.value for st in LegacyStage}