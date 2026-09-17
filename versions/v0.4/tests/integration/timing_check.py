# SPDX-License-Identifier: Apache-2.0
"""Check reported containment against actual accepted-step wall time."""
import json
import math


def check_timing(path, active):
    rows = [json.loads(line)['payload'] for line in path.read_text().splitlines()]
    assert rows
    for row in rows:
        wall = row['seconds']
        assert math.isfinite(wall) and wall > 0
        physics = row['physics_modules']
        assert physics['scope'] == 'max_rank_all_attempts_inclusive_modules'
        assert physics['containment'] == 'inside_advance_overlaps_cn_phases'
        values = physics['seconds']
        assert set(values) == {'reaction_sources', 'mean_reaction', 'esf_reaction', 'tcr_statistics'}
        assert all(math.isfinite(v) and 0 <= v <= wall for v in values.values()), row
        for name in values:
            assert (values[name] > 0) == (name in active), (name, values)
        comm = row['communication_observations']
        assert comm['scope'] == 'max_rank_all_attempts_named_operations'
        assert comm['containment'] == 'inside_advance_overlaps_modules_and_cn_phases'
        assert set(comm['seconds']) == {'structured_wait', 'structured_control', 'linear_reductions'}
        assert all(math.isfinite(v) and 0 <= v <= wall for v in comm['seconds'].values()), row
    return rows
