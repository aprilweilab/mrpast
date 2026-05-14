import sys
import json
import os

filename = os.path.join("aa5", "aa5.bestlhoods")
with open(filename) as f:
    columns = None
    values = None
    for i, line in enumerate(f):
        if not line.strip():
            continue
        if i == 0:
            columns = line.strip().split("\t")
        else:
            assert i == 1
            values = line.strip().split("\t")
    assert len(columns) == len(values)
    del values[columns.index("MaxObsLhood")]
    del values[columns.index("MaxEstLhood")]
    del columns[columns.index("MaxObsLhood")]
    del columns[columns.index("MaxEstLhood")]
    assert len(columns) == len(values)
    print(json.dumps({c: float(v) for c, v in zip(columns, values)}, indent=2))
