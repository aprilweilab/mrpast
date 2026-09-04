import sys
import json
import os

# NOTE: the model outputs still contain N1, but just like mrpast and gLike
# that value is unused in downstream comparison, because for all methods
# we let them implicitly calculate it according to the exponential growth.
# I.e., N1 result is invalid and also ignored, downstream.

for rep in range(50):
    filename = os.path.join(f"rep{rep}", "1t12", "1t12.bestlhoods")
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
        print(json.dumps({c: float(v) for c, v in zip(columns, values)}))
