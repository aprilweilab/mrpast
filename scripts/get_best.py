# Given a list of JSON filenames (output from the mrpast solver), find the one with
# the highest likelihood (lowest negative log-likehood).
#   Example usage: python get_best.py *.json
import sys
import json
from typing import List, Optional


def get_best(filenames: List[str]) -> Optional[str]:
    bestLL = 2**64
    best = None
    for fn in filenames:
        with open(fn) as f:
            data = json.load(f)
        negLL = data["negLL"]
        if negLL is None:
            negLL = 2**64
        if negLL < bestLL:
            bestLL = negLL
            best = fn
    return best


if __name__ == "__main__":
    best = get_best(sys.argv[1:])
    print(best)
