import tskit
import sys
import math

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            "Convert a hapmap-style recombination map to an msprime/SINGER style rate map"
        )
        print("Usage: make_rate_map.py <hapmap file>")
        exit(1)

    map = tskit.RateMap.read_hapmap(sys.argv[1], rate_col=2)
    for left, right, rate in zip(map.left, map.right, map.rate):
        assert not math.isnan(left)
        assert not math.isnan(right)
        if math.isnan(rate):
            continue
        print(f"{left} {right} {rate}")
