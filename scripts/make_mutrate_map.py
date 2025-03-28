# Simple script to make a rate map that SINGER can understand, for mutation rate.
#   Example usage: python make_mutrate_map.py 1e-8
import msprime
import sys
import math

MAX_CHROM_LENGTH = 500_000_000

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: make_mutrate_map.py <rate>")
        exit(1)
    rate = float(sys.argv[1])
    map = msprime.RateMap(position=[0, MAX_CHROM_LENGTH], rate=[rate])
    for left, right, rate in zip(map.left, map.right, map.rate):
        assert not math.isnan(left)
        assert not math.isnan(right)
        if math.isnan(rate):
            continue
        print(f"{left} {right} {rate}")
