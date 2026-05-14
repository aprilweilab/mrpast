import pandas as pd
import os
import sys

THRESHOLD = 0.8

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            "Usage: bscore_modify.py <BED file> <rate map> <chromosome>",
            file=sys.stderr,
        )
        exit(1)

    bed_file = sys.argv[1]
    rate_map = sys.argv[2]
    chrom = int(sys.argv[3])

    b_scores = pd.read_csv(bed_file, sep="\t")
    b_scores = b_scores[b_scores["chrom"] == f"chr{chrom}"]
    print(f"Saw {len(b_scores)} scores", file=sys.stderr)

    b_scores_belowT = b_scores[b_scores["score"] <= THRESHOLD]
    print(f"Saw {len(b_scores_belowT)} below threshold", file=sys.stderr)

    # This is pretty inefficient: mark all BP in the regions w/ selection
    inflate_position = set()
    for i, row in b_scores_belowT.iterrows():
        for j in range(row["start"], row["end"]):
            inflate_position.add(j)

    # Now process the rate map
    with open(rate_map) as f:
        for line in f:
            row = line.strip().split(" ")
            start = int(float(row[0]))
            end = int(float(row[1]))
            if start in inflate_position or end in inflate_position:
                rate = 1e-4
            else:
                rate = float(row[2])
            print(" ".join([row[0], row[1], str(rate)]))
