"""
Some very simple tests for end-to-end functionality.
"""

import os
import subprocess
import tempfile
import unittest
from typing import List, Dict

THISDIR = os.path.dirname(os.path.realpath(__file__))

MODEL_5D1E = os.path.join(THISDIR, "..", "examples", "5deme1epoch.yaml")
MODEL_OOA3 = os.path.join(THISDIR, "..", "examples", "ooa_3g09.yaml")
MODEL_OOA2 = os.path.join(THISDIR, "..", "examples", "ooa_2t12.yaml")
MODEL_OOA4 = os.path.join(THISDIR, "..", "examples", "ooa_4j17.yaml")

try:
    MRPAST = [subprocess.check_output(["which", "mrpast"]).decode("utf-8").strip()]
except subprocess.CalledProcessError:
    MRPAST = ["python", "mrpast/main.py"]


class EndToEndTests(unittest.TestCase):
    def simulate(self, model_file: str, tmpdirname: str):
        arg_prefix = os.path.join(tmpdirname, "test_arg")
        command = list(
            map(
                str,
                MRPAST
                + [
                    "simulate",
                    "--replicates",
                    1,
                    "--seq-len",
                    10_000_000,
                    "--individuals",
                    5,
                    "--debug-demo",
                    model_file,
                    arg_prefix,
                ],
            )
        )
        subprocess.check_call(command)
        outfile = arg_prefix + "_0-0.trees"
        self.assertTrue(os.path.isfile(outfile))
        return outfile

    def process(
        self,
        model_file: str,
        tmpdirname: str,
        leave_out: List[int] = [],
        pop_map: Dict[int, int] = {},
    ):
        arg_prefix = os.path.join(tmpdirname, "test_arg")
        extra_args = []
        if leave_out:
            extra_args.extend(["--leave-out", ",".join(map(str, leave_out))])
        if pop_map:
            extra_args.extend(
                ["--map-pops", ",".join([f"{k}:{v}" for k, v in pop_map.items()])]
            )
        command = list(
            map(
                str,
                MRPAST
                + [
                    "process",
                    "--out-dir",
                    tmpdirname,
                ]
                + extra_args
                + [
                    model_file,
                    arg_prefix,
                ],
            )
        )
        subprocess.check_call(command)

    def test_demes_5d1e(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.simulate(MODEL_5D1E, tmpdirname)

    def test_demes_ooa3(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.simulate(MODEL_OOA3, tmpdirname)

            # Now we process this simulated data three separate ways:
            # 1. With the same model that generated it.
            self.process(MODEL_OOA3, tmpdirname)
            # 2. With the same model that generated it, but treat pop2 as unsampled.
            self.process(MODEL_OOA3, tmpdirname, leave_out=[2])
            # 3. With a model that contains _fewer_ populations, here we have to drop
            #    one of the populations and make sure the mapping is correct.
            self.process(MODEL_OOA2, tmpdirname, leave_out=[2])
            # 4. With a model that contains _more_ populations, here we have to specify
            #    the population mapping, which can include mapping a non-existing ARG
            #    population to the last population (which treats it as unsampled).
            self.process(MODEL_OOA4, tmpdirname, pop_map={0: 0, 1: 1, 2: 2, 3: 3})


if __name__ == "__main__":
    unittest.main()
