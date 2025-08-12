"""
Some very simple tests for end-to-end functionality.
"""

import os
import subprocess
import tempfile
import unittest
import glob

THISDIR = os.path.dirname(os.path.realpath(__file__))

MODEL_5D1E = os.path.join(THISDIR, "..", "examples", "5deme1epoch.yaml")
MODEL_OOA3 = os.path.join(THISDIR, "..", "examples", "ooa_3g09.yaml")

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

    def test_demes_5d1e(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.simulate(MODEL_5D1E, tmpdirname)

    def test_demes_ooa3(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.simulate(MODEL_OOA3, tmpdirname)


if __name__ == "__main__":
    unittest.main()
