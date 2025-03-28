from mrpast.model import validate_model, SymbolicMatrices, load_model_config
from mrpast.simulate import build_demography
from mrpast.from_demes import convert_from_demes
from mrpast.helpers import dump_model_yaml
import demes
import unittest
import tempfile
import os

THISDIR = os.path.dirname(os.path.realpath(__file__))

MODEL_5D1E = os.path.join(THISDIR, "..", "examples", "5deme1epoch.yaml")


class ModelTests(unittest.TestCase):
    # This test round-trips one of our example models through the Demes to/from conversion
    # and then verifies that the resulting parameters are the same.
    def test_demes_integration(self):
        validate_model(MODEL_5D1E)
        demography, _ = build_demography(MODEL_5D1E)
        demes_model = demography.to_demes()
        with tempfile.TemporaryDirectory() as tmpdirname:
            demes_file = os.path.join(tmpdirname, "testing.demes.yaml")
            demes.dump(demes_model, demes_file)
            roundtrip_model = convert_from_demes(demes_file)

            rt_model_file = os.path.join(tmpdirname, "roundtrip.yaml")
            with open(rt_model_file, "w") as fout:
                dump_model_yaml(roundtrip_model, fout)
            validate_model(rt_model_file)

            config_orig = load_model_config(MODEL_5D1E)
            coal_orig = SymbolicMatrices.from_config(
                config_orig["coalescence"], is_vector=True
            )
            mig_orig = SymbolicMatrices.from_config(config_orig["migration"])

            config_copy = load_model_config(rt_model_file)
            coal_copy = SymbolicMatrices.from_config(
                config_copy["coalescence"], is_vector=True
            )
            mig_copy = SymbolicMatrices.from_config(config_copy["migration"])

            for m_orig, m_copy in zip(coal_orig.matrices, coal_copy.matrices):
                for idx_orig, idx_copy in zip(m_orig, m_copy):
                    assert idx_orig != 0 or (idx_copy == idx_orig)
                    if idx_orig == 0:
                        continue
                    orig = coal_orig.get_parameter(idx_orig)
                    copy = coal_copy.get_parameter(idx_copy)
                    self.assertAlmostEqual(orig.ground_truth, copy.ground_truth, 6)
                    self.assertAlmostEqual(orig.lower_bound, copy.lower_bound, 6)
                    self.assertAlmostEqual(orig.upper_bound, copy.upper_bound, 6)

            for m_orig, m_copy in zip(mig_orig.matrices, mig_copy.matrices):
                for row_orig, row_copy in zip(m_orig, m_copy):
                    for idx_orig, idx_copy in zip(row_orig, row_copy):
                        assert idx_orig != 0 or (idx_copy == idx_orig)
                        if idx_orig == 0:
                            continue
                        orig = mig_orig.get_parameter(idx_orig)
                        copy = mig_copy.get_parameter(idx_copy)
                        self.assertAlmostEqual(orig.ground_truth, copy.ground_truth, 6)
                        self.assertAlmostEqual(orig.lower_bound, copy.lower_bound, 6)
                        self.assertAlmostEqual(orig.upper_bound, copy.upper_bound, 6)


if __name__ == "__main__":
    unittest.main()
