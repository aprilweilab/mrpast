from mrpast.model import UserModel, ParamRef
from mrpast.simulate import build_demography
from mrpast.from_demes import convert_from_demes
from mrpast.helpers import dump_model_yaml
import demes
import unittest
import tempfile
import os

THISDIR = os.path.dirname(os.path.realpath(__file__))

MODEL_5D1E = os.path.join(THISDIR, "..", "examples", "5deme1epoch.yaml")
MODEL_6D2ES = os.path.join(THISDIR, "..", "examples", "6deme2epoch.split.yaml")
MODEL_OOA2 = os.path.join(THISDIR, "..", "examples", "ooa_2t12.yaml")
MODEL_OOA3 = os.path.join(THISDIR, "..", "examples", "ooa_3g09.yaml")


class ModelTests(unittest.TestCase):
    # This test round-trips one of our example models through the Demes to/from conversion
    # and then verifies that the resulting parameters are the same.
    def run_demes_integration(self, model_name: str):
        in_model = UserModel.from_file(model_name)
        demography, _ = build_demography(in_model)
        demes_model = demography.to_demes()
        with tempfile.TemporaryDirectory() as tmpdirname:
            demes_file = os.path.join(tmpdirname, "testing.demes.yaml")
            demes.dump(demes_model, demes_file)
            roundtrip_model = convert_from_demes(demes_file)

            rt_model_file = os.path.join(tmpdirname, "roundtrip.yaml")
            with open(rt_model_file, "w") as fout:
                dump_model_yaml(roundtrip_model, fout)
            UserModel.from_file(rt_model_file)

            config_orig = UserModel.from_file(model_name)
            config_copy = UserModel.from_file(rt_model_file)

            # We want deme names, _not_ indexes, so that we can compare. NOTE: this restricts
            # what we can do with these objects, since many methods require resolved names.
            config_orig.unresolve_names()
            config_copy.unresolve_names()

            self.assertEqual(
                len(config_orig.coalescence.entries),
                len(config_copy.coalescence.entries),
            )
            coals_orig = sorted(
                config_orig.coalescence.entries, key=lambda e: (e.epoch, e.deme)
            )
            coals_copy = sorted(
                config_copy.coalescence.entries, key=lambda e: (e.epoch, e.deme)
            )
            for e_orig, e_copy in zip(coals_orig, coals_copy):
                self.assertEqual(e_orig.epoch, e_copy.epoch)
                self.assertEqual(e_orig.deme, e_copy.deme)
                if isinstance(e_orig.rate, ParamRef):
                    self.assertTrue(isinstance(e_copy.rate, ParamRef))
                    orig = config_orig.coalescence.get_parameter(e_orig.rate.param)
                    copy = config_copy.coalescence.get_parameter(e_copy.rate.param)
                    self.assertAlmostEqual(orig.ground_truth, copy.ground_truth, 6)
                    # The bounds don't exist in the demes model, so we can't test this.
                    # self.assertAlmostEqual(orig.lb, copy.lb, 6)
                    # self.assertAlmostEqual(orig.ub, copy.ub, 6)

                self.assertEqual(
                    len(config_orig.migration.entries),
                    len(config_copy.migration.entries),
                )
            migs_orig = sorted(
                config_orig.migration.entries, key=lambda e: (e.epoch, e.source, e.dest)
            )
            migs_copy = sorted(
                config_copy.migration.entries, key=lambda e: (e.epoch, e.source, e.dest)
            )
            for e_orig, e_copy in zip(migs_orig, migs_copy):
                self.assertEqual(e_orig.epoch, e_copy.epoch)
                self.assertEqual(e_orig.source, e_copy.source)
                self.assertEqual(e_orig.dest, e_copy.dest)
                if isinstance(e_orig.rate, ParamRef):
                    self.assertTrue(isinstance(e_copy.rate, ParamRef))
                    orig = config_orig.migration.get_parameter(e_orig.rate.param)
                    copy = config_copy.migration.get_parameter(e_copy.rate.param)
                    self.assertAlmostEqual(orig.ground_truth, copy.ground_truth, 6)
                    # The bounds don't exist in the demes model, so we can't test this.
                    # self.assertAlmostEqual(orig.lb, copy.lb, 6)
                    # self.assertAlmostEqual(orig.ub, copy.ub, 6)

    def test_demes_5d1e(self):
        self.run_demes_integration(MODEL_5D1E)

    def test_demes_ooa2(self):
        self.run_demes_integration(MODEL_OOA2)

    def test_demes_ooa3(self):
        self.run_demes_integration(MODEL_OOA3)

    def test_demes_6d2es(self):
        self.run_demes_integration(MODEL_6D2ES)


if __name__ == "__main__":
    unittest.main()
