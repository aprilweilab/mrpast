from mrpast.arginfer import get_vcf_stats
import unittest
import tempfile
import os


class ArgInferTests(unittest.TestCase):
    def test_vcf_stats(self):
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            fp = os.path.join(tmp_dir_name, "tmp.vcf")
            with open(fp, "w") as f:
                f.write(
                    """# Ignore
# Ignored because of leading #
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1
1\t12345\tv1\tA\tG\t1\tNONE\tBLAH\tGT\t0|0\t1|1
"""
                )
            (first, last), sites, individuals = get_vcf_stats(fp)
        self.assertEqual(first, last)
        self.assertEqual(sites, 1)
        self.assertEqual(individuals, 2)


if __name__ == "__main__":
    unittest.main()
