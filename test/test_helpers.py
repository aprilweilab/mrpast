from mrpast.helpers import remove_ext, count_lines
import unittest
import tempfile
import os


class HelperTests(unittest.TestCase):
    def test_remove_ext(self):
        result = remove_ext("testing.1234.dfsd.foobar.foobar", "foobar")
        self.assertEqual(result, "testing.1234.dfsd.foobar")

    def test_count_lines(self):
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            fn = os.path.join(tmp_dir_name, "tmp.txt")
            with open(fn, "w") as f:
                for i in range(100):
                    f.write("test\n")
            result = count_lines(fn)
            self.assertEqual(result, 100)


if __name__ == "__main__":
    unittest.main()
