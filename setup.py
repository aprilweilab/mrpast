import os
import subprocess
import sys
import shutil
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

PACKAGE_NAME = "mrpast"
VERSION = "0.1"


THIS_DIR = os.path.dirname(os.path.realpath(__file__))

env_debug = int(os.environ.get("MRPAST_DEBUG", 0))
env_native = int(os.environ.get("MRPAST_ENABLE_NATIVE", 0))

# This is for handling `python setup.py bdist_wheel`, etc.
extra_cmake_args = []
build_type = "Release"

if env_debug:
    build_type = "Debug"
if env_native:
    extra_cmake_args.append("-DENABLE_NATIVE=ON")


class CMakeExtension(Extension):
    def __init__(
        self, name, cmake_lists_dir=".", sources=[], extra_executables=[], **kwa
    ):
        Extension.__init__(self, name, sources=sources, **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)
        self.extra_executables = extra_executables


class CMakeBuild(build_ext):
    def get_source_files(self):
        return [
            "CMakeLists.txt",
            "src/args.hxx",
            "src/common.h",
            "src/derivatives.h",
            "src/mrp_eval.cpp",
            "src/mrp_solver.cpp",
            "src/objective.cpp",
            "src/objective.h",
            "src/solve.cpp",
            "src/solve.h",
        ]

    def build_extensions(self):
        assert len(self.extensions) == 1
        ext = self.extensions[0]
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        try:
            subprocess.check_call(["cmake", "--version"])
        except OSError:
            raise RuntimeError("Cannot find CMake executable")

        cmake_args = [
            "-DCMAKE_BUILD_TYPE=%s" % build_type,
            "-DENABLE_CHECKERS=OFF",
            "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{}={}".format(build_type.upper(), "."),
            "-DNLOPT_PYTHON=OFF",
            "-DNLOPT_OCTAVE=OFF",
            "-DNLOPT_MATLAB=OFF",
            "-DNLOPT_GUILE=OFF",
            "-DNLOPT_SWIG=OFF",
        ] + extra_cmake_args
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        print(f"Building with args: {cmake_args}")

        # Config and build the executable
        subprocess.check_call(
            ["cmake", ext.cmake_lists_dir] + cmake_args,
            cwd=self.build_temp,
            stdout=sys.stdout,
        )
        subprocess.check_call(
            ["cmake", "--build", ".", "--config", build_type, "--", "-j"],
            cwd=self.build_temp,
            stdout=sys.stdout,
        )
        for executable in ext.extra_executables:
            shutil.copy2(
                os.path.join(self.build_temp, executable),
                os.path.join(extdir, executable),
            )


with open(os.path.join(THIS_DIR, "requirements.txt")) as f:
    requires = list(map(str.strip, f))
with open(os.path.join(THIS_DIR, "README.md")) as f:
    long_description = f.read()

SOLVER_NAME = "mrp-solver"
EVAL_NAME = "mrp-eval"

setup(
    name=PACKAGE_NAME,
    packages=[PACKAGE_NAME],
    version=VERSION,
    description="Demographic inference from Ancestral Recombination Graphs.",
    author="Drew DeHaas, April Wei",
    author_email="",
    url="https://aprilweilab.github.io/",
    zip_safe=False,
    ext_modules=[CMakeExtension(".", extra_executables=[SOLVER_NAME, EVAL_NAME])],
    cmdclass={"build_ext": CMakeBuild},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
    ],
    entry_points={
        "console_scripts": ["mrpast=mrpast.main:main"],
    },
    install_requires=requires,
    long_description=long_description,
    long_description_content_type="text/markdown",
)
