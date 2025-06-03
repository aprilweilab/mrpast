import os
import subprocess
import sys
import shutil
from pprint import pprint
from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
BUILD_DIR = "py_build"

env_debug = int(os.environ.get("MRPAST_DEBUG", 0))
env_native = int(os.environ.get("MRPAST_ENABLE_NATIVE", 0))

# This is for handling `python setup.py bdist_wheel`, etc.
extra_cmake_args = []
build_type = "Release"

if env_debug:
    build_type = "Debug"
if env_native:
    extra_cmake_args.append("-DENABLE_NATIVE=ON")


# This is for handling `pip install -e --install-option="..."`, etc.
class CommandBase:
    def initialize_options(self):
        super().initialize_options()

    def finalize_options(self):
        super().finalize_options()

    def run(self):
        build_solver()
        super().run()


class InstallCommand(CommandBase, install):
    pass


class DevelopCommand(CommandBase, develop):
    pass


def build_solver():
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
    if not os.path.exists(BUILD_DIR):
        os.makedirs(BUILD_DIR)

    print(f"Building with args: {cmake_args}")

    # Config and build the extension
    subprocess.check_call(
        ["cmake", THIS_DIR] + cmake_args, cwd=BUILD_DIR, stdout=sys.stdout
    )
    subprocess.check_call(
        ["cmake", "--build", ".", "--config", build_type, "--", "-j"],
        cwd=BUILD_DIR,
        stdout=sys.stdout,
    )
    shutil.copy(f"{BUILD_DIR}/{SOLVER_NAME}", f"{PACKAGE_NAME}/{SOLVER_NAME}")
    shutil.copy(f"{BUILD_DIR}/{EVAL_NAME}", f"{PACKAGE_NAME}/{EVAL_NAME}")


with open(os.path.join(THIS_DIR, "requirements.txt")) as f:
    requires = list(map(str.strip, f))
with open(os.path.join(THIS_DIR, "README.md")) as f:
    long_description = f.read()

PACKAGE_NAME = "mrpast"
SOLVER_NAME = "mrp-solver"
EVAL_NAME = "mrp-eval"
version = "0.1"
setup(
    cmdclass={
        "install": InstallCommand,
        "develop": DevelopCommand,
    },
    name=PACKAGE_NAME,
    packages=[PACKAGE_NAME],
    version=version,
    description="Demographic inference from Ancestral Recombination Graphs.",
    author="Drew DeHaas, April Wei",
    author_email="",
    url="https://aprilweilab.github.io/",
    zip_safe=False,
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
    package_data={
        PACKAGE_NAME: [SOLVER_NAME, EVAL_NAME],
    },
    include_package_data=True,
    entry_points={
        "console_scripts": ["mrpast=mrpast.main:main"],
    },
    install_requires=requires,
    long_description=long_description,
    long_description_content_type="text/markdown",
)
