.. _installation:

Installing mrpast
=================

The simplest way to install mrpast is via pip:

::

	pip install mrpast

The resulting installation has two parts:

  1. The ``mrpast`` command, used for the majority of mrpast's functionality. See ``mrpast --help``.
  2. The Python API, accessed via ``import mrpast``. This API is primarily used for plotting and
     interpreting results from mrpast - see the `Python API reference <python_api.html>`_ for more
     details.

On MacOS platforms or non-standard Linux platforms, the installation via ``pip`` may actually compile
mrpast from source. In this case, you'll need `CMake <https://cmake.org>`_ at least version 3.10, and
a g++ or Clang version that supports C++17.

It is best practice to use a `virtual environment <https://docs.python.org/3/library/venv.html>`_
when installing mrpast, since a lot of additional packages may be installed as well (numpy, tskit, etc.).

All default dependencies of mrpast will be installed automatically, but there is some optional functionality
that you may need to install extra dependencies for:

- The documentation section on `Relate <relate.html>`_ describes relevant installation steps.
- The documentation section on `tsinfer <tsinfer.html>`_ describes relevant installation steps.
- The documentation section on `SINGER <singer.html>`_ describes relevant installation steps.
- Viewing of mrpast models and conversion to/from `Demes format <https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html>`_ requires ``pip install demes networkx matplotlib``
- Some of the examples in the documentation make use of `demesdraw <https://github.com/grahamgower/demesdraw>`_ which can be obtained via ``pip install demesdraw``

Advanced: Installing from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only developers should need to install from source.  Installing from source:

::

  git clone --recursive https://github.com/aprilweilab/mrpast.git
  cd mrpast
  # Compiles C++, copies the binaries to the right directories, and install the Python code in
  # an editable way (so changes are reflected immediately)
  MRPAST_ENABLE_NATIVE=1 pip install -v -e .

There are other environment variables that control the behavior of the mrpast build as well:

- ``MRPAST_DEBUG=1``: Build the C++ code in debug mode.