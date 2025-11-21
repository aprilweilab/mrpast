# Releasing new mrpast versions

This outlines the steps for releasing a new version of mrpast. This is just for the wrapper
that posts mrpast to PyPi, there are other ways to install (build manually, or use the
Docker image).

## Version numbering

Bump the version in setup.py whenever a release is done on GitHub. Major version number is
incremented on breaking changes, minor version number is incremented for backwards-
compatible changes.

Version number needs to be updated in `setup.py` and `conda/meta.yaml`.

## Packaging for PyPi

Build the package distributions for PyPi. We build a source dist and then Linux binary distributions for recent Python versions. The container is from the [manylinux](https://github.com/pypa/manylinux) project.

```
# Remove dist/ to start fresh
rm -rf dist

# Pull the container for manylinux:
docker pull quay.io/pypa/manylinux_2_28_x86_64

# Run the packaging inside the container
docker run -v $PWD/../:/io -it quay.io/pypa/manylinux_2_28_x86_64 /io/mrpast/package.sh

# Fix file permissions from Docker
sudo chown -R ddehaas dist/
sudo chgrp -R ddehaas dist/

# Copy the source wheel to wheelhouse
cp ./dist/*.tar.gz ./dist/wheelhouse/

```

Test the source distribution:
```
pip install --force-reinstall ./dist/wheelhouse/mrpast-*.tar.gz
```


To upload to PyPi:
```
python3 -m twine upload dist/wheelhouse/*
```
