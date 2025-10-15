#!/bin/bash
set -ev

./format-check.sh

flake8 mrpast --count --select=E9,F63,F7,F82,F401 --show-source --statistics

black mrpast/ setup.py test/ scripts/ --exclude "third-party" --check

mypy mrpast --no-namespace-packages --ignore-missing-imports

pytest test/

./build/mrpast_test
