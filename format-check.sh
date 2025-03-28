#!/bin/bash

set -ev

clang-format --dry-run -Werror src/*.cpp
clang-format --dry-run -Werror src/*.h
