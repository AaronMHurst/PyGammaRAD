#!/bin/bash

# Script to remove folders created during package-building process.
# Uninstalls PyGammaRAD package, then re-installs it.
# Run this script after editing source modules in the PyGammaRAD package.

rm -rf build
rm -rf dist
rm -rf PyGammaRAD.egg-info

rm -f *~
rm -f PyGammaRAD/*~
rm -f PyGammaRAD/data/*~
rm -f tests/*~
rm -rf PyGammaRAD/__pycache__
rm -rf __pycache__
rm -rf .tox

pip uninstall PyGammaRAD<<EOF
y
EOF

python setup.py install
