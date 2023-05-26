#! /bin/bash
# this script generates single source rst files

set -e

rst_include include .readme.rst.template readme.rst
rst_include include docs/.manual.rst.template docs/manual.rst