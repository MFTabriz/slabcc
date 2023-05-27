#! /bin/bash
# this script generates single source rst files. run this script in the root of repository to generate single source files.
## $> utils/generate_docs.sh

set -e

rst_include include .readme.rst.template readme.rst
rst_include include docs/.manual.rst.template docs/manual.rst