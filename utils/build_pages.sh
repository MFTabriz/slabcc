#! /bin/bash

set -e

color_red="\033[31;49;1m"
#color_green="\033[32;49;1m"
color_reset="\033[0m"

if [ -z "$STAGE_PATH" ]; then
    printf "%bERROR: \$STAGE_PATH is not set! Skipping the deployment...%b\n" "${color_red}" "${color_reset}"
    exit 1
fi

mkdir -p "$STAGE_PATH"
INDEX_HTML=docs/manual.rst
SOURCE_LINK="https://codeberg.org/meisam/slabcc/src/commit/${CI_COMMIT_SHA}/${INDEX_HTML}"

echo "SOURCE LINK: $SOURCE_LINK"
tail -n +2 < "$INDEX_HTML" | rst2html.py --title="SLABCC" -g --source-url="$SOURCE_LINK" > "$STAGE_PATH"/index.html
