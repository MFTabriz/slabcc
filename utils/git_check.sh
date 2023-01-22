#! /bin/bash

# This script runs "git add" on all of the files in current git repository.
# If any of these files are filtered by .gitignore, it will return error.

color_red="\033[31;49;1m"
color_green="\033[32;49;1m"
color_reset="\033[0m"

errfile=$(mktemp)
git_repo_err=0

pushd "$(git rev-parse --show-toplevel)" > /dev/null || exit 1
git config --local advice.addIgnoredFile false

files_list=$(find . -not -path './.git/*' -type f)
IFS=$'\n'
# shellcheck disable=SC2068
for file in $files_list; do
    git add "$file" 2> "$errfile"
    done=$?

    if [ "$done" -ne 0 ]; then
        printf "Checking: %s \n%bERROR: " "$file" "${color_red}"
        cat "$errfile"
        printf "%b\n\n" "${color_reset}"
        git_repo_err=1
    fi
done

if [ "$git_repo_err" -eq 0 ]; then
    printf "%bPASS: no file in current git repository is ignored by your .gitignore files.%b\n" "${color_green}" "${color_reset}"
fi

popd > /dev/null || exit

exit "$git_repo_err"