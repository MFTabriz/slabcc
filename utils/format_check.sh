#! /bin/bash

set -e

color_red="\033[31;49;1m"
color_green="\033[32;49;1m"
color_reset="\033[0m"

function finish() {
    echo "Cleaning up..."
	mapfile -t CHILDREN_PID_LIST < <(ps --ppid $$ -o pid= || true)
        for child_pid in "${CHILDREN_PID_LIST[@]}"; do
            kill "$child_pid" &>/dev/null || true
            wait "$child_pid" &>/dev/null || true
        done
	rm -fr "$TMP_DIR"
}
trap finish EXIT


if [[ -z "$CLANG_FORMAT_EXEC" ]]; then CLANG_FORMAT_EXEC='clang-format'; fi;

if [[ -x "$(command -v "$CLANG_FORMAT_EXEC")" ]]; then
	"$CLANG_FORMAT_EXEC" --version
else
	echo -e "${color_red}ERROR: $CLANG_FORMAT_EXEC was not found!${color_reset}"
	exit 1
fi

TMP_DIR=$(mktemp -d)

cpp_sources="$(find src/* -maxdepth 0 -name '*.cpp' -or -name '*.hpp')"
for source_file in ${cpp_sources}; do
	styled_source_file="$TMP_DIR/$source_file"
	mkdir -p "${styled_source_file%/*}"
	"$CLANG_FORMAT_EXEC" --verbose -style=file "$source_file" > "$styled_source_file" || fail=1
	if [[ -s "$styled_source_file" ]]; then
		diff=$(diff -q "$source_file" "$styled_source_file" || true)
		if [[ -n "$diff" ]]; then
			printf "\n%bERROR: File is not formatted correctly: %s%b\n\n" "${color_red}" "$source_file" "${color_reset}"
			diff "$source_file" "$styled_source_file" || true
			fail=1
		fi
	else
		printf "\n%bERROR: clang-format failed to generate output file!%b\n" "${color_red}" "${color_reset}"
		fail=1
	fi
done
	

if [[ "$fail" == "1" ]]; then
	echo "CLANG_FORMAT_EXEC: $CLANG_FORMAT_EXEC"
	echo "default llvm config:" && "$CLANG_FORMAT_EXEC" -style=llvm -dump-config
	exit 1
else
	printf "\n%bPASS: all files are styled correctly!%b\n" "${color_green}" "${color_reset}"
fi
