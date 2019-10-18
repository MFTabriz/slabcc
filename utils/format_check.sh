#! /bin/bash
files_to_check="$(find src/* -maxdepth 0 -name '*.cpp' -or -name '*.hpp')"
for f in ${files_to_check}; do
	d=$(diff -u "$f" <(clang-format -style=file "$f") || true)
	if ! [ -z "$d" ]; then
		printf ":\n%s\n" "$d"
		fail=1
	fi
done

if [ "$fail" = 1 ]; then
		exit 1
fi

exit 0