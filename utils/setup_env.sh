#! /bin/bash

# reads $TOOLCHAIN from environment, installs packages (if required) and exports following variables to .env file:
## CC
## CXX
## MKL

set -e


if [[ -z $TOOLCHAIN ]]; then
    # shellcheck disable=SC2016
    echo '
ERROR: $TOOLCHAIN is not set!

$TOOLCHAIN string syntax:
    "CONTAINER:tag C_compiler_package[:C_compiler_exec] C++_compiler_package[:C++_compiler_exec] [other packages]"
    (defining "C/C++_compiler_exec" can be omitted if compiler package name is same as the compiler executable) '
    exit 1
fi

_container="${TOOLCHAIN%% *}"
_pkgs_compilers_list="${TOOLCHAIN##"$_container" }"
read -ra _pkgs_compilers_array <<< "$_pkgs_compilers_list"

_c_pkg_compiler=${_pkgs_compilers_array[0]}
_cxx_pkg_compiler=${_pkgs_compilers_array[1]}

_mkl=0
if [[ "${_pkgs_compilers_array[2]}" == "intel-oneapi-mkl-devel" ]]; then
    _mkl=1
fi

for pkg in "${_pkgs_compilers_array[@]}"; do
    _pkg_name=${pkg%%:*}
    _pkgs_array=("${_pkgs_array[@]}" "$_pkg_name")
done

_c_pkg=${_c_pkg_compiler%:*}
_cxx_pkg=${_cxx_pkg_compiler%:*}
_CC=${_c_pkg_compiler#*:}
_CC_version=$($_CC --version | head -n 1)
_CXX=${_cxx_pkg_compiler#*:}
_CXX_version=$($_CXX --version | head -n 1)

echo "TOOLCHAIN            : $TOOLCHAIN"
echo "Container            : $_container"
echo "Packages             : ${_pkgs_array[*]}"
echo "C compiler package   : $_c_pkg"
echo "C++ compiler package : $_cxx_pkg"
echo "CC                   : $_CC ($_CC_version)"
echo "CXX                  : $_CXX ($_CXX_version)"
echo "MKL                  : $_mkl"

_distro=${_container%:*}
if [[ "$_distro" == 'intel/oneapi-basekit' ]]; then
    apt update && apt install -y numdiff
fi

echo "export CC=$_CC && export CXX=$_CXX && export MKL=$_mkl" > .env
chmod +x .env