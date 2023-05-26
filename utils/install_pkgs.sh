#! /bin/bash

# reads $TOOLCHAIN from environment, installs packages and exports following variables to .env file:
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

install_numdiff () {
    curl -o numdiff-5.9.0.tar.gz https://de.freedif.org/savannah/numdiff/numdiff-5.9.0.tar.gz
    tar -xvf numdiff-5.9.0.tar.gz
    cd numdiff-5.9.0
    mkdir -p bin && cd bin
    ../configure
    make
    cd ../..
}
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
_CXX=${_cxx_pkg_compiler#*:}

echo "TOOLCHAIN            : $TOOLCHAIN"
echo "Container            : $_container"
echo "Packages             : ${_pkgs_array[*]}"
echo "C compiler package   : $_c_pkg"
echo "C++ compiler package : $_cxx_pkg"
echo "CC                   : $_CC"
echo "CXX                  : $_CXX"
echo "MKL                  : $_mkl"

_distro=${_container%:*}
if [[ "$_distro" == 'ubuntu' ]]; then 
    apt update
    # shellcheck disable=SC2068
    apt install -y make numdiff ${_pkgs_array[@]}

elif [[ "$_distro" == 'almalinux' ]]; then
    yum install -y dnf dnf-plugins-core && dnf config-manager --set-enabled powertools
    # shellcheck disable=SC2068
    dnf install -y make ${_pkgs_array[@]}
    [[ -x "numdiff" ]] || install_numdiff
    
elif [[ "$_distro" == 'opensuse/leap' ]]; then

    # shellcheck disable=SC2068
    zypper install -y make gzip ${_pkgs_array[@]}
    zypper addrepo https://download.opensuse.org/repositories/Base:System/standard/Base:System.repo
    zypper --gpg-auto-import-keys ref
    zypper install -y tar 
    [[ -x "numdiff" ]] || install_numdiff
elif [[ "$_distro" == 'intel/oneapi-basekit' ]]; then
    echo "Using MKL from OneAPI basekit..."
    apt update && apt install -y numdiff
else
    echo "ERROR: unsupported environment: $_distro"
    exit 1
fi

echo "export CC=$_CC && export CXX=$_CXX && export MKL=$_mkl && export PATH=$(pwd)/numdiff-5.9.0/bin:$PATH" > .env && chmod +x .env

$_CC --version | head -n 1
$_CXX --version | head -n 1