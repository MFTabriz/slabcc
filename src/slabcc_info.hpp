#pragma once

constexpr auto SLABCC_VERSION_MAJOR = 0;
constexpr auto SLABCC_VERSION_MINOR = 8;
constexpr auto SLABCC_VERSION_PATCH = 5;


#include "targetver.h"

#include <future>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>

#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>

#include <algorithm>

#ifdef MKL
#include "mkl_service.h"
#include "mkl_types.h"
#endif