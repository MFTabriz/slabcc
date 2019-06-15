#pragma once

constexpr auto SLABCC_VERSION_MAJOR = 0;
constexpr auto SLABCC_VERSION_MINOR = 8;
constexpr auto SLABCC_VERSION_PATCH = 1;

//comment the next line while developing the code! Uncommenting will disable the range checks in Armadillo. Makefile may override this!
//#define ARMA_NO_DEBUG


#include "targetver.h"

#include <future>

#include <stdio.h>
#include <iomanip> 
#include <iostream> 
#include <fstream> 
#include <sstream>

#include <chrono>
#include <vector>  
#include <unordered_map>
#include <string>  

#include <algorithm> 


#ifdef MKL
#include "mkl_service.h"
#include "mkl_types.h"
#endif