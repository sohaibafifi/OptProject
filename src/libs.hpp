#pragma once
#include <fstream>
#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <limits>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#if defined(WIN32) || defined(_WIN32)
#include <sys/sysinfo.h>
#endif
#ifdef MIPSOLVER
#pragma GCC system_header
#include <ilcplex/ilocplex.h>
#endif
#ifdef GMIPSOLVER
#pragma GCC system_header
#include "gurobi_c++.h"
#endif
#define PI 3.14159265
using namespace std;
#define IN
#define OUT
