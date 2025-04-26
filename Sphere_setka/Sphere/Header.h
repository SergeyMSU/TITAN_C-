#pragma once

#define kv(x) ((x) * (x))

#include <algorithm>
#include <omp.h>
#include <mutex>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

const double pi = acos(-1);

double polar_angle(const double& x, const double& y);

#include "point.h"
#include "Cell.h"
#include "Setka.h"


