#pragma once


// Описание всех классов
class Geo_param; 
class Phys_param;
class Setka;
class Yzel;
class Luch;
class Cell;
class Gran;
class Surfaces;


#define kv(x) ((x) * (x))
#define kyb(x) ((x) * (x) * (x))
#define kvv(x, y, z) ((x) * (x) + (y) * (y) + (z) * (z))
#define norm2(x, y, z) (sqrt(kvv(x, y, z)))
#define whach(x) cout << #x <<": " << (x) << endl


#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <map>
#include <math.h>
#include <cmath>
#include <limits>
#include <iterator>
#include <cstdlib>
#include <Eigen/Dense>
#include <optional>
#include <omp.h>
#include <chrono>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>




const double const_pi = 3.14159265358979323846;
const double cpi4 = 4.0 * const_pi;
const double cpi8 = 8.0 * const_pi;
const double spi4 = sqrt(cpi4);
const double eps = 1E-12;
const double epsb = 1E-4;
const double eps_p = 1E-6;
const double eps_d = 1E-3;
const double MF_meDmp = (1.0 / 1836.15);  // Отношения массы электрона к массе протона

// Boost библиотека (надо подключать к компилятору отдельно)
#include "boost/multi_array.hpp"
#include <boost/parameter.hpp>

using namespace std;

#include "Phys_param.h"
#include "Geo_param.h"
#include "Setka.h"
#include "Help_funk.h"
#include "Luch.h"
#include "Gran.h"
#include "Cell.h"
#include "Yzel.h"
#include "Surfaces.h"






// Описание всех типов комметариев

// &INIT&  -  переменная инициализации (возможно нужно будет изменить)
