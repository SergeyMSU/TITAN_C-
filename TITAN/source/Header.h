#pragma once


// �������� ���� �������
class Geo_param;
class Setka;
class Yzel;
class Luch;
class Cell;
class Gran;
class Surfaces;


const double const_pi= 3.14159265358979323846;

#define kv(x) ((x) * (x))
#define kvv(x, y, z) ((x) * (x) + (y) * (y) + (z) * (z))

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


// Boost ���������� (���� ���������� � ����������� ��������)
#include "boost/multi_array.hpp"

#include "Setka.h"
#include "Help_funk.h"
#include "Geo_param.h"
#include "Luch.h"
#include "Gran.h"
#include "Cell.h"
#include "Yzel.h"
#include "Surfaces.h"


using namespace std;





// �������� ���� ����� �����������

// &INIT&  -  ���������� ������������� (�������� ����� ����� ��������)