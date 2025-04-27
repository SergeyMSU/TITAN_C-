#pragma once


// Описание всех классов
class Geo_param;
class Yzel;
class Luch;
class Cell;
class Gran;


const double const_pi= 3.14159265358979323846;

#define kv(x) ((x) * (x))
#define kvv(x, y, z) ((x) * (x) + (y) * (y) + (z) * (z))

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <math.h>



#include "Help_funk.h"
#include "Geo_param.h"
#include "Luch.h"
#include "Yzel.h"
#include "Setka.h"
#include "Cell.h"
#include "Gran.h"

using namespace std;





// Описание всех типов комметариев

// &INIT&  -  переменная инициализации (возможно нужно будет изменить)