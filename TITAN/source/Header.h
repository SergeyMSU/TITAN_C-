#pragma once


// �������� ���� �������
class Geo_param;
class Yzel;
class Luch;
class Cell;


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

using namespace std;





// �������� ���� ����� �����������

// &INIT&  -  ���������� ������������� (�������� ����� ����� ��������)