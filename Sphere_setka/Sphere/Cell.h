#pragma once
using namespace std;
#include "Header.h"
#include <string>
#include <vector>

class point;

class Cell
{
public:
	vector<point*> points;    // точки €чейки расположены по кругу
	vector<Cell*> soseds;   
	vector<int> int_points;
	double xc;             //  оординаты центра €чейки
	double yc;
	int number;

	Cell(const int& a1, const int& a2, const int& a3, const int& a4);
	Cell(point* a1, point* a2, point* a3, point* a4);

	void Set_center();            // ѕосчитать координаты центра €чейки
};

