#pragma once
using namespace std;
#include "Header.h"
#include <string>
#include <vector>
class Cell;


class point
{

public:
	double x, y, z;
	double Vx, Vy, Vz;
	int Vnum;
	int id;
	int number;
	vector<point*> point_soseds;    
	vector<Cell*> cells;


	point(const double & x, const double& y, const double& z);


	double get_radius(void);

	friend double Pdistance(const point* A, const point* B);
};

