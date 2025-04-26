#include "Cell.h"

Cell::Cell(const int& a1, const int& a2, const int& a3, const int& a4)
{
	this->int_points.reserve(4);
	this->int_points.push_back(a1);
	this->int_points.push_back(a2);
	this->int_points.push_back(a3);
	this->int_points.push_back(a4);
}

Cell::Cell(point* a1, point* a2, point* a3, point* a4)
{
	this->points.reserve(4);
	this->points.push_back(a1);
	this->points.push_back(a2);
	this->points.push_back(a3);
	this->points.push_back(a4);
}

void Cell::Set_center()
{
	xc = 0.0;
	yc = 0.0;

	for (auto& i : this->points)
	{
		xc = xc + i->x;
		yc = yc + i->y;
	}

	xc /= this->points.size();
	yc /= this->points.size();
}
