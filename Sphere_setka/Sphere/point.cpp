#include "point.h"

point::point(const double& x, const double& y, const double& z)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->id = -1;
	this->Vx = 0.0;
	this->Vy = 0.0;
	this->Vz = 0.0;
	this->Vnum = 0;

}

double point::get_radius(void)
{
	return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

double Pdistance(const point* A, const point* B)
{

	return sqrt( (A->x - B->x) * (A->x - B->x) + (A->y - B->y) * (A->y - B->y) + (A->z - B->z) * (A->z - B->z));
}
