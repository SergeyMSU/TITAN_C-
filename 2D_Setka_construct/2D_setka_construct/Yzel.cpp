#include "Yzel.h"

Yzel::Yzel()
{
	this->coord[0][0] = this->coord[0][1] = this->coord[0][2] = 0.0;
	this->coord[1][0] = this->coord[1][1] = this->coord[1][2] = 0.0;
	this->luch = nullptr;
	this->zavisimost = 0;
}

Yzel::Yzel(const double& a, const double& b, const double& c)
{
	this->coord[0][0] = a;
	this->coord[0][1] = b;
	this->coord[0][2] = c;

	this->coord[1][0] = a;
	this->coord[1][1] = b;
	this->coord[1][2] = c;

	this->zavisimost = 0;

	this->luch = nullptr;
}
