#pragma once
#include <math.h>
#include "Header.h"

double polar_angle(const double& x, const double& y);

double triangleArea3D(
	const double& x1, const double& y1, const double& z1,
	const double& x2, const double& y2, const double& z2,
	const double& x3, const double& y3, const double& z3);

void crossProductFast(
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz,
	double& rx, double& ry, double& rz);

double scalarProductFast(
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz);

double tetrahedronVolume(
	const double& x1, const double& y1, const double& z1,  // Вершина A
	const double& x2, const double& y2, const double& z2,  // Вершина B
	const double& x3, const double& y3, const double& z3,  // Вершина C
	const double& x4, const double& y4, const double& z4);  // Вершина D

