#pragma once
#include <math.h>
#include "Header.h"

double rbfKernel(double r, double epsilon = 1.0);

double polar_angle(const double& x, const double& y);

void dekard_skorost(const double& z, const double& x, const double& y,
	const double& Vr, const double& Vphi, const double& Vtheta,
	double& Vz, double& Vx, double& Vy);

double triangleArea3D(
	const double& x1, const double& y1, const double& z1,
	const double& x2, const double& y2, const double& z2,
	const double& x3, const double& y3, const double& z3);

void crossProductFast(
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz,
	double& rx, double& ry, double& rz);

void get_bazis(const Eigen::Vector3d& n, Eigen::Vector3d& t, Eigen::Vector3d& m);
// �� ���������� ������� n ������� ��� �������, ���������� � ��� ������ ����������������� ������

double scalarProductFast(
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz);

double tetrahedronVolume(
	const double& x1, const double& y1, const double& z1,  // ������� A
	const double& x2, const double& y2, const double& z2,  // ������� B
	const double& x3, const double& y3, const double& z3,  // ������� C
	const double& x4, const double& y4, const double& z4);  // ������� D

