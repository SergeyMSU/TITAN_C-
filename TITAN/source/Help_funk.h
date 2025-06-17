#pragma once
#include <math.h>
#include "Header.h"

void Sootnosheniya(const double& rho, const double& p, const double& rho_He,
	const double& rho_Pui, const double& T_Pui,
	const short int& zone,
	double& rho_Th, double& rho_E, double& p_Th, double& p_Pui,
	double& T_Th, double& T_E);

short int signum(const double& x);

double minmod(const double& x, const double& y);

double linear(const double& x1, const double& t1, const double& x2,
	const double& t2, const double& x3, const double& t3,
	const double& y);

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

