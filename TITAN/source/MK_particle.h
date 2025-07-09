#pragma once
#include "Header.h"

class MK_particle
{
public:
	double Vel[3];    // �������� �������
	double coord[3];  // ��������� �������
	short int sort;         // ���� �������  1, 2, 3, 4
	double mu;              // ��� �������
	double KSI;              // ��� ���������� ����� ���������� ��������
	double I_do;             // ��� ��������� ����� ���������� �������
	Cell* cel;              // ������, � ������� ��������� �������

	MK_particle();

	void AddVel(const double& a, const double& b, const double& c);
	void AddVel(const Eigen::Vector3d& a);
	void Addcoord(const double& a, const double& b, const double& c);
	void Addcoord(const Eigen::Vector3d& a);

	void Move(const Eigen::Vector3d& a); // ����������� ����� �� ������ a

	double Vel_norm(void);

};

