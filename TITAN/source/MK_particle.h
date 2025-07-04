#pragma once
#include "Header.h"

class MK_particle
{
public:
	Eigen::Vector3d Vel;    // �������� �������
	Eigen::Vector3d coord;  // ��������� �������
	short int sort;         // ���� �������  1, 2, 3, 4
	double mu;              // ��� �������
	double KSI;              // ��� ���������� ����� ���������� ��������
	double I_do;             // ��� ��������� ����� ���������� �������
	Cell* cel;              // ������, � ������� ��������� �������

	MK_particle();

};

