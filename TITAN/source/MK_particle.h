#pragma once
#include "Header.h"

class MK_particle
{
public:
	Eigen::Vector3d Vel;    // �������� �������
	Eigen::Vector3d coord;  // ��������� �������
	short int sort;         // ���� �������
	double mu;              // ��� �������
	Cell* cel;              // ������, � ������� ��������� �������

};

