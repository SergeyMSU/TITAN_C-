#pragma once
#include "Header.h"

class Yzel
{
public:
	double coord[2][3];   // [��������� ���� ���������, ��� ����������]
	class Luch* luch;           // ��� �� ������� ���������� ���� (������ ���� ����� ������ ������ �� ����� ����)
	int number = 0;                // ������ ���������� � �������

	Yzel();
	Yzel(const double& a, const double& b, const double& c);

	double func_R(int i_time); // ���������� �� ������ ���������
};

