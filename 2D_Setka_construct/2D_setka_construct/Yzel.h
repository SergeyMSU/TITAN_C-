#pragma once
#include "Header.h"

class Yzel
{
public:
	double coord[2][3];   // ��������� ���� ���������, ��� ���������� � ��� ������� �������
	Luch* luch;           // ��� �� ������� ���������� ����
	vector<Cell*> Cells;  // ������, ������� ����������� ����
	int number = 0;

	int zavisimost;       // ������� ����������� ���� �� ��������� ����� (� �������� ����������� = 0
	vector<Yzel*> zav_yzels;   // ����, �� ������� ���� �����������
	vector<double> zav_koeff;  // ������������ ����������� ����� (�� ����� ����� �������)

	Yzel();
	Yzel(const double& a, const double& b, const double& c);
};

