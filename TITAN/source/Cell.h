#pragma once
#include"Header.h"

class Cell
{
public:
	vector<Yzel*> yzels;
	int number = 0;                // ������ ���������� � �������

	double normal[2][3];           // ������� ����� (����� � ���������� � ��������� ����� �������)
	Cell* sosed[2];                // ��� �����-������ ����� (����������, ����� ������� �������� �� ������ �� ������ ������)
};

