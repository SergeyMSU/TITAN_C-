#pragma once
#include"Header.h"


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	int number = 0;                // ������ ���������� � �������

	Cell* sosed[2];                // ��� �����-������ ����� (����������, ����� ������� �������� �� ������ �� ������ ������)

};

