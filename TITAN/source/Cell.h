#pragma once
#include"Header.h"


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	int number = 0;                // ������ ���������� � �������

	Cell* sosed[2];                // ��� �����-������ ����� (����������, ����� ������� �������� �� ������ �� ������ ������)
	
	double center[2][3];           // ����� ����� (����� � ���������� � ������o�� ������ �������)
	double volume[2];

	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 0);

};

