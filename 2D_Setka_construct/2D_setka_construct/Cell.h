#pragma once

#include "Header.h"


class Cell
{
public:
	double xc;             // ���������� ������ ������
	double yc;
	double zc;
	int number = 0;          // ��������� ���� ����� ���������� � 1

	vector<Yzel*> Yzels;          // ��� ���� ������
	vector<Cell*> Soseds;          // ������ - ������

	void Set_center();            // ��������� ���������� ������ ������
};

