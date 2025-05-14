#pragma once
#include "Header.h"
#include <vector>

class Gran;

class Yzel
{
public:
	double coord[2][3];   // [��������� ���� ���������, ��� ����������]
	//Luch* luch;           // ��� �� ������� ���������� ���� (������ ���� ����� ������ ������ �� ����� ����)
	int number = 0;                // ������ ���������� � �������
	vector<Gran*> grans;


	Yzel();
	Yzel(const double& a, const double& b, const double& c);

	double func_R(int i_time); // ���������� �� ������ ���������
	double func_Ryz(int i_time); // ���������� �� ������ ���������

	friend double Yzel_distance(Yzel* A, Yzel* B, int time); // ���������� ����� ����� ������ � ������������
	friend double Yzel_distance_x(Yzel* A, Yzel* B, int time); // ���������� ����� ����� ������ � ������������ ��� ����� � ����������
};
