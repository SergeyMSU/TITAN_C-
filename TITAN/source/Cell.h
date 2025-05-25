#pragma once
#include"Header.h"


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	int number = 0;                // ������ ���������� � �������
	bool is_inner = false;         // �������� �� ������ ���������� (������� ��������� ��������)

	std::array < unordered_map<string, double>, 2> parameters;   // ��������� � ������ (��� ����� ���� ��������
	// ���������� �����, �������� ����������� � �.�.
	// "rho", 'p', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz'
	// 'Q' - ������ ��� �������� � ����������� HP

	Cell* Get_Sosed(Gran* gr); // �������� ������ ����� ������ �����
	// ���������� nullptr, ���� ������ ����� ������ ����� ���, �.�. ����� ���������
	// ��� ������ ������� ������ ����� ������ ���� ���������

	friend Cell* Get_Sosed(Cell* C, Gran* gr);
	// ��� ����������, �� ��������� ��� �� ���������, �� ������� ������� �����

	double center[2][3];           // ����� ����� (����� � ���������� � ������o�� ������ �������)
	double volume[2];

	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 1);

	double func_R(unsigned short int i_time); // ���������� �� ������ �� ������ ���������

};

