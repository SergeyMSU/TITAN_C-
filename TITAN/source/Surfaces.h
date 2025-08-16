#pragma once
#include "Header.h"


class Surfaces // ���� ���������� ����� �� ���������� ������
// ������ ��� ������� ������� (�� ��������� MIK �� ��������)
{
public:
	string name;
	map<string, int> parameters;

	vector<double> phi_angle;

	// ����� ������ ���������� �������� - �������������� ���� phi (��� ��� ��� ����� 
	// � ������ ��������� �������� �������������� ���������

	// �������� ������� - ��� ����������� ���������
	boost::multi_array<double, 2> the_angle;
	boost::multi_array<double, 2> TS_radial;
	boost::multi_array<double, 2> HP_radial;
	boost::multi_array<double, 2> BS_radial;

	// TS � ������ ����� - ��� ���������
	boost::multi_array<double, 2> the_angle2; 
	boost::multi_array<double, 2> TS_radial2;

	// HP � ������ ����� - ��� ��������������
	boost::multi_array<double, 2> x_cilindr;  // x < 0
	boost::multi_array<double, 2> HP_cilindr;

	~Surfaces();

	void Read_old(string nam);

	double Get_TS(const double& phi, const double& the);
	double Get_HP(const double& phi, const double& the, int met);
	double Get_BS(const double& phi, const double& the);

	// ������������� �����������

	void Print_TS(); // �������� ����� TS
	void Print_HP(); // �������� ����� TS
	void Print_BS(); // �������� ����� TS
};

