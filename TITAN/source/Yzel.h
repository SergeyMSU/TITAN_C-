#pragma once
#include "Header.h"
#include <vector>

class Gran;

enum class Type_yzel{
	TS,    // 0    ���� �� TS
	HP, // 1       ���� �� HP
	BS, // 2       ���� �� BS
	Us, // 3       ������� ���� (�� ���������������), 
	//����� ���� �� ������� ����� ������ (���� ��� ��� ����������� �������)
	// ����� ���� � ������ ��� ����� ��������� HP � BS
	Zone_1,  // 4   ���� �� TS
	Zone_2,  // 5   ���� ����� TS � HP
	Zone_3,  // 6   ���� ����� HP � BS
	Zone_4  // 7   ���� ����� BS
};

class Yzel
{
public:
	double coord[2][3];   // [��������� ���� ���������, ��� ����������]
	//Luch* luch;           // ��� �� ������� ���������� ���� (������ ���� ����� ������ ������ �� ����� ����)
	int number = 0;                // ������ ���������� � �������
	vector<Gran*> grans;
	Type_yzel type = Type_yzel::Us;  // �� ��������� ������ ������� ����
	bool is_inner = false;

	double velocity[3];
	int num_velocity = 0;
	mutex mut;

	unordered_map<string, Yzel*> Yzel_sosed_sglag;
	// ���� - ������ ��� ����������� �����������
	// AA33  AA3  (AA)  AA1  AA11 - ����� HP  ����������� � ������
	// A22  AA2  (AA)  AA4  AA44  - �������� HP
	std::unordered_set<Yzel*> Yzel_sosed_sglag2;

	unordered_map<string, double> parameters;   // ��������� � ���� (��� ������������ ����)

	Yzel();
	Yzel(const double& a, const double& b, const double& c);

	double func_R(int i_time); // ���������� �� ������ ���������
	double func_Ryz(int i_time); // ���������� �� ������ ���������

	friend double Yzel_distance(Yzel* A, Yzel* B, int time); // ���������� ����� ����� ������ � ������������
	friend double Yzel_distance_x(Yzel* A, Yzel* B, int time); // ���������� ����� ����� ������ � ������������ ��� ����� � ����������
};
