#pragma once
#include "Header.h"

enum class Type_Gran {
	Us,    // 0    ������� �����
	Inner_Hard,      // 1       ��������� ����� �� ���������� ����� � ������� ���������� ���������
	Outer_Hard,      // 2       ��������� ����� �� ������� ������� � ������� ���������� ���������
	Outer_Soft      // 3       ��������� ����� �� ������� ������� � ������� ���������� ���������
};

class Gran
{
public:
	vector<Yzel*> yzels;  // ���� ����� (���� ������ ���� ����������� �� �����)
	vector<Cell*> cells;  // ��� ������ ����� (������� ������ �������� �� ������ �� ������
	// ������������ ����� �� ����� ���������� ��� ����, ����� ����� ��������� ������ ��� ����� 
	// ������������ ����� �����
	vector<Cell*> cells_TVD;  // ������ ��� ����� ��� ��������� �� �����

	int number = 0;                // ������ ���������� � �������
	Type_Gran type = Type_Gran::Us;  // �� ��������� ������ �������

	unordered_map<string, double> parameters;   // ��������� �� ����� (��� ����� ���� ��������
	// ���������� ����� (������� �������� � �.�.), ����� ���� �������� �������

	double normal[2][3];           // ������� ����� (����� � ���������� � ������o�� ������ �������)
	double center[2][3];           // ����� ����� (����� � ���������� � ������o�� ������ �������)
	double area[2];           // ������� ����� (����� � ���������� � ������o�� ������ �������)

	void Culc_measure(unsigned short int st_time);
	// ��������� normal, center, area
	// ��� ���������� ������ ����������� ������� ������ ����� ������ ���� ���������

	double func_R(unsigned short int i_time); // ���������� �� ������ �� ������ ���������

	// ������� ��������� �� ������ �����
	friend bool areCellsEqual(const Gran& cell1, const Gran& cell2);
	friend bool areCellsEqual(const Gran* cell1, const Gran* cell2);
	friend bool areCellsEqual_my(const Gran* cell1, const Gran* cell2);

	Gran();
};

