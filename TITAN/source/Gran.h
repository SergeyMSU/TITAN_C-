#pragma once
#include "Header.h"


class Gran
{
public:
	vector<Yzel*> yzels;  // ���� ����� (���� ������ ���� ����������� �� �����)
	vector<Cell*> cells;  // ��� ������ ����� (������� ������ �������� �� ������ �� ������
	// ������������ ����� �� ����� ���������� ��� ����, ����� ����� ��������� ������ ��� ����� 
	// ������������ ����� �����
	int number = 0;                // ������ ���������� � �������

	double normal[2][3];           // ������� ����� (����� � ���������� � ������o�� ������ �������)
	double center[2][3];           // ����� ����� (����� � ���������� � ������o�� ������ �������)
	double area[2];           // ������� ����� (����� � ���������� � ������o�� ������ �������)

	void Culc_measure(unsigned short int st_time);
	// ��������� normal, center, area

	// ������� ��������� �� ������ �����
	friend bool areCellsEqual(const Gran& cell1, const Gran& cell2);
	friend bool areCellsEqual(const Gran* cell1, const Gran* cell2);
	friend bool areCellsEqual_my(const Gran* cell1, const Gran* cell2);

	Gran();
};

