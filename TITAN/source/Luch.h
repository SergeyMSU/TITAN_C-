#pragma once

#include"Header.h"

using namespace std;


class Luch
{
public:
	string type;        // "A_Luch", "A2_Luch", .....
	static Geo_param* geo;             // ������� ��� ������������� ������ Setka
	map<string, double> parameters;    // ��������� ���� (����� ���� ���� ��� ������� ����)
	// "the" - ���� � ����������� �� �� ��� � - ��� ����������� ������ - �������� ����
	// "phi" - ���� � ����������� �� - ������������ ����
	//Yzel* start;            
	//Yzel* end;              
	vector<Yzel*> Yzels;
	vector<Yzel*> Yzels_opor;
	//vector<Luch*> Luch_soseds;    // ����-������ (�� ��� ����� ������ �������� �����)

	void dvigenie(int i_time);

	Yzel* get_yzel_near_opor(int num_opor, int shift);
	// ����� ���� ����� ������-�� �������� � �����-�� �������
};
