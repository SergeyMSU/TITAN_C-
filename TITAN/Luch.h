#pragma once

#include"Header.h"

using namespace std;


class Luch
{
public:
	string type;
	static Geo_param* geo;
	map<string, double> parameters;    // ��������� ���� (����� ���� ���� ��� ������� ����)
	// "the" - ���� � ����������� �� �� ��� � - ��� ����������� ������ - �������� ����
	// "phi" - ���� � ����������� �� - ������������ ����
	//Yzel* start;            
	//Yzel* end;              
	vector<Yzel*> Yzels;
	vector<Yzel*> Yzels_opor;
	//vector<Luch*> Luch_soseds;    // ����-������ (�� ��� ����� ������ �������� �����)

	void dvigenie(int i_time);
};
