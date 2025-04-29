#pragma once
#include "Header.h"


class Surfaces // Файл считывания сетки из записанных файлов
// первоначально сделан для старого формата (из программы MIK на фортране)
{
public:
	string name;
	map<string, int> parameters;

	vector<double> phi_angle;

	// Далее первая координата массивов - цилинтрический угол phi (так как вся сетка 
	// в старой программе получена цилиндрическим вращением

	// Головная область - все поверхности ражиальны
	boost::multi_array<double, 2> the_angle;
	boost::multi_array<double, 2> TS_radial;
	boost::multi_array<double, 2> HP_radial;
	boost::multi_array<double, 2> BS_radial;

	// TS в задней части - она радиальна
	boost::multi_array<double, 2> the_angle2; 
	boost::multi_array<double, 2> TS_radial2;

	// HP в задней части - она цилиндрическая
	boost::multi_array<double, 2> x_cilindr;  // x < 0
	boost::multi_array<double, 2> HP_cilindr;

	~Surfaces();

	void Read_old(string nam);

	double Get_TS(const double& phi, const double& the);
	double Get_HP(const double& phi, const double& the, int met);
};

