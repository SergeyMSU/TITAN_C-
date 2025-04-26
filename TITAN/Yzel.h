#pragma once
#include "Header.h"

class Yzel
{
public:
	double coord[2][3];   // ¬ременной слой координат, три координаты
	class Luch* luch;           // Ћуч на котором расположен узел
	int number;                // номера начинаютс€ с единицы

	Yzel();
	Yzel(const double& a, const double& b, const double& c);

	double func_R(int i_time); // –ассто€ние до начала координат
};

