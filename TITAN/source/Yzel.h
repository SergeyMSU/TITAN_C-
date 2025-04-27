#pragma once
#include "Header.h"

class Yzel
{
public:
	double coord[2][3];   // [Временной слой координат, три координаты]
	class Luch* luch;           // Луч на котором расположен узел (каждый узел может лежать только на одном луче)
	int number = 0;                // номера начинаются с единицы

	Yzel();
	Yzel(const double& a, const double& b, const double& c);

	double func_R(int i_time); // Расстояние до начала координат
};

