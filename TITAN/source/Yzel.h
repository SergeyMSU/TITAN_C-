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
	double func_Ryz(int i_time); // Расстояние до начала координат

	friend double Yzel_distance(Yzel* A, Yzel* B, int time); // расстояние между двумя узлами в пространстве
	friend double Yzel_distance_x(Yzel* A, Yzel* B, int time); // расстояние между двумя узлами в пространстве без учёта х координаты
};
