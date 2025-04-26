#pragma once
#include "Header.h"

class Yzel
{
public:
	double coord[2][3];   // Временной слой координат, три координаты в два момента времени
	Luch* luch;           // Луч на котором расположен узел
	vector<Cell*> Cells;  // Ячейки, которым принадлежит узел
	int number = 0;

	int zavisimost;       // Степень зависимости узла от остальных узлов (у исходных зависимость = 0
	vector<Yzel*> zav_yzels;   // Узлы, от которых есть зависимость
	vector<double> zav_koeff;  // Коэффициенты зависимости узлов (их сумма равна единице)

	Yzel();
	Yzel(const double& a, const double& b, const double& c);
};

