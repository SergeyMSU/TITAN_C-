#pragma once

#include "Header.h"


class Cell
{
public:
	double xc;             // Координаты центра ячейки
	double yc;
	double zc;
	int number = 0;          // Нумерация всех ячеек начинается с 1

	vector<Yzel*> Yzels;          // Все узлы ячейки
	vector<Cell*> Soseds;          // Ячейки - соседи

	void Set_center();            // Посчитать координаты центра ячейки
};

