#pragma once
#include "Header.h"

class MK_particle
{
public:
	Eigen::Vector3d Vel;    // Скорость частицы
	Eigen::Vector3d coord;  // Положение частицы
	short int sort;         // Сорт частицы
	double mu;              // Вес частицы
	Cell* cel;              // Ячейка, в которой находится частица

};

