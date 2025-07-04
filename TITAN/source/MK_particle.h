#pragma once
#include "Header.h"

class MK_particle
{
public:
	Eigen::Vector3d Vel;    // Скорость частицы
	Eigen::Vector3d coord;  // Положение частицы
	short int sort;         // Сорт частицы  1, 2, 3, 4
	double mu;              // Вес частицы
	double KSI;              // Для вычисления длины свободного провбега
	double I_do;             // Для розыгрыша длины свободного пробега
	Cell* cel;              // Ячейка, в которой находится частица

	MK_particle();

};

