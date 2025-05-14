#pragma once
#include "Header.h"


class Gran
{
public:
	vector<Yzel*> yzels;  // Узлы грани (узлы должны быть расположены по кругу)
	vector<Cell*> cells;  // две ячейки грани (нормаль должна смотреть от первой ко второй
	// Расположение узлов по кругу необходимо для того, чтобы любые выбранные подряд две точки 
	// образовывали ребро грани
	int number = 0;                // номера начинаются с единицы

	double normal[2][3];           // Нормаль грани (также в предыдущий и следуюoий момент времени)
	double center[2][3];           // Центр грани (также в предыдущий и следуюoий момент времени)
	double area[2];           // Площадь грани (также в предыдущий и следуюoий момент времени)

	void Culc_measure(unsigned short int st_time);
	// вычисляет normal, center, area

	// Функция сравнения по набору узлов
	friend bool areCellsEqual(const Gran& cell1, const Gran& cell2);
	friend bool areCellsEqual(const Gran* cell1, const Gran* cell2);
	friend bool areCellsEqual_my(const Gran* cell1, const Gran* cell2);

	Gran();
};

