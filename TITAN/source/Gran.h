#pragma once
#include "Header.h"

enum class Type_Gran {
	Us,    // 0    обычная грань
	Inner_Hard,      // 1       граничная грань на внутренней сфере с жёсткими граничными условиями
	Outer_Hard,      // 2       граничная грань на внешней границе с жёсткими граничными условиями
	Outer_Soft      // 3       граничная грань на внешней границе с мягкими граничными условиями
};

class Gran
{
public:
	vector<Yzel*> yzels;  // Узлы грани (узлы должны быть расположены по кругу)
	vector<Cell*> cells;  // две ячейки грани (нормаль должна смотреть от первой ко второй
	// Расположение узлов по кругу необходимо для того, чтобы любые выбранные подряд две точки 
	// образовывали ребро грани
	vector<Cell*> cells_TVD;  // Ячейки для сноса ТВД процедуоы на грань

	int number = 0;                // номера начинаются с единицы
	Type_Gran type = Type_Gran::Us;  // по умолчанию создаём обычный

	unordered_map<string, double> parameters;   // Параметры на грани (тут могут быть значения
	// плазменных полей (большие величины и т.д.), могут быть значения потоков

	double normal[2][3];           // Нормаль грани (также в предыдущий и следуюoий момент времени)
	double center[2][3];           // Центр грани (также в предыдущий и следуюoий момент времени)
	double area[2];           // Площадь грани (также в предыдущий и следуюoий момент времени)

	void Culc_measure(unsigned short int st_time);
	// вычисляет normal, center, area
	// для корректной работы определения нормали центры ячеек должны быть актуальны

	double func_R(unsigned short int i_time); // Расстояние от центра до начала координат

	// Функция сравнения по набору узлов
	friend bool areCellsEqual(const Gran& cell1, const Gran& cell2);
	friend bool areCellsEqual(const Gran* cell1, const Gran* cell2);
	friend bool areCellsEqual_my(const Gran* cell1, const Gran* cell2);

	Gran();
};

