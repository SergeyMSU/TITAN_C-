#pragma once
#include"Header.h"

class Cell
{
public:
	vector<Yzel*> yzels;
	int number = 0;                // номера начинаютс€ с единицы

	double normal[2][3];           // Ќормаль грани (также в предыдущий и следуюзий момнт времени)
	Cell* sosed[2];                // ƒве €чеки-соседи грани (необходимо, чтобы нормаль смотрела от первой ко второй €чейке)
};

