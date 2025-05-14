#pragma once
#include"Header.h"


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	int number = 0;                // номера начинаютс€ с единицы

	Cell* sosed[2];                // ƒве €чеки-соседи грани (необходимо, чтобы нормаль смотрела от первой ко второй €чейке)

};

