#pragma once
#include"Header.h"


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	int number = 0;                // номера начинаютс€ с единицы

	Cell* sosed[2];                // ƒве €чеки-соседи грани (необходимо, чтобы нормаль смотрела от первой ко второй €чейке)
	
	double center[2][3];           // ÷ентр грани (также в предыдущий и следуюoий момент времени)
	double volume[2];

	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 0);

};

