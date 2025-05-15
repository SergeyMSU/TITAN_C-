#pragma once
#include"Header.h"


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	int number = 0;                // номера начинаютс€ с единицы

	unordered_map<string, double> parameters;   // ѕараметры в €чейке (тут могут быть значени€
	// плазменных полей, значени€ дивергенций и т.д.
	// "rho", 'p', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz'
	// 'Q' - ма€чок дл€ переноса и определени€ HP

	Cell* sosed[2];                // ƒве €чеки-соседи грани (необходимо, чтобы нормаль смотрела от первой ко второй €чейке)
	
	double center[2][3];           // ÷ентр грани (также в предыдущий и следуюoий момент времени)
	double volume[2];

	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 0);

	double func_R(unsigned short int i_time); // –ассто€ние от центра до начала координат

};

