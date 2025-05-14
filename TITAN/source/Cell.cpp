#include "Cell.h"

void Cell::Culc_center(unsigned short int st_time)
{
	double xc, yc, zc;
	int Ny = this->yzels.size();
	xc = 0.0;
	yc = 0.0;
	zc = 0.0;

	// Вычисляем центр грани
	for (auto& i : this->yzels)
	{
		xc += i->coord[st_time][0];
		yc += i->coord[st_time][1];
		zc += i->coord[st_time][2];
	}
	xc /= Ny;
	yc /= Ny;
	zc /= Ny;
	this->center[st_time][0] = xc;
	this->center[st_time][1] = yc;
	this->center[st_time][2] = zc;
}

void Cell::Culc_volume(unsigned short int st_time, unsigned short int method)
{
	// 0 - быстрый вариант (нужны посчанные площади граней и их центры)
	// 1 - медленный вариант (нужны посчитанные только центры граней)
	// для ровных граней оба методы работают одинаково, но чем более кривая грань, 
	// тем большее расхождения получается.
	// медленный вариант более правильный (но это не точно). 
	// Думаю всегда можно выбирать вариант "0"
	if (method == 0)
	{
		double xg, yg, zg, h;
		double V = 0.0;
		for (auto& i : this->grans)
		{
			xg = i->center[st_time][0];
			yg = i->center[st_time][1];
			zg = i->center[st_time][2];
			h = fabs(scalarProductFast(i->normal[st_time][0], i->normal[st_time][1],
				i->normal[st_time][2],
				this->center[st_time][0] - xg, this->center[st_time][1] - yg,
				this->center[st_time][2] - zg));
			V += h * i->area[st_time];
		}
		this->volume[st_time] = V / 3.0; 
		return;
	}
	if (method == 1)
	{
		int j1;
		double V = 0.0;
		for (auto& i : this->grans)
		{
			for (int j = 0; j < i->yzels.size(); j++)
			{
				j1 = j + 1;
				if (j1 >= i->yzels.size()) j1 = 0;
				V += tetrahedronVolume(this->center[st_time][0], this->center[st_time][1],
					this->center[st_time][2],
					i->center[st_time][0], i->center[st_time][1], i->center[st_time][2],
					i->yzels[j]->coord[st_time][0], i->yzels[j]->coord[st_time][1],
					i->yzels[j]->coord[st_time][2],
					i->yzels[j1]->coord[st_time][0], i->yzels[j1]->coord[st_time][1],
					i->yzels[j1]->coord[st_time][2]);
			}
		}
		this->volume[st_time] = V;
		return;
	}
	else
	{
		cout << "Error  0967452122" << endl;
		return;
	}
}
