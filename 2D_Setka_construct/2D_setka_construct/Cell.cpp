#include "Cell.h"

void Cell::Set_center()
{
	xc = 0.0;
	yc = 0.0;
	zc = 0.0;

	for (auto& i : this->Yzels)
	{
		xc = xc + i->coord[0][0];
		yc = yc + i->coord[0][1];
		zc = zc + i->coord[0][2];
	}

	xc /= this->Yzels.size();
	yc /= this->Yzels.size();
	zc /= this->Yzels.size();
}
