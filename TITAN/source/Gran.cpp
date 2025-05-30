#include "Gran.h"

void Gran::Culc_measure(unsigned short int st_time)
{
	double xc, yc, zc;
	int Ny = this->yzels.size();
	xc = 0.0;
	yc = 0.0;
	zc = 0.0;

	// ¬ычисл€ем центр грани
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

	// ¬ычисл€ем площадь грани
	double S = 0.0;
	int i2;
	for (int i = 0; i < this->yzels.size(); i++)
	{
		i2 = i + 1;
		if (i2 >= this->yzels.size()) i2 = 0;
		S += triangleArea3D(this->yzels[i]->coord[st_time][0], this->yzels[i]->coord[st_time][1],
			this->yzels[i]->coord[st_time][2],
			this->yzels[i2]->coord[st_time][0], this->yzels[i2]->coord[st_time][1],
			this->yzels[i2]->coord[st_time][2],
			xc, yc, zc);
	}
	this->area[st_time] = S;

	// ¬ычисл€ем нормаль грани
	if (Ny != 4) // дл€ других случаем следующий блок может работать не правильно, 
		// нужно отдельно их рассматривать
	{
		cout << "Error  0906764104" << endl;
	}

	double x1, y1, z1;
	double x2, y2, z2;
	double x3, y3, z3;

	x1 = (this->yzels[0]->coord[st_time][0] + this->yzels[1]->coord[st_time][0]) / 2.0;
	y1 = (this->yzels[0]->coord[st_time][1] + this->yzels[1]->coord[st_time][1]) / 2.0;
	z1 = (this->yzels[0]->coord[st_time][2] + this->yzels[1]->coord[st_time][2]) / 2.0;

	x2 = (this->yzels[1]->coord[st_time][0] + this->yzels[2]->coord[st_time][0]) / 2.0;
	y2 = (this->yzels[1]->coord[st_time][1] + this->yzels[2]->coord[st_time][1]) / 2.0;
	z2 = (this->yzels[1]->coord[st_time][2] + this->yzels[2]->coord[st_time][2]) / 2.0;

	x3 = (this->yzels[2]->coord[st_time][0] + this->yzels[3]->coord[st_time][0]) / 2.0;
	y3 = (this->yzels[2]->coord[st_time][1] + this->yzels[3]->coord[st_time][1]) / 2.0;
	z3 = (this->yzels[2]->coord[st_time][2] + this->yzels[3]->coord[st_time][2]) / 2.0;

	double n1, n2, n3, nn;
	crossProductFast(x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1, z3 - z1, n1, n2, n3);
	nn = norm2(n1, n2, n3);
	n1 /= nn;
	n2 /= nn;
	n3 /= nn;

	// Ќужно чтобы нормаль смотрела от первой €чейки ко второй
	if (scalarProductFast(n1, n2, n3, this->cells[0]->center[st_time][0] - xc,
		this->cells[0]->center[st_time][1] - yc, this->cells[0]->center[st_time][2] - zc) > 0.0)
	{
		n1 = -n1;
		n2 = -n2;
		n3 = -n3;
	}

	this->normal[st_time][0] = n1;
	this->normal[st_time][1] = n2;
	this->normal[st_time][2] = n3;
}

Gran::Gran()
{
	this->yzels.reserve(4);
	this->cells.reserve(2);
	this->cells_TVD.reserve(2);
	this->area[0] = this->area[1] = 0.0;
}

double Gran::culc_velosity(short int now1)
{
	int now2 = (now1 + 1) % 2;
	Eigen::Vector3d n, dx;
	n << this->normal[now1][0], this->normal[now1][1], this->normal[now1][2];
	double ddot = 0.0;

	for (auto& i : this->yzels)
	{
		dx(0) = i->coord[now2][0] - i->coord[now1][0];
		dx(1) = i->coord[now2][1] - i->coord[now1][1];
		dx(2) = i->coord[now2][2] - i->coord[now1][2];
		ddot = ddot + dx.dot(n);
		//V = V + dx.dot(n) * n;
	}

	ddot = ddot / this->yzels.size();
	return ddot;
}

double Gran::func_R(unsigned short int i_time)
{
	return norm2(this->center[i_time][0], this->center[i_time][1], this->center[i_time][2]);
}