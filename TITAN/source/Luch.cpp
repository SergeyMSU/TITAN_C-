#include "Luch.h"

// Инициализация статического поля (обязательно вне класса!)
Geo_param* Luch::geo = nullptr;  // Можно инициализировать nullptr

void Luch::dvigenie(int i_time)
// i_time - 0 или 1 - какую временную кардинату меняем? Считаем, что для опорных точек эта координата уже поменяна
{
	if (this->type == "A_Luch" || this->type == "A2_Luch")
	{
		double R0 = this->geo->R0;
		double R1 = this->geo->R1;
		double R2 = this->Yzels_opor[1]->func_R(i_time);
		double R3 = this->Yzels_opor[2]->func_R(i_time);
		double R4 = this->Yzels_opor[3]->func_R(i_time);
		double R5 = this->geo->R5;

		int M0 = this->geo->M0;
		int M1 = this->geo->M1;
		int M2 = this->geo->M2;
		int M3 = this->geo->M3;
		int M4 = this->geo->M4;

		if (this->parameters.find("phi") == this->parameters.end()) cout << "9573676612  ERROR" << endl;
		if (this->parameters.find("the") == this->parameters.end()) cout << "0452198906  ERROR" << endl;

		double phi = this->parameters["phi"];
		double the = this->parameters["the"];
		double r;
		int num = 0;

		for (int j = 0; j < M0; j++)
		{
			r = R0 + j * (R1 - R0) / M0;
			this->Yzels[j]->coord[i_time][0] = r * cos(the);
			this->Yzels[j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M0 + 1;

		for (int j = 0; j < M1; j++)
		{
			r = R1 + (j + 1) * (R2 - R1) / (M1 + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M1 + 1;

		for (int j = 0; j < M2; j++)
		{
			r = R2 + (j + 1) * (R3 - R2) / (M2 + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M2 + 1;

		for (int j = 0; j < M3; j++)
		{
			r = R3 + (j + 1) * (R4 - R3) / (M3 + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M3 + 1;

		for (int j = 0; j < M4; j++)
		{
			r = R4 + (j + 1) * (R5 - R4) / (M4 + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}

		//this->Yzels[num + M4]->coord[i_time][0] = R5 * cos(the);
		//this->Yzels[num + M4]->coord[i_time][1] = R5 * sin(the) * cos(phi);
		//this->Yzels[num + M4]->coord[i_time][2] = R5 * sin(the) * sin(phi);

	}
	else if (this->type == "C_Luch" || this->type == "C2_Luch")
	{
		double R0 = this->geo->R0;
		double R1 = this->geo->R1;
		double R2 = this->Yzels_opor[1]->func_R(i_time);
		double L6 = this->geo->L6;
		double L7 = this->geo->L7;

		int M0 = this->geo->M0;
		int M1 = this->geo->M1;
		int M11 = this->geo->M11;
		int N4 = this->geo->N4;
		int N5 = this->geo->N5;
		int dd5 = this->geo->dd5;
		int dd6 = this->geo->dd6;

		if (this->parameters.find("phi") == this->parameters.end()) cout << "0099116454  ERROR" << endl;
		if (this->parameters.find("the") == this->parameters.end()) cout << "1667675420  ERROR" << endl;

		double phi = this->parameters["phi"];
		double the = this->parameters["the"];
		//cout << phi << " " << the << endl;
		double r, dr;
		int num = 0;
		r = 0.0;

		for (int j = 0; j < M0; j++)
		{
			r = R0 + j * (R1 - R0) / M0;
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M0 + 1;


		for (int j = 0; j < M1; j++)
		{
			r = R1 + (j + 1) * (R2 - R1) / (M1 + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		dr = (R2 - R1) / (M1 + 1);
		num += M1 + 1;

		for (int j = 0; j < M11; j++)
		{
			r = R2 + dr * (j + 1); // Здесь не правильное 
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M11;

		double x0 = r * cos(the); // Это координаты с последней итерации предыдущего цикла
		double y0 = r * sin(the);
		double x1, y1, t1, t2, tt1, tt2, y, z;

		t1 = (cos(the) - 1.0) / 2.0 * dd5;
		tt1 = sin(the) * dd5;

		t2 = -1.0 * dd6;
		tt2 = 0.0;
		x1 = L6;
		y1 = y0;

		// Строим сплайн (на концах заданы две точки и два вектора)
		double a1, b1, c1, d1, a2, b2, c2, d2;

		a1 = x0;
		b1 = t1;
		c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
		d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

		a2 = y0;
		b2 = tt1;
		c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
		d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

		for (int j = 0; j < N4; j++)
		{
			double s = 1.0 * (j + 1) / (N4 + 1);
			this->Yzels[num + j]->coord[i_time][0] = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
			y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
			z = y * sin(phi);
			y = y * cos(phi);
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}

		num += N4;
		for (int j = 0; j < N5; j++)
		{
			double s = 1.0 * (j + 1) / (N4 + 1);
			this->Yzels[num + j]->coord[i_time][0] = L6 + (j) * (L7 - L6) / (N5 - 1);
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}

	}
}
