#include "Luch.h"

// ������������� ������������ ���� (����������� ��� ������!)
Geo_param* Luch::geo = nullptr;  // ����� ���������������� nullptr

void Luch::dvigenie(int i_time)
// i_time - 0 ��� 1 - ����� ��������� ��������� ������? �������, ��� ��� ������� ����� ��� ���������� ��� ��������
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
		int M11 = this->geo->M11;
		int M2 = this->geo->M2;
		int M3 = this->geo->M3;
		int M4 = this->geo->M4;

		double da1 = this->geo->da1;
		double da2 = this->geo->da2;

		if (this->parameters.find("da1") != this->parameters.end())
		{
			da1 = this->parameters["da1"];
		}

		if (this->parameters.find("da2") != this->parameters.end()) da2 = this->parameters["da2"];


		//cout << "da1 = " << da1 << endl;


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
		double dr = (R2 - R1) / (M1 + 1);

		for (int j = 0; j < M11; j++)
		{
			r = R2 + dr * (j + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M11;

		// ������ �������� ����
		double x0, y0, x1, y1, t1, t2, tt1, tt2;

		t1 = cos(the) * da1;
		tt1 = sin(the) * da1;
		t2 = cos(the) * da2;
		tt2 = sin(the) * da2;
		x0 = r * cos(the); // ��� ���������� � ��������� �������� ����������� �����
		y0 = r * sin(the);
		x1 = R3 * cos(the); // ��� ���������� � ��������� �������� ����������� �����
		y1 = R3 * sin(the);

		double a1, b1, c1, d1, a2, b2, c2, d2, x, y, z;

		a1 = x0;
		b1 = t1;
		c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
		d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

		a2 = y0;
		b2 = tt1;
		c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
		d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;


		for (int j = 0; j < M2 - M11; j++)
		{
			double s = 1.0 * (j + 1) / (M2 - M11 + 1);
			x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
			y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
			z = y * sin(phi);
			y = y * cos(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;

			/*r = R2 + dr * M11 + (j + 1) * (R3 - R2 - dr * M11) / (M2 - M11 + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);*/
		}
		num += (M2 - M11) + 1;

		double da3 = this->geo->da3;
		double da4 = this->geo->da4;

		if (this->parameters.find("da3") != this->parameters.end())
		{
			da3 = this->parameters["da3"];
		}

		if (this->parameters.find("da4") != this->parameters.end())
		{
			da4 = this->parameters["da4"];
		}

		t1 = cos(the) * da3;
		tt1 = sin(the) * da3;
		t2 = cos(the) * da4;
		tt2 = sin(the) * da4;
		x0 = R3 * cos(the); // ��� ���������� � ��������� �������� ����������� �����
		y0 = R3 * sin(the);
		x1 = R4 * cos(the); // ��� ���������� � ��������� �������� ����������� �����
		y1 = R4 * sin(the);

		a1 = x0;
		b1 = t1;
		c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
		d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

		a2 = y0;
		b2 = tt1;
		c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
		d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;


		for (int j = 0; j < M3; j++)
		{
			double s = 1.0 * (j + 1) / (M3 + 1);
			x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
			y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
			z = y * sin(phi);
			y = y * cos(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;

			/*r = R3 + (j + 1) * (R4 - R3) / (M3 + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);*/
		}
		num += M3 + 1;

		for (int j = 0; j < M4; j++)
		{
			r = R4 + (j + 1) * (R5 - R4) / (M4 + 1);
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}

		this->Yzels[num + M4]->coord[i_time][0] = R5 * cos(the);
		this->Yzels[num + M4]->coord[i_time][1] = R5 * sin(the) * cos(phi);
		this->Yzels[num + M4]->coord[i_time][2] = R5 * sin(the) * sin(phi);

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
			r = R2 + dr * (j + 1); // ����� �� ���������� 
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M11;

		double x0 = r * cos(the); // ��� ���������� � ��������� �������� ����������� �����
		double y0 = r * sin(the);
		double x1, y1, t1, t2, tt1, tt2, y, z;

		t1 = (cos(the) - 1.0) / 2.0 * dd5;
		tt1 = sin(the) * dd5;

		t2 = -1.0 * dd6;
		tt2 = 0.0;
		x1 = L6;
		y1 = y0;

		// ������ ������ (�� ������ ������ ��� ����� � ��� �������)
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
	else if (this->type == "B_Luch")
	{
		double R0 = this->geo->R0;
		double R1 = this->geo->R1;
		double R5 = this->geo->R5;
		double R2 = this->Yzels_opor[1]->func_R(i_time); // TS
		double R3 = sqrt(kv(this->Yzels_opor[2]->coord[i_time][1]) + 
			kv(this->Yzels_opor[2]->coord[i_time][2])); // HP (���������� �� HP �� ��� ��������)
		//double x3 = this->Yzels_opor[2]->coord[i_time][0];

		double R4 = sqrt(kv(this->Yzels_opor[3]->coord[i_time][1]) +
			kv(this->Yzels_opor[3]->coord[i_time][2])); // BS ������ ��� B ����� ��� ��������� BS
		// ��� ��� ��� �� ����������, ������� ��������� ��� ���������� ���� ������� � ������ ������������


		int M0 = this->geo->M0;
		int M1 = this->geo->M1;
		int M11 = this->geo->M11;
		int M2 = this->geo->M2;
		int M3 = this->geo->M3;
		int M4 = this->geo->M4;

		if (this->parameters.find("phi") == this->parameters.end()) cout << "5690890909  ERROR" << endl;
		if (this->parameters.find("the") == this->parameters.end()) cout << "1237130806  ERROR" << endl;

		double phi = this->parameters["phi"];
		double the = this->parameters["the"];
		double r, x, y, z, dr;
		int num = 0;

		for (int j = 0; j < M0; j++)
		{
			r = R0 + j * (R1 - R0) / M0;
			this->Yzels[j + num]->coord[i_time][0] = r * cos(the);
			this->Yzels[j + num]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[j + num]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M0;

		this->Yzels[num]->coord[i_time][0] = R1 * cos(the);
		this->Yzels[num]->coord[i_time][1] = R1 * sin(the) * cos(phi);
		this->Yzels[num]->coord[i_time][2] = R1 * sin(the) * sin(phi);
		num += 1;

		for (int j = 0; j < M1; j++)
		{
			r = R1 + (j + 1) * (R2 - R1) / (M1 + 1);
			this->Yzels[j + num]->coord[i_time][0] = r * cos(the);
			this->Yzels[j + num]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[j + num]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M1 + 1; // (+1 ��� ������� �����) �� ������� ����� ����� ������� �� �����������, ��� ������ ���� ��� ���������
		dr = (R2 - R1) / (M1 + 1);

		for (int j = 0; j < M11; j++)
		{
			r = R2 + dr * (j + 1); 
			this->Yzels[num + j]->coord[i_time][0] = r * cos(the);
			this->Yzels[num + j]->coord[i_time][1] = r * sin(the) * cos(phi);
			this->Yzels[num + j]->coord[i_time][2] = r * sin(the) * sin(phi);
		}
		num += M11;

		double x3 = r * cos(the); // ���������
		this->Yzels_opor[2]->coord[i_time][0] = x3;
		this->Yzels_opor[3]->coord[i_time][0] = x3;

		// ������ �������� ����
		double x0, y0, x1, y1, t1, t2, tt1, tt2;

		double dd3 = this->geo->dd3;
		if (this->parameters.find("dd3") != this->parameters.end()) dd3 = this->parameters["dd3"];

		double dd4 = this->geo->dd4;
		if (this->parameters.find("dd4") != this->parameters.end()) dd4 = this->parameters["dd4"];
		

		t1 = cos(the) * dd3;
		tt1 = sin(the) * dd3;
		t2 = 0.0;
		tt2 = 1.0 * dd4;
		x0 = r * cos(the); // ��� ���������� � ��������� �������� ����������� �����
		y0 = r * sin(the);
		x1 = x3;
		y1 = R3;

		double a1, b1, c1, d1, a2, b2, c2, d2;

		a1 = x0;
		b1 = t1;
		c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
		d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

		a2 = y0;
		b2 = tt1;
		c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
		d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

		for (int j = 0; j < M2 - M11; j++)
		{
			double s = 1.0 * (j + 1) / (M2 - M11 + 1);
			x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
			y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
			z = y * sin(phi);
			y = y * cos(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}
		num += (M2 - M11) + 1;

		// ����� �� R3 (HP) �� R4 (BS)
		for (int j = 0; j < M3; j++)
		{
			r = R3 + (j + 1) * (R4 - R3) / (M3 + 1);
			//x = x;  // ��������� ���������� ����������� �����
			y = r * cos(phi);
			z = r * sin(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}
		num += M3 + 1;

		for (int j = 0; j < M4; j++)
		{
			r = R4 + (j + 1) * (R5 - R4) / (M4 + 1);
			//x = x;
			y = r * cos(phi);
			z = r * sin(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}
		num += M4;

		y = R5 * cos(phi);
		z = R5 * sin(phi);
		this->Yzels[num]->coord[i_time][0] = x;
		this->Yzels[num]->coord[i_time][1] = y;
		this->Yzels[num]->coord[i_time][2] = z;
	}
	else if (this->type == "D_Luch")
	{	// D ���� ��������� ����������� ����� C - �����, ��� ��� ��� ��������� �� ���
		int M2 = this->geo->M2;
		int M11 = this->geo->M11;
		int M3 = this->geo->M3;
		int M4 = this->geo->M4;

		double R5 = this->geo->R5;

		if (this->parameters.find("phi") == this->parameters.end()) cout << "0990967451  ERROR" << endl;

		double phi = this->parameters["phi"];
		int num = 0;
		double x, y, z, r;

		double x0 = this->Yzels_opor[0]->coord[i_time][0];
		double y0 = sqrt(kv(this->Yzels_opor[0]->coord[i_time][1]) + 
			kv(this->Yzels_opor[0]->coord[i_time][2]));

		double x1 = x0;// this->Yzels_opor[1]->coord[i_time][0];
		this->Yzels_opor[1]->coord[i_time][0] = x0;
		double y1 = sqrt(kv(this->Yzels_opor[1]->coord[i_time][1]) +
			kv(this->Yzels_opor[1]->coord[i_time][2]));

		double x2 = x0;// this->Yzels_opor[2]->coord[i_time][0];
		this->Yzels_opor[2]->coord[i_time][0] = x0;
		double y2 = sqrt(kv(this->Yzels_opor[2]->coord[i_time][1]) +
			kv(this->Yzels_opor[2]->coord[i_time][2]));

		x = x0;

		// ����� �� R2 (TS) �� R3 (HP)
		num += 1;
		for (int j = 0; j < M2 - M11; j++)
		{
			//x = x0 + (x1 - x0) * (j + 1) / (M2 - M11 + 1);
			y = y0 + (y1 - y0) * (j + 1) / (M2 - M11 + 1);
			z = y * sin(phi);
			y = y * cos(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}
		num += (M2 - M11) + 1;

		for (int j = 0; j < M3; j++)
		{
			r = y1 + (j + 1) * (y2 - y1) / (M3 + 1);
			//x = x1;
			y = r * cos(phi);
			z = r * sin(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}
		num += M3 + 1;

		for (int j = 0; j < M4; j++)
		{
			r = y2 + (j + 1) * (R5 - y2) / (M4 + 1);
			//x = x1;
			y = r * cos(phi);
			z = r * sin(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}
		num += M4;

		r = R5;
		//x = x1;
		y = r * cos(phi);
		z = r * sin(phi);
		this->Yzels[num]->coord[i_time][0] = x;
		this->Yzels[num]->coord[i_time][1] = y;
		this->Yzels[num]->coord[i_time][2] = z;
	}
	else if (this->type == "E_Luch")
	{	// E ���� ��������� ����������� ����� C - �����, ��� ��� ��� ��������� �� ���
		int M3 = this->geo->M3;
		int M4 = this->geo->M4;

		double R5 = this->geo->R5;

		if (this->parameters.find("phi") == this->parameters.end()) cout << "4533212122  ERROR" << endl;

		double phi = this->parameters["phi"];
		int num = 0;
		double x, y, z, r;

		//double x0 = this->Yzels_opor[1]->coord[i_time][0];
		double y0 = sqrt(kv(this->Yzels_opor[1]->coord[i_time][1]) +
			kv(this->Yzels_opor[1]->coord[i_time][2]));

		//double x1 = this->Yzels_opor[2]->coord[i_time][0];
		double y1 = sqrt(kv(this->Yzels_opor[2]->coord[i_time][1]) +
			kv(this->Yzels_opor[2]->coord[i_time][2]));

		x = this->Yzels_opor[0]->coord[i_time][0];

		this->Yzels[num]->coord[i_time][0] = x;

		num += 1;
		for (int j = 0; j < M3; j++)
		{
			r = y0 + (j + 1) * (y1 - y0) / (M3 + 1);
			y = r * cos(phi);
			z = r * sin(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}
		num += M3;

		this->Yzels[num]->coord[i_time][0] = x;
		
		num += 1;

		for (int j = 0; j < M4; j++)
		{
			r = y1 + (j + 1) * (R5 - y1) / (M4 + 1);
			y = r * cos(phi);
			z = r * sin(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}
		num += M4;

		r = R5;
		y = r * cos(phi);
		z = r * sin(phi);
		this->Yzels[num]->coord[i_time][0] = x;
		this->Yzels[num]->coord[i_time][1] = y;
		this->Yzels[num]->coord[i_time][2] = z;

	}
	else if (this->type == "H_Luch")
	{	
		int M3 = this->geo->M3;
		int M4 = this->geo->M4;

		double R5 = this->geo->R5;

		if (this->parameters.find("phi") == this->parameters.end()) cout << "7635497867  ERROR" << endl;

		double phi = this->parameters["phi"];
		int num = 0;
		double x, y, z, r;

		// ������ �������� ����
		double x0, y0, x1, y1, t1, t2, tt1, tt2;

		auto yz3 = this->Yzels_opor[0];
		auto yz1 = this->Yzels_opor[1];
		auto yz4 = this->Yzels_opor[2];
		auto yz2 = this->Yzels_opor[3];

		double dd7 = this->geo->dd7;
		double dd8 = this->geo->dd8;

		t1 = (yz1->coord[i_time][0] - yz3->coord[i_time][0]) * dd7;
		tt1 = (yz1->func_Ryz(i_time) - yz3->func_Ryz(i_time)) * dd7;
		t2 = (yz4->coord[i_time][0] - yz2->coord[i_time][0]) * dd8;
		tt2 = (yz4->func_Ryz(i_time) - yz2->func_Ryz(i_time)) * dd8;
		x0 = yz1->coord[i_time][0];       
		y0 = yz1->func_Ryz(i_time);
		x1 = yz2->coord[i_time][0];        
		y1 = yz2->func_Ryz(i_time);

		double a1, b1, c1, d1, a2, b2, c2, d2;

		a1 = x0;
		b1 = t1;
		c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
		d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

		a2 = y0;
		b2 = tt1;
		c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
		d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

		num += 1;
		for (int j = 0; j < 4; j++)
		{
			double s = 1.0 * (j + 1) / (5);
			x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
			y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
			z = y * sin(phi);
			y = y * cos(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}

	}
	else if (this->type == "G_Luch")
	{
		// ��� G - ����� ����� ������������ ��������� �����������, �� ��� ����� ������ � ��������� ������������
		// ��������� ���� ����� (��� ������� ��������� �����)
		// �������� ��� ������� ���� ������� ������� ��������� ���������� (��� ��� ����� ������� �� ���� phi �����������)
		auto yz1 = this->Yzels_opor[1];
		auto yz2 = this->Yzels_opor[2];
		auto yz3 = this->Yzels_opor[0];
		double dd1 = this->geo->dd1;
		double dd2 = this->geo->dd2;

		int M1 = this->geo->M1;
		int M11 = this->geo->M11;
		int MF = this->geo->MF;
		int M2 = this->geo->M2;

		if (this->parameters.find("phi") == this->parameters.end()) cout << "0989878874  ERROR" << endl;

		double phi = this->parameters["phi"];
		int num = 0;

		// ������ �������� ����
		double x0, y0, x1, y1, t1, t2, tt1, tt2;

		double a1, b1, c1, d1, a2, b2, c2, d2, x, y, z;

		t1 = (yz1->coord[i_time][0] - yz3->coord[i_time][0]) * dd1;
		tt1 = (yz1->func_Ryz(i_time) - yz3->func_Ryz(i_time)) * dd1;
		t2 = 0.0;
		tt2 = 1.0 * dd2;
		x0 = yz1->coord[i_time][0];       
		y0 = yz1->func_Ryz(i_time);
		x1 = yz2->coord[i_time][0];
		y1 = yz2->func_Ryz(i_time);

		//cout << "GGG  " << yz1->number << " " << x0 << " " << y0 << " " << t1 << " " << tt1 << endl;

		a1 = x0;
		b1 = t1;
		c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
		d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

		a2 = y0;
		b2 = tt1;
		c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
		d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

		// ����� �� R2 (TS) �� R3 (HP)
		num = 1;
		for (int j = 0; j < M2 - M11 - MF + 1; j++)
		{
			double s = 1.0 * (j + 1.0) / (M2 - M11 - MF + 1 + 1.0);
			x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
			y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
			z = y * sin(phi);
			y = y * cos(phi);
			this->Yzels[num + j]->coord[i_time][0] = x;
			this->Yzels[num + j]->coord[i_time][1] = y;
			this->Yzels[num + j]->coord[i_time][2] = z;
		}

	}
	
}


Yzel* Luch::get_yzel_near_opor(int num_opor, int shift)
{
	// ��� �������� �� ����� �� ������� !!!
	auto A = this->Yzels_opor[num_opor];

	// ������� ������� �������� ��������
	vector<Yzel*>::iterator target_it = std::find(this->Yzels.begin(), this->Yzels.end(), A);

	// ��������
	if (target_it == this->Yzels.end()) 
	{
		cout << "1546767586 �������� �� ������� � �������" << endl;
		return nullptr;
	}

	if (target_it == this->Yzels.begin()) 
	{
		cout << "9867656565 ��� ��������� ����� ������� ���������" << endl;
		return nullptr;
	}

	return *(target_it + shift);
}
