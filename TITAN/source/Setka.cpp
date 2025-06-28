#include "Setka.h"
#include <algorithm>

#define  macros1(n, x, y, z) A->yzels[n]->coord[0][0] = x; \
	A->yzels[n]->coord[0][1] = y; \
	A->yzels[n]->coord[0][2] = z


// ���� ������ ��������� Gran_TS   Gran_HP    Gran_BS
// ��� ����������� ����� ����, ��� ����� ������������ �� Type_yzel::
#define  macros2(name) for (auto& i : this->All_Gran)\
	{ \
		b1 = true; \
		for (auto& j : i->yzels) \
		{ \
			if (j->type != Type_yzel::name) \
			{ \
				b1 = false; \
				break; \
			} \
		} \
		if (b1 == true) \
		{ \
			this->Gran_##name.push_back(i); \
		} \
	}


#define  macros3(name, proc) if ((100.0 - d2 * 100.0 / d1) > proc) \
		{\
			k++;\
			izmen = true;\
			if (i->parameters.find(#name) != i->parameters.end())\
			{\
				i->parameters[#name] *= (1.0 + procent / 100.0);\
			}\
			else\
			{\
				i->parameters[#name] = this->geo->name * (1.0 + procent / 100.0);\
			}\
		}\
		else if ((100.0 - d2 * 100.0 / d1) < -proc)\
		{\
			k++;\
			izmen = true;\
			if (i->parameters.find(#name) != i->parameters.end())\
			{\
				i->parameters[#name] *= (1.0 - procent / 100.0);\
			}\
			else\
			{\
				i->parameters[#name] = this->geo->name * (1.0 - procent / 100.0);\
			}\
		}

// ���� ������ ��������� �����-������ � Gran_TS   Gran_HP    Gran_BS
#define  macros4(name) for (auto& i : this->Gran_##name)\
	{\
	auto A = i->cells[0];\
	for (auto& j : A->grans)\
	{\
		if (j->number == i->number) continue;\
		auto B = A->Get_Sosed(j);\
		for (auto& k : B->grans)\
		{\
			if (k->type2 == Type_Gran_surf::name && k->number != i->number)\
			{\
				i->grans_surf.push_back(k);\
			}\
		}\
	}\
	}

using namespace std;

Setka::Setka()
{
	this->Surf1 = nullptr;
	this->geo = new Geo_param();
	this->phys_param = new Phys_param();
	Luch::geo = this->geo;
	this->geo->Nphi = 60;     // &INIT&
	// ! ��� ����� ������� �� ����� ������������ ����� (������� ��� �������� �� ���� ��� �������������)
	//cout << "AAAAAAAAAAAAAAAAAAAAA " << endl;
	this->All_name_luch["A_Luch"] = &this->A_Luch;
	this->All_name_luch["B_Luch"] = &this->B_Luch;
	this->All_name_luch["C_Luch"] = &this->C_Luch;
	this->All_name_luch["D_Luch"] = &this->D_Luch;
	this->All_name_luch["E_Luch"] = &this->E_Luch;
	this->All_name_luch["H_Luch"] = &this->H_Luch;
	this->All_name_luch["G_Luch"] = &this->G_Luch;
	this->name_luch.reserve(7);
	this->name_luch.push_back("A_Luch");
	this->name_luch.push_back("B_Luch");
	this->name_luch.push_back("C_Luch");
	this->name_luch.push_back("D_Luch");
	this->name_luch.push_back("E_Luch");
	this->name_luch.push_back("H_Luch");
	this->name_luch.push_back("G_Luch");

	this->Cell_Center = new Cell;

	this->New_initial();                     // ��������� �������� ����� � �����
	this->New_connect();                     // ��������� �������� ������ � ���������� �� � ��������
	// ����� ��������� ������� ��� ������ �����

	this->New_append_surfaces();
	// �������� ������ �� TS, HP, BS 
	// ���������� ��������� ������ � ��������

	this->Renumerate();

	// ��������� ���������� ����� �� ������ ��������� ����
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}

	this->Tecplot_print_all_lush_in_2D();

}

Setka::~Setka()
{
	for (auto& i : this->All_Yzel)
	{
		delete i;
	}
	this->All_Yzel.clear();
	cout << "A1" << endl;
	for (auto& i : this->All_Cell)
	{
		delete i;
	}
	this->All_Cell.clear();

	cout << "A2" << endl;
	for (auto& i : this->All_Gran)
	{
		delete i;
	}
	this->All_Gran.clear();

	cout << "A3" << endl;
	for (auto& i : this->All_Luch)
	{
		delete i;
	}
	this->All_Luch.clear();

	cout << "A4" << endl;
	delete this->Cell_Center;
	cout << "A5" << endl;
	delete this->geo;
	cout << "A6" << endl;

	//delete this->phys_param;
	cout << "A7" << endl;

}

void Setka::Algoritm(short int alg)
{
	// 2 - �����-�����

	cout << "Start Algoritm " << alg << endl;

	if (alg == 2)
	{
		// ��������� ���� ��� ��
		this->Set_MK_Zone();

		//�������� ����
		/*this->Tecplot_print_gran_with_condition(0);
		this->Tecplot_print_gran_with_condition(1);
		this->Tecplot_print_gran_with_condition(2);
		this->Tecplot_print_gran_with_condition(3);
		this->Tecplot_print_gran_with_condition(4);
		this->Tecplot_print_gran_with_condition(5);
		this->Tecplot_print_gran_with_condition(6);*/

		// ������� MK, �������, ������� �������������
		this->MK_prepare(7);
		this->MK_go(7);

		this->MK_delete(7);
	}
}

Cell* Setka::Find_cell_point(const double& x, const double& y, const double& z, short int now, Cell*& previos)
{
	// now - ����� ���������� ������ ��������� (0 �� ���������)

	if (previos == nullptr)
	{
		previos = this->All_Cell[50];
	}

	Cell* next;
	Cell* next2;
	bool b, b2;
	unsigned int kkk = 0;

	while (true)
	{
		kkk++;
		if (kkk > 1000000)
		{
			//cout << "Ne naydeno! A " << endl;
			return nullptr;
		}


		b = previos->Belong_point(x, y, z, now, true, next);

		if (b == true || next == nullptr)
		{
			next2 = next;
			b2 = previos->Belong_point(x, y, z, now, false, next);
			if (b2 == true)
			{
				return previos;
			}
			else
			{
				for (auto& gr : previos->grans)
				{
					auto C = previos->Get_Sosed(gr);
					if (C == nullptr) continue;
					b2 = C->Belong_point(x, y, z, now, false, next);
					if(b2 == true) return C;
				}

				if(next2 == nullptr) cout << "Ne naydeno! B1    " << x << " " << y << " " << z << endl;
				if(next2 != nullptr) cout << "Ne naydeno! B2    " << endl;

				return nullptr;
			}
		}
		else
		{
			previos = next;
			continue;
		}

	}
}

void Setka::Renumerate(void)
{
	int kkk = 1;
	for (auto& i : this->All_Yzel)
	{
		i->number = kkk;
		kkk++;
	}

	kkk = 1;
	for (auto& i : this->All_Cell)
	{
		i->number = kkk;
		kkk++;
	}

	kkk = 1;
	for (auto& i : this->All_Gran)
	{
		i->number = kkk;
		kkk++;
	}
}

bool Setka::Test_geometr(void)
{
	cout << endl;
	cout << "Testing" << endl;

	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "Test 1:" << endl; // ����� ��������� ���� �� 1, ���� �� 2 ������
	// ������� � ����� �� ������ ������ �� ������
	for (auto& i : this->All_Gran)
	{
		if (i->cells.size() != 1 && i->cells.size() != 2)
		{
			cout << "Failure: i->cells.size()  = " << i->cells.size() << endl;
			exit(-1);
		}

		if (i->cells.size() == 2)
		{
			double n1 = i->normal[0][0];
			double n2 = i->normal[0][1];
			double n3 = i->normal[0][2];

			double a1 = i->cells[1]->center[0][0] - i->cells[0]->center[0][0];
			double a2 = i->cells[1]->center[0][1] - i->cells[0]->center[0][1];
			double a3 = i->cells[1]->center[0][2] - i->cells[0]->center[0][2];

			if (scalarProductFast(n1, n2, n3, a1, a2, a3) < 0.0)
			{
				cout << "Failure: 5678675634   = " << scalarProductFast(n1, n2, n3, a1, a2, a3) << endl;
				exit(-1);
			}
		}


	}
	cout << "Success" << endl;


	cout << "Test 2:" << endl;  // � ������ 6 ������ � 8 �����
	for (auto& i : this->All_Cell)
	{
		if (i->grans.size() != 6)
		{
			cout << "Failure: i->grans.size()  = " << i->grans.size() << endl;
			exit(-1);
		}

		if (i->yzels.size() != 8)
		{
			cout << "Failure: i->yzels.size()  = " << i->yzels.size() << endl;
			exit(-1);
		}
	}
	cout << "Success" << endl;

	cout << "Test 3:" << endl;  // ����� ������� ����� �� ���� � (�2) ����� ����������
	int nnn = this->A_Luch[0][0]->Yzels_opor.size();
	int nnn2 = this->C_Luch[0][0]->Yzels_opor.size();
	for (auto& i : this->All_Luch)
	{
		if (i->type == "A_Luch" || i->type == "A2_Luch")
		{
			if (i->Yzels_opor.size() != nnn)
			{
				cout << "Failure: i->Yzels_opor.size()  = " << i->Yzels_opor.size() << "   " << nnn << 
					"   " << i->type << endl;
				exit(-1);
			}
		}

		if (i->type == "C_Luch" || i->type == "C2_Luch")
		{
			if (i->Yzels_opor.size() != nnn2)
			{
				cout << "Failure: i->Yzels_opor.size()  = " << i->Yzels_opor.size() << "   " << nnn2 <<
					"   " << i->type << endl;
				exit(-1);
			}
		}
	}
	cout << "Success" << endl;


	cout << "Test 5:" << endl; // ��� ������� ����� ��������� � � ������� �������� �������������� ����� ����� ��������
	for (auto& i : this->All_Cell)
	{
		double V1 = 0.0;
		double V2 = 0.0;
		double V3 = 0.0;

		for (auto& j : i->grans)
		{
			double n1 = j->normal[0][0];
			double n2 = j->normal[0][1];
			double n3 = j->normal[0][2];
			if (fabs(norm2(n1, n2, n3) - 1.0) > 0.00001)
			{
				cout << "Error 5654354876" << endl;
				whach(norm2(n1, n2, n3));
				exit(-1);
			}

			if (j->cells[0] != i)
			{
				n1 *= -1;
				n2 *= -1;
				n3 *= -1;
			}
			V1 += n1 * j->area[0];
			V2 += n2 * j->area[0];
			V3 += n3 * j->area[0];
		}

		/*if (norm2(V1, V2, V3)/i->volume[0] > 0.1)
		{
			cout << "Error 8768978654" << endl;
			whach(V1);
			whach(V2);
			whach(V3);
			whach(i->volume[0]);
			whach(norm2(V1, V2, V3) / i->volume[0]);
			exit(-1);
		}*/
	}
	cout << "Success" << endl << endl;

	return true;
}

void Setka::Set_luch_parametr()
{
	// ��������� ����� ����������� ��������� 
	// (��������, ���� the � phi ��� ���������� ����� A, B, C �����) - ��� ������� �������, ��� ��� 
	// ��� �������� ����������� ��� ����, � ����� ��������� �� ������������� ��������
	for (auto& i : this->A2_Luch)
	{
		auto A = i->Yzels[0];
		double x = A->coord[0][0];
		double y = A->coord[0][1];
		double z = A->coord[0][2];
		double r = std::sqrt(x * x + y * y + z * z);
		i->parameters["the"] = std::acos(x / r);
		i->parameters["phi"] = polar_angle(y, z);
	}

	for (auto& j : this->A_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->B_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->D_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->E_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->H_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->G_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& i : this->C2_Luch)
	{
		auto A = i->Yzels[0];
		double x = A->coord[0][0];
		double y = A->coord[0][1];
		double z = A->coord[0][2];
		double r = std::sqrt(x * x + y * y + z * z);
		i->parameters["the"] = std::acos(x / r);
		i->parameters["phi"] = polar_angle(y, z);
	}

	for (auto& j : this->C_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}
}

void Setka::Read_old_surface(string name)
{
	cout << "Start Read_old_surface: " << name << endl;
	if (this->Surf1 != nullptr)
	{
		this->Surf1->~Surfaces();
		delete this->Surf1;
	}

	this->Surf1 = new Surfaces();
	this->Surf1->Read_old(name);
	cout << "End Read_old_surface: " << name << endl;
}

void Setka::Move_to_surf(Surfaces* Surf)
{
	double x, y, z, phi, the, r, rr;
	cout << "Start: Move_to_surf" << endl;
	double R_BS; // ��������� ������� ������� ����� (��� �������� B E D �����)
	for (int st = 0; st < 15; st++)
	{
		for (auto& i : this->A_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);

				

				rr = Surf->Get_TS(phi, the);

				j->Yzels_opor[1]->coord[0][0] *= rr / r;
				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "0989898653   errjr" << endl;
				}


				x = j->Yzels_opor[2]->coord[0][0];
				y = j->Yzels_opor[2]->coord[0][1];
				z = j->Yzels_opor[2]->coord[0][2];
				r = sqrt(kvv(x, y, z));

				rr = Surf->Get_HP(phi, the, 0);

				j->Yzels_opor[2]->coord[0][0] *= rr / r;
				j->Yzels_opor[2]->coord[0][1] *= rr / r;
				j->Yzels_opor[2]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "9443563295   errjr" << endl;
				}

				x = j->Yzels_opor[3]->coord[0][0];
				y = j->Yzels_opor[3]->coord[0][1];
				z = j->Yzels_opor[3]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);


				rr = Surf->Get_BS(phi, the);

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "5794671565   errjr" << endl;
				}

				j->Yzels_opor[3]->coord[0][0] *= rr / r;
				j->Yzels_opor[3]->coord[0][1] *= rr / r;
				j->Yzels_opor[3]->coord[0][2] *= rr / r;

			}
		}

		for (auto& j : this->A2_Luch)
		{
			x = j->Yzels_opor[1]->coord[0][0];
			y = j->Yzels_opor[1]->coord[0][1];
			z = j->Yzels_opor[1]->coord[0][2];
			r = sqrt(kvv(x, y, z));
			the = polar_angle(x, sqrt(kv(y) + kv(z)));
			phi = polar_angle(y, z);



			rr = Surf->Get_TS(phi, the);

			j->Yzels_opor[1]->coord[0][0] *= rr / r;
			j->Yzels_opor[1]->coord[0][1] *= rr / r;
			j->Yzels_opor[1]->coord[0][2] *= rr / r;

			if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
			{
				cout << "0989898653   errjr" << endl;
			}


			x = j->Yzels_opor[2]->coord[0][0];
			y = j->Yzels_opor[2]->coord[0][1];
			z = j->Yzels_opor[2]->coord[0][2];
			r = sqrt(kvv(x, y, z));

			rr = Surf->Get_HP(phi, the, 0);

			j->Yzels_opor[2]->coord[0][0] *= rr / r;
			j->Yzels_opor[2]->coord[0][1] *= rr / r;
			j->Yzels_opor[2]->coord[0][2] *= rr / r;

			if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
			{
				cout << "9443563295   errjr" << endl;
			}

			x = j->Yzels_opor[3]->coord[0][0];
			y = j->Yzels_opor[3]->coord[0][1];
			z = j->Yzels_opor[3]->coord[0][2];
			r = sqrt(kvv(x, y, z));
			the = polar_angle(x, sqrt(kv(y) + kv(z)));
			phi = polar_angle(y, z);


			rr = Surf->Get_BS(phi, the);

			if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
			{
				cout << "5794671565   errjr" << endl;
			}

			j->Yzels_opor[3]->coord[0][0] *= rr / r;
			j->Yzels_opor[3]->coord[0][1] *= rr / r;
			j->Yzels_opor[3]->coord[0][2] *= rr / r;

		}

		for (auto& i : this->B_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);

				rr = Surf->Get_TS(phi, the);

				j->Yzels_opor[1]->coord[0][0] *= rr / r;
				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "0989898653   errjr" << endl;
				}


				x = j->Yzels_opor[2]->coord[0][0];
				y = j->Yzels_opor[2]->coord[0][1];
				z = j->Yzels_opor[2]->coord[0][2];

				rr = Surf->Get_HP(phi, x, 1);
				r = sqrt(kvv(0.0, y, z));
				j->Yzels_opor[2]->coord[0][1] *= rr / r;
				j->Yzels_opor[2]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "6510292073   errjr" << endl;
				}

			}
		}

		for (auto& i : this->C_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);



				rr = Surf->Get_TS(phi, the);

				j->Yzels_opor[1]->coord[0][0] *= rr / r;
				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "0989898653   errjr" << endl;
				}
			}
		}

		for (auto& j : this->C2_Luch)
		{
			x = j->Yzels_opor[1]->coord[0][0];
			y = j->Yzels_opor[1]->coord[0][1];
			z = j->Yzels_opor[1]->coord[0][2];
			r = sqrt(kvv(x, y, z));
			the = polar_angle(x, sqrt(kv(y) + kv(z)));
			phi = polar_angle(y, z);



			rr = Surf->Get_TS(phi, the);

			j->Yzels_opor[1]->coord[0][0] *= rr / r;
			j->Yzels_opor[1]->coord[0][1] *= rr / r;
			j->Yzels_opor[1]->coord[0][2] *= rr / r;

			if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
			{
				cout << "0989898653   errjr" << endl;
			}
		}

		for (auto& i : this->D_Luch)
		{
			x = i[0]->Yzels_opor[1]->coord[0][0];
			y = i[0]->Yzels_opor[1]->coord[0][1];
			z = i[0]->Yzels_opor[1]->coord[0][2];
			r = sqrt(kvv(0.0, y, z));
			phi = polar_angle(y, z);
			rr = Surf->Get_HP(phi, x, 1);

			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				//phi = polar_angle(y, z);
				//rr = Surf->Get_HP(phi, x, 1);

				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

			}
		}

		for (auto& i : this->E_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				phi = polar_angle(y, z);
				rr = Surf->Get_HP(phi, x, 1);

				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "6510292073   errjr" << endl;
				}
			}
		}
	
		// ������� BS ��� B E D �����
		for (int i = 0; i < this->B_Luch.size(); i++)
		{
			R_BS = this->A_Luch[i].back()->Yzels_opor[3]->func_R(0);
			for (auto& j : this->B_Luch[i])
			{
				y = j->Yzels_opor[3]->coord[0][1];
				z = j->Yzels_opor[3]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));

				j->Yzels_opor[3]->coord[0][1] *= R_BS / r;
				j->Yzels_opor[3]->coord[0][2] *= R_BS / r;
			}

			for (auto& j : this->E_Luch[i])
			{
				y = j->Yzels_opor[2]->coord[0][1];
				z = j->Yzels_opor[2]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));

				j->Yzels_opor[2]->coord[0][1] *= R_BS / r;
				j->Yzels_opor[2]->coord[0][2] *= R_BS / r;
			}

			for (auto& j : this->D_Luch[i])
			{
				y = j->Yzels_opor[2]->coord[0][1];
				z = j->Yzels_opor[2]->coord[0][2];
				r = sqrt(kvv(0.0, y, z));

				j->Yzels_opor[2]->coord[0][1] *= R_BS / r;
				j->Yzels_opor[2]->coord[0][2] *= R_BS / r;
			}
		}

		for (auto& i : this->All_Luch)
		{
			i->dvigenie(0);
		}

	}

	cout << "End: Move_to_surf" << endl;
}

void Setka::New_initial()
{
	cout << "---START New_initial---" << endl;

	int n, i, j, n2, n3;
	double x, y;
	vector<Yzel*> Yz_;

	// ��������� ���� 2� �����
	ifstream ffin("SDK1_2D_Setka.bin", ios::binary | ios::in);
	if (!ffin)
	{
		cout << "Net takogo fajla (fajl 2D setki)" << endl;
		exit(-1);
	}

	ffin.read((char*)&n, sizeof n);  // ������ �����
	cout << "Versiya fajla 2D setki: " << n << endl;

	// ��������� �������������� ���������
	if (true)
	{
		ffin.read((char*)&x, sizeof x);  this->geo->tetta0 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->tetta1 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->tetta2 = x;

		ffin.read((char*)&n, sizeof n);  this->geo->N1 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->N2 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->N3 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->N4 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->N5 = n;

		ffin.read((char*)&n, sizeof n);  this->geo->M0 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M1 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M11 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M2 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M3 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M4 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->MF = n;

		ffin.read((char*)&x, sizeof x);  this->geo->R0 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R1 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R2 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R3 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R4 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R5 = x;

		ffin.read((char*)&x, sizeof x);  this->geo->L6 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->L7 = x;

		ffin.read((char*)&x, sizeof x);  this->geo->dd1 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd2 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd3 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd4 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd5 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd6 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd7 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd8 = x;
	}

	// ������ �������������� ���������, ������� �����
	this->geo->dd3 = 19.0;
	this->geo->dd4 = 15.0;
	this->geo->dd5 = 19.0;
	this->geo->dd6 = 14.0;


	this->geo->dd7 = 0.7;


	// ��������� ����
	int sdvig_yz = 0;
	if (true)
	{
		ffin.read((char*)&n, sizeof n); // ����� �����
		for (int i = 0; i < n; i++)
		{
			ffin.read((char*)&x, sizeof x);
			ffin.read((char*)&y, sizeof y);
			auto yzel = new Yzel(x, y, 0.0);  // ������ ����
			this->All_Yzel.push_back(yzel);
			Yz_.push_back(yzel);
			sdvig_yz++;
		}

		Yzel_2D.push_back(Yz_);  // � ���� ������ ���������� ����������� ���������. ������ �������� ������������, �.�. ����� ������ Yz_, ��� �� �������� �� Yzel_2D
		Yz_.clear();
	}

	// ��������� ����
	if (true)
	{
		vector<Luch*> Luch_;
		//for (auto& [name, vec_ptr] : this->All_name_luch)
		for(auto& name : this->name_luch)
		{
			auto vec_ptr = this->All_name_luch[name];
			cout << "Zopolnyaem  " << name << endl;
			ffin.read((char*)&n, sizeof n); // ����� �����
			for (i = 0; i < n; i++)
			{
				auto Luch_0 = new Luch();  // ������ ���
				All_Luch.push_back(Luch_0);
				Luch_0->type = name;
				Luch_.push_back(Luch_0);
				ffin.read((char*)&n2, sizeof n2); // ����� �����
				for (j = 0; j < n2; j++)
				{
					ffin.read((char*)&n3, sizeof n3); // ����� ����
					if (Yzel_2D[0].size() < static_cast<size_t>(n3))
					{
						cout << "86740576536 PROBLEM   " << Yzel_2D[0].size() << " " << n3;
					}
					Luch_0->Yzels.push_back(Yzel_2D[0][n3 - 1]);
				}

				ffin.read((char*)&n2, sizeof n2); // ����� ������� �����
				//cout << name << "  Opor yzel = " << n2 << endl;
				for (j = 0; j < n2; j++)
				{
					ffin.read((char*)&n3, sizeof n3); // ����� ����
					Luch_0->Yzels_opor.push_back(Yzel_2D[0][n3 - 1]);
				}
			}
			vec_ptr->push_back(Luch_);
			cout << "Dobavili luchey:  " << Luch_.size() << " " << n << endl;
			Luch_.clear();
		}
	}

	// ������ ���� ������� ����� � 3�, ���������� ����� ����� � ����� � ���

	// ����������� ����� �����
	if (true)
	{
		for (int iii1 = 1; iii1 < this->geo->Nphi; iii1++)
		{
			Yz_.clear();
			for (auto& i : Yzel_2D[0])
			{
				auto yzel = new Yzel(i->coord[0][0], i->coord[0][1], i->coord[0][2]);  // ������ ����
				Yz_.push_back(yzel);
				this->All_Yzel.push_back(yzel);
			}
			this->Yzel_2D.push_back(Yz_);

			// ������� ���������� (������������ ���������)
			double tet = iii1 * (2.0 * const_pi) / this->geo->Nphi;
			for (auto& i : Yz_)
			{
				double newy = cos(tet) * i->coord[0][1] - sin(tet) * i->coord[0][2];
				double newz = sin(tet) * i->coord[0][1] + cos(tet) * i->coord[0][2];
				i->coord[0][1] = newy;
				i->coord[0][2] = newz;
			}

			Yz_.clear();
		}
	}

	// ����������� ����� �����
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}

		vector<Luch*> Luch_;
		for (auto& [name, vec_ptr] : this->All_name_luch)
		{
			vec_ptr->resize(this->geo->Nphi);
			cout << "Yvelichivaim luchi = " << name << endl;
			for (int iii1 = 1; iii1 < this->geo->Nphi; iii1++)
			{
				//cout << "cds1 = " << (*vec_ptr).size() << endl;
				//cout << "cds2 = " << (*vec_ptr)[0].size() << endl;
				//cout << "cds3 = " << (*vec_ptr)[1].size() << endl;
				//cout << "cds2 = " << (*vec_ptr)[0][0]->Yzels.size() << endl;

				(*vec_ptr)[iii1].resize((*vec_ptr)[0].size());
				int kk = 0;
				for (auto& ii : (*vec_ptr)[0])
				{
					auto Luch_0 = new Luch();  // ������ ���
					All_Luch.push_back(Luch_0);
					Luch_0->type = name;

					(*vec_ptr)[iii1][kk] = Luch_0;
					kk++;

					for (auto& j : ii->Yzels)
					{
						Luch_0->Yzels.push_back(this->All_Yzel[j->number - 1 + iii1 * (sdvig_yz)]);
					}
					for (auto& j : ii->Yzels_opor)
					{
						Luch_0->Yzels_opor.push_back(this->All_Yzel[j->number - 1 + iii1 * (sdvig_yz)]);
					}
				}
			}
		}
	}

	// ��������� ������
	if (true)
	{
		this->Cell_2D.resize(this->geo->Nphi);

		ffin.read((char*)&n, sizeof n); // ����� �����
		for (int i = 0; i < this->geo->Nphi; i++)
		{
			this->Cell_2D[i].resize(n);
		}
		this->All_Cell.reserve(this->geo->Nphi* n * 2);

		for (int i = 0; i < n; i++) // ��������� ������ ������
		{
			ffin.read((char*)&n2, sizeof n2);
			auto cell_ = new Cell();  // ������ ������
			this->All_Cell.push_back(cell_);
			this->Cell_2D[0][i] = cell_;

			for (int j = 0; j < n2; j++)
			{
				ffin.read((char*)&n3, sizeof n3);
				cell_->yzels.push_back(Yzel_2D[0][n3 - 1]);
			}
		}

		for (int iii1 = 1; iii1 < this->geo->Nphi; iii1++)
		{
			for (int i = 0; i < n; i++) // ��������� ������ ������
			{
				n2 = this->Cell_2D[0][i]->yzels.size();
				auto cell_ = new Cell();  // ������ ������
				this->All_Cell.push_back(cell_);
				this->Cell_2D[iii1][i] = cell_;

				for (int j = 0; j < n2; j++)
				{
					n3 = this->Cell_2D[0][i]->yzels[j]->number;
					cell_->yzels.push_back(Yzel_2D[iii1][n3 - 1]);
				}
			}
		}
	}

	// ������ ���� � ������ ������ �������� ���� � �������� ���� �����
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}

		for (int i = 0; i < this->Cell_2D.size(); i++)
		{
			j = i + 1;
			if (j > this->Cell_2D.size() - 1) j = 0;
			for (int k = 0; k < this->Cell_2D[0].size(); k++)
			{
				for (int kk = 0; kk < 4; kk++)
				{
					int nn = this->Cell_2D[j][k]->yzels[kk]->number;
					this->Cell_2D[i][k]->yzels.push_back(All_Yzel[nn - 1]);
				}
			}
		}
	}

	// ��������� ���� �����
	ifstream fin("SDK1_krug_setka.bin", ios::binary | ios::in);
	if (!fin)
	{
		cout << "Net takogo fajla (fajl setki v krugu)" << endl;
		exit(-1);
	}

	fin.read((char*)&n, sizeof n);  // ������ �����
	cout << "Versiya fajla setki v krugu: " << n << endl;

	// ��������� ���� (��������� �� � ������ ���� �������� �����)
	if (true) 
	{
		this->Krug_Yzel.resize(this->A_Luch[0][0]->Yzels.size());
		this->Krug_Yzel_2.resize(this->C_Luch[0][0]->Yzels.size());

		fin.read((char*)&n, sizeof n); // ����� �����

		for (int i = 0; i < this->Krug_Yzel.size(); i++)
		{
			this->Krug_Yzel[i].resize(n);
		}

		for (int i = 0; i < this->Krug_Yzel_2.size(); i++)
		{
			this->Krug_Yzel_2[i].resize(n);
		}
		
		for (int i = 0; i < n; i++)
		{
			fin.read((char*)&x, sizeof x);
			fin.read((char*)&y, sizeof y);

			for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
			{
				auto yzel = new Yzel(x, y, 0.0);
				//this->All_Yzel.push_back(yzel);
				yzel->number = -1;
				this->Krug_Yzel[ii][i] = yzel;
			}

			for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
			{
				auto yzel = new Yzel(x, y, 0.0);
				yzel->number = -1;
				//this->All_Yzel.push_back(yzel);
				this->Krug_Yzel_2[ii][i] = yzel;
			}
		}

		double Rmax = 1.0;
		double thetamin = 4.0;
		// ������� ������ �����

		//Rmax = sqrt(kv(this->Krug_Yzel[0][0]->coord[0][0]) + kv(this->Krug_Yzel[0][0]->coord[0][1]));
		//cout << "R max: " << Rmax << endl;

		Yzel* AAA = nullptr;
		for (int i = 0; i < this->Krug_Yzel[0].size(); i++)
		{
			if (this->Krug_Yzel[0][i]->func_R(0) > 1.0 - 0.00001)
			{
				AAA = this->Krug_Yzel[0][i];
				break;
			}
		}

		thetamin = polar_angle(AAA->coord[0][0], AAA->coord[0][1]);
		cout << "Ugol vrashcheniya setki: " << thetamin * 180.0 / const_pi << endl;
		cout << "this->geo.tetta0: " << this->geo->tetta0 << endl;
		cout << "this->geo.tetta1: " << this->geo->tetta1 << endl;

		// ������� ���������� ����� 1 �� �����
		for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel[ii])
			{
				i->coord[0][0] /= Rmax;
				i->coord[0][1] /= Rmax;
				i->coord[0][2] = i->coord[0][1];
				i->coord[0][1] = i->coord[0][0];
				i->coord[0][0] = 1.0 / tan(this->geo->tetta0);// 1.0 / tan(this->geo.tetta0 - (const_pi / 2.0 - this->geo.tetta0) / (this->geo.N1 - 1));
				double r = sqrt(kv(i->coord[0][0]) + kv(i->coord[0][1]) + kv(i->coord[0][2]));
				i->coord[0][0] *= this->geo->R0 / r * (ii + 1);
				i->coord[0][1] *= this->geo->R0 / r * (ii + 1);
				i->coord[0][2] *= this->geo->R0 / r * (ii + 1);
			}
		}

		// ������� ���������� ��� ���������� � ������� ������ �����
		for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel[ii])
			{
				double newy = cos(-thetamin) * i->coord[0][1] - sin(-thetamin) * i->coord[0][2];
				double newz = sin(-thetamin) * i->coord[0][1] + cos(-thetamin) * i->coord[0][2];
				i->coord[0][1] = newy;
				i->coord[0][2] = newz;
			}
		}

		// ������� ���������� ����� 2 �� �����
		for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel_2[ii])
			{
				i->coord[0][0] /= Rmax;
				i->coord[0][1] /= Rmax;
				i->coord[0][2] = i->coord[0][1];
				i->coord[0][1] = i->coord[0][0];
				i->coord[0][0] = 1.0 / tan(this->geo->tetta1);// 1.0 / tan(this->geo.tetta0 - (const_pi / 2.0 - this->geo.tetta0) / (this->geo.N1 - 1));
				double r = sqrt(kv(i->coord[0][0]) + kv(i->coord[0][1]) + kv(i->coord[0][2]));
				i->coord[0][0] *= this->geo->R0 / r * (ii + 1);
				i->coord[0][1] *= this->geo->R0 / r * (ii + 1);
				i->coord[0][2] *= this->geo->R0 / r * (ii + 1);

			}
		}

		// ������� ���������� ��� ���������� � ������� ������ �����
		for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel_2[ii])
			{
				double newy = cos(-thetamin) * i->coord[0][1] - sin(-thetamin) * i->coord[0][2];
				double newz = sin(-thetamin) * i->coord[0][1] + cos(-thetamin) * i->coord[0][2];
				i->coord[0][1] = newy;
				i->coord[0][2] = newz;
			}
		}

	}

	// ���������������� ��� ����
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}
	}

	// ���� ������� ������ ���� (������ ���������)
	// ��������� ���� ���� ���������
	// ��� ����� ���� �� �� ������, ���� �� ����� � ����� � ����� ������ � ���������� ������� (������� ���������)
	cout << "Yzlov do dobavleniya  " << this->All_Yzel.size() << endl;
	for (int i = 0; i < this->A_Luch.size(); i++)
	{
		auto A = this->A_Luch[i][0]->Yzels[0];
		for (int j = 0; j < this->Krug_Yzel[0].size(); j++)
		{
			auto B = this->Krug_Yzel[0][j];


			double d = sqrt(kv(A->coord[0][0] - B->coord[0][0]) + kv(A->coord[0][1] - B->coord[0][1]) +
				kv(A->coord[0][2] - B->coord[0][2]));
			if (d < this->geo->R0 / 1000.0)
			{
				//cout << "Nashol " << j << endl;

				for (int k = 0; k < this->Krug_Yzel.size(); k++)
				{
					delete this->Krug_Yzel[k][j];
					this->Krug_Yzel[k][j] = this->A_Luch[i][0]->Yzels[k];
				}

				break;
			}
		}
		
	}

	for (int i = 0; i < this->C_Luch.size(); i++)
	{
		auto A = this->C_Luch[i].back()->Yzels[0];
		for (int j = 0; j < this->Krug_Yzel_2[0].size(); j++)
		{
			auto B = this->Krug_Yzel_2[0][j];


			double d = sqrt(kv(A->coord[0][0] - B->coord[0][0]) + kv(A->coord[0][1] - B->coord[0][1]) +
				kv(A->coord[0][2] - B->coord[0][2]));
			if (d < this->geo->R0 / 1000.0)
			{
				//cout << "Nashol " << j << endl;

				for (int k = 0; k < this->Krug_Yzel_2.size(); k++)
				{
					delete this->Krug_Yzel_2[k][j];
					this->Krug_Yzel_2[k][j] = this->C_Luch[i].back()->Yzels[k];
				}

				break;
			}
		}
	}

	// ��������� ���� � ����� ������ �����
	for (int i = 0; i < this->Krug_Yzel.size(); i++)
	{
		for (int j = 0; j < this->Krug_Yzel[0].size(); j++)
		{
			if (this->Krug_Yzel[i][j]->number < 0)
			{
				this->All_Yzel.push_back(this->Krug_Yzel[i][j]);
			}
		}
	}

	for (int i = 0; i < this->Krug_Yzel_2.size(); i++)
	{
		for (int j = 0; j < this->Krug_Yzel_2[0].size(); j++)
		{
			if (this->Krug_Yzel_2[i][j]->number < 0)
			{
				this->All_Yzel.push_back(this->Krug_Yzel_2[i][j]);
			}
		}
	}

	cout << "Yzlov posle dobavleniya  " << this->All_Yzel.size() << endl;

	// ���������������� ��� ����
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}
	}

	// ���� ������� ���� �� ����� �� ������ (A2 � C2 �����)  // ---------------------------------------------------------------------------------------------
	if (true)
	{
		this->A2_Luch.reserve(this->Krug_Yzel[0].size() - 60);
		this->C2_Luch.reserve(this->Krug_Yzel_2[0].size() - 60); // ����� ������� ����� �� �������, ��� ���
			// ��� ��� ������ � ���� A

		for (int j = 0; j < this->Krug_Yzel[0].size(); j++) 
		{
			int nf = this->Krug_Yzel[0][j]->number;

			bool hasN = std::any_of(  // �������� �� ������ �������, ��������������� �������
				this->A_Luch.begin(),
				this->A_Luch.end(),
				[nf](auto x) { return x[0]->Yzels[0]->number == nf; }
			);

			if (hasN) continue;

			auto Luch_0 = new Luch();  // ������ ���
			All_Luch.push_back(Luch_0);
			Luch_0->type = "A2_Luch";
			this->A2_Luch.push_back(Luch_0);
			for (int i = 0; i < this->Krug_Yzel.size(); i++)
			{
				Luch_0->Yzels.push_back(this->Krug_Yzel[i][j]);
			}


			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 + this->geo->M1 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M2 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M2 + 1 + this->geo->M3 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M2 + 1 + this->geo->M3 + 1 + this->geo->M4 + 1 - 1][j]);
		}


		for (int j = 0; j < this->Krug_Yzel_2[0].size(); j++)
		{
			int nf = this->Krug_Yzel_2[0][j]->number;

			bool hasN = std::any_of(  // �������� �� ������ �������, ��������������� �������
				this->C_Luch.begin(),
				this->C_Luch.end(),
				[nf](auto x) { return x.back()->Yzels[0]->number == nf; }
			);

			if (hasN) continue;

			auto Luch_0 = new Luch();  // ������ ���
			All_Luch.push_back(Luch_0);
			Luch_0->type = "C2_Luch";
			this->C2_Luch.push_back(Luch_0);
			for (int i = 0; i < this->Krug_Yzel_2.size(); i++)
			{
				Luch_0->Yzels.push_back(this->Krug_Yzel_2[i][j]);
			}


			Luch_0->Yzels_opor.push_back(this->Krug_Yzel_2[this->geo->M0 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel_2[this->geo->M0 + 1 + this->geo->M1 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel_2[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M11 + this->geo->N4 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel_2[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M11 + this->geo->N4 + this->geo->N5 - 1][j]);
		}
		
	}


	// ������ ���� ������� ���� � ������������ � ��������� ��������
	if (true) 
	{
		// SET_PARAMETER   &INIT&
		this->geo->R1 = 7.5;
		this->geo->R2 = 20.0;
		this->geo->R3 = 30.0;
		this->geo->R4 = 140.0;
		this->geo->R5 = 400.0;
		this->geo->L6 = -60.0;
		this->geo->L7 = -160.0;
		
		this->Set_luch_parametr();
		// ���������� ������� ����� � �2-����� � ���������� �������
		double x, y, z, r;

		for (auto& i : this->All_Luch)
		{
			if (i->type == "A_Luch")
			{
				auto A = i->Yzels_opor[0];
				x = A->coord[0][0];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(x, y, z));
				A->coord[0][0] = x / r * this->geo->R1;
				A->coord[0][1] = y / r * this->geo->R1;
				A->coord[0][2] = z / r * this->geo->R1;

				A = i->Yzels_opor[1];
				x = A->coord[0][0];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(x, y, z));
				A->coord[0][0] = x / r * this->geo->R2;
				A->coord[0][1] = y / r * this->geo->R2;
				A->coord[0][2] = z / r * this->geo->R2;

				A = i->Yzels_opor[2];
				x = A->coord[0][0];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(x, y, z));
				A->coord[0][0] = x / r * this->geo->R3;
				A->coord[0][1] = y / r * this->geo->R3;
				A->coord[0][2] = z / r * this->geo->R3;

				A = i->Yzels_opor[3];
				x = A->coord[0][0];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(x, y, z));
				A->coord[0][0] = x / r * this->geo->R4;
				A->coord[0][1] = y / r * this->geo->R4;
				A->coord[0][2] = z / r * this->geo->R4;

				A = i->Yzels_opor[4];
				x = A->coord[0][0];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(x, y, z));
				A->coord[0][0] = x / r * this->geo->R5;
				A->coord[0][1] = y / r * this->geo->R5;
				A->coord[0][2] = z / r * this->geo->R5;
			}
			else if (i->type == "C_Luch" || i->type == "B_Luch")
			{
				auto A = i->Yzels_opor[0];
				x = A->coord[0][0];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(x, y, z));
				A->coord[0][0] = x / r * this->geo->R1;
				A->coord[0][1] = y / r * this->geo->R1;
				A->coord[0][2] = z / r * this->geo->R1;

				A = i->Yzels_opor[1];
				x = A->coord[0][0];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(x, y, z));
				A->coord[0][0] = x / r * this->geo->R2;
				A->coord[0][1] = y / r * this->geo->R2;
				A->coord[0][2] = z / r * this->geo->R2;
			}
		}

		for (auto& i : this->All_Luch)
		{
			if (i->type == "B_Luch")
			{
				auto A = i->Yzels_opor[2];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				A->coord[0][0] = i->Yzels_opor[1]->coord[0][0];
				A->coord[0][1] = y / r * this->geo->R3;
				A->coord[0][2] = z / r * this->geo->R3;

				A = i->Yzels_opor[3];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				A->coord[0][0] = i->Yzels_opor[1]->coord[0][0];
				A->coord[0][1] = y / r * this->geo->R4;
				A->coord[0][2] = z / r * this->geo->R4;
			}

			if (i->type == "D_Luch")
			{
				auto A = i->Yzels_opor[1];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				A->coord[0][0] = i->Yzels_opor[0]->coord[0][0];
				A->coord[0][1] = y / r * this->geo->R3;
				A->coord[0][2] = z / r * this->geo->R3;

				A = i->Yzels_opor[2];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				A->coord[0][0] = i->Yzels_opor[0]->coord[0][0];
				A->coord[0][1] = y / r * this->geo->R4;
				A->coord[0][2] = z / r * this->geo->R4;
			}
		}


		for (auto& i : this->A2_Luch)
		{
			auto A = i->Yzels_opor[0];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2]; 
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R1;
			A->coord[0][1] = y / r * this->geo->R1;
			A->coord[0][2] = z / r * this->geo->R1;

			A = i->Yzels_opor[1];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R2;
			A->coord[0][1] = y / r * this->geo->R2;
			A->coord[0][2] = z / r * this->geo->R2;

			A = i->Yzels_opor[2];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R3;
			A->coord[0][1] = y / r * this->geo->R3;
			A->coord[0][2] = z / r * this->geo->R3;

			A = i->Yzels_opor[3];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R4;
			A->coord[0][1] = y / r * this->geo->R4;
			A->coord[0][2] = z / r * this->geo->R4;

			A = i->Yzels_opor[4];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R5;
			A->coord[0][1] = y / r * this->geo->R5;
			A->coord[0][2] = z / r * this->geo->R5;
		}

		for (auto& i : this->C2_Luch)
		{
			auto A = i->Yzels_opor[0];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R1;
			A->coord[0][1] = y / r * this->geo->R1;
			A->coord[0][2] = z / r * this->geo->R1;

			A = i->Yzels_opor[1];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R2;
			A->coord[0][1] = y / r * this->geo->R2;
			A->coord[0][2] = z / r * this->geo->R2;
		}

		// ���� �� ������������� ������ �����, ����� ������� �� � ���������� �������
		std::sort(this->All_Luch.begin(), this->All_Luch.end(), [](Luch* a, Luch* b) {
			if (a->type == "A_Luch" && b->type != "A_Luch")
			{
				return true;
			}

			if (a->type == "B_Luch" && b->type != "A_Luch" && b->type != "B_Luch")
			{
				return true;
			}

			if (a->type == "C_Luch" && b->type != "A_Luch" && 
				b->type != "B_Luch" && b->type != "C_Luch")
			{
				return true;
			}
			if (a->type == "D_Luch" && b->type != "A_Luch" &&
				b->type != "B_Luch" && b->type != "C_Luch" &&
				b->type != "D_Luch")
			{
				return true;
			}

			if (a->type == "E_Luch" && b->type != "A_Luch" &&
				b->type != "B_Luch" && b->type != "C_Luch" && 
				b->type != "D_Luch" && b->type != "E_Luch")
			{
			return true;
			}

			if (a->type == "H_Luch" && b->type != "A_Luch" &&
				b->type != "B_Luch" && b->type != "C_Luch" &&
				b->type != "D_Luch" && b->type != "E_Luch" &&
				b->type != "H_Luch")
			{
			return true;
			}

			return false;
			});

		// ������� � E-���� ������� ������ �� C-����
		for (int i = 0; i < this->E_Luch.size(); i++)
		{
			for (int j = 0; j < this->E_Luch[0].size(); j++)
			{
				this->E_Luch[i][j]->Yzels_opor.insert(this->E_Luch[i][j]->Yzels_opor.begin(), 
					this->C_Luch[i][0]->Yzels[this->geo->M0 + this->geo->M1 + this->geo->M11 + 1 + j]);
			}
		}


		for (auto& i : this->All_Luch)
		{
			if (i->type == "E_Luch")
			{
				auto A = i->Yzels_opor[1];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				A->coord[0][0] = i->Yzels_opor[0]->coord[0][0];
				A->coord[0][1] = y / r * this->geo->R3;
				A->coord[0][2] = z / r * this->geo->R3;

				A = i->Yzels_opor[2];
				y = A->coord[0][1];
				z = A->coord[0][2];
				r = sqrt(kvv(0.0, y, z));
				A->coord[0][0] = i->Yzels_opor[0]->coord[0][0];
				A->coord[0][1] = y / r * this->geo->R4;
				A->coord[0][2] = z / r * this->geo->R4;
			}
		}

		for (auto& i : this->All_Luch)
		{
			//cout << i->type << endl;
			i->dvigenie(0);
		}
		//for (auto& nn : this->name_luch) // ���-�� ���������� ���� � ���������� �������
		//{
		//	auto A = *(All_name_luch[nn]);
		//	for (auto& j : A)
		//	{
		//		for (auto& i : j)
		//		{
		//			i->dvigenie(0);
		//		}
		//	}
		//}
	}

	// ��������� ������ (��� ����� � ������)
	if (true)
	{
		this->Cell_layer_head.resize(this->Krug_Yzel.size() - 1);
		this->Cell_layer_tail.resize(this->Krug_Yzel_2.size() - 1);

		int nn, ny;
		fin.read((char*)&n, sizeof n); // ����� �����

		for (int i = 0; i < this->Cell_layer_head.size(); i++)
		{
			this->Cell_layer_head[i].resize(n);
		}

		for (int i = 0; i < this->Cell_layer_tail.size(); i++)
		{
			this->Cell_layer_tail[i].resize(n);
		}

		for (int i = 0; i < this->Cell_layer_head.size(); i++)
		{
			for (int j = 0; j < this->Cell_layer_head[0].size(); j++)
			{
				auto cell = new Cell();
				this->Cell_layer_head[i][j] = cell;
				this->All_Cell.push_back(cell);
			}
		}

		for (int i = 0; i < this->Cell_layer_tail.size(); i++)
		{
			for (int j = 0; j < this->Cell_layer_tail[0].size(); j++)
			{
				auto cell = new Cell();
				this->Cell_layer_tail[i][j] = cell;
				this->All_Cell.push_back(cell);
			}
		}


		for (int i = 0; i < n; i++)
		{
			fin.read((char*)&nn, sizeof nn); // ����� ����� � ������
			for (int j = 0; j < nn; j++)
			{
				fin.read((char*)&ny, sizeof ny); // ����� ����

				for (int k = 0; k < this->Cell_layer_head.size(); k++)
				{
					this->Cell_layer_head[k][i]->yzels.push_back(Krug_Yzel[k][ny - 1]);
					this->Cell_layer_head[k][i]->yzels.push_back(Krug_Yzel[k + 1][ny - 1]);
				}

				for (int k = 0; k < this->Cell_layer_tail.size(); k++)
				{
					this->Cell_layer_tail[k][i]->yzels.push_back(Krug_Yzel_2[k][ny - 1]);
					this->Cell_layer_tail[k][i]->yzels.push_back(Krug_Yzel_2[k + 1][ny - 1]);
				}
			}
			
		}

		// ������ ������� ����, ����� ��� ��� � ���������� �������
		if (true)
		{
			for (int i = 0; i < this->Cell_layer_head.size(); i++)
			{
				for (int j = 0; j < this->Cell_layer_head[0].size(); j++)
				{
					auto a1 = this->Cell_layer_head[i][j]->yzels[0];
					auto a2 = this->Cell_layer_head[i][j]->yzels[1];
					auto a3 = this->Cell_layer_head[i][j]->yzels[2];
					auto a4 = this->Cell_layer_head[i][j]->yzels[3];
					auto a5 = this->Cell_layer_head[i][j]->yzels[4];
					auto a6 = this->Cell_layer_head[i][j]->yzels[5];
					auto a7 = this->Cell_layer_head[i][j]->yzels[6];
					auto a8 = this->Cell_layer_head[i][j]->yzels[7];

					this->Cell_layer_head[i][j]->yzels[1] = a3;
					this->Cell_layer_head[i][j]->yzels[2] = a5;
					this->Cell_layer_head[i][j]->yzels[3] = a7;
					this->Cell_layer_head[i][j]->yzels[4] = a2;
					this->Cell_layer_head[i][j]->yzels[5] = a4;
					this->Cell_layer_head[i][j]->yzels[6] = a6;
					this->Cell_layer_head[i][j]->yzels[7] = a8;
				}
			}

			for (int i = 0; i < this->Cell_layer_tail.size(); i++)
			{
				for (int j = 0; j < this->Cell_layer_tail[0].size(); j++)
				{
					auto a1 = this->Cell_layer_tail[i][j]->yzels[0];
					auto a2 = this->Cell_layer_tail[i][j]->yzels[1];
					auto a3 = this->Cell_layer_tail[i][j]->yzels[2];
					auto a4 = this->Cell_layer_tail[i][j]->yzels[3];
					auto a5 = this->Cell_layer_tail[i][j]->yzels[4];
					auto a6 = this->Cell_layer_tail[i][j]->yzels[5];
					auto a7 = this->Cell_layer_tail[i][j]->yzels[6];
					auto a8 = this->Cell_layer_tail[i][j]->yzels[7];

					this->Cell_layer_tail[i][j]->yzels[1] = a3;
					this->Cell_layer_tail[i][j]->yzels[2] = a5;
					this->Cell_layer_tail[i][j]->yzels[3] = a7;
					this->Cell_layer_tail[i][j]->yzels[4] = a2;
					this->Cell_layer_tail[i][j]->yzels[5] = a4;
					this->Cell_layer_tail[i][j]->yzels[6] = a6;
					this->Cell_layer_tail[i][j]->yzels[7] = a8;
				}
			}
		}

		
	}

	// �� ��������� ������ ������ �������, ������ ���������� �� ������� (�������� �������)!
	// ������� �����, ����, ������� ��
	// ������� ����������� ������ �� �����

	fin.close();
	ffin.close();
	cout << "---END New_initial---" << endl;
}

void Setka::New_connect()
{
	// ���� ���������: ������� ������ ��� ������ ������ 6 ������
	// ����� ������� �������, ������� ������ � ��������� ����� � ��������
	// ��� ��������� �������� � �.�. ������� �� ������ � �������� ��������� ������, ����� �������� �� � ������

	cout << "START: New_connect" << endl;
	vector<Gran*> All_gran_;
	std::unordered_set<int> number_gran_for_delete;

	this->Renumerate();

	// ������ ��� �����
	for (auto& i : this->All_Cell)
	{
		auto G = new Gran();
		G->yzels.push_back(i->yzels[0]);
		G->yzels.push_back(i->yzels[3]);
		G->yzels.push_back(i->yzels[2]);
		G->yzels.push_back(i->yzels[1]);
		All_gran_.push_back(G);
		G->cells.push_back(i);

		G = new Gran();
		G->yzels.push_back(i->yzels[4]);
		G->yzels.push_back(i->yzels[5]);
		G->yzels.push_back(i->yzels[6]);
		G->yzels.push_back(i->yzels[7]);
		All_gran_.push_back(G);
		G->cells.push_back(i);

		G = new Gran();
		G->yzels.push_back(i->yzels[0]);
		G->yzels.push_back(i->yzels[1]);
		G->yzels.push_back(i->yzels[5]);
		G->yzels.push_back(i->yzels[4]);
		All_gran_.push_back(G);
		G->cells.push_back(i);

		G = new Gran();
		G->yzels.push_back(i->yzels[2]);
		G->yzels.push_back(i->yzels[3]);
		G->yzels.push_back(i->yzels[7]);
		G->yzels.push_back(i->yzels[6]);
		All_gran_.push_back(G);
		G->cells.push_back(i);

		G = new Gran();
		G->yzels.push_back(i->yzels[1]);
		G->yzels.push_back(i->yzels[2]);
		G->yzels.push_back(i->yzels[6]);
		G->yzels.push_back(i->yzels[5]);
		All_gran_.push_back(G);
		G->cells.push_back(i);

		G = new Gran();
		G->yzels.push_back(i->yzels[0]);
		G->yzels.push_back(i->yzels[4]);
		G->yzels.push_back(i->yzels[7]);
		G->yzels.push_back(i->yzels[3]);
		All_gran_.push_back(G);
		G->cells.push_back(i);
	}

	// ������� ��������� ������

	// ����������� �� ���� ������ � ��������� ���� � �������
	for (auto& i : All_gran_)
	{
		for (auto& j : i->yzels)
		{
			j->grans.push_back(i);
		}
	}

	// �������� ��� ����� (��� ����, ����� ���������, ����� ���� �������)
	int kkk = 1;
	for (auto& i : All_gran_)
	{
		i->number = kkk;
		kkk++;
	}

	//����������� �� ���� ����� � ��������� �����, ������� ����� �������
	kkk = 1;
	for (auto& i : this->All_Yzel)
	{
		/*if (kkk % 5000 == 0)
		{
			cout << kkk << "   from  " << this->All_Yzel.size() << endl;
		}*/
		for (int j = 0; j < i->grans.size(); j++)
		{
			for (int k = j + 1; k < i->grans.size(); k++)
			{
				if (areCellsEqual_my(i->grans[j], i->grans[k]) == true)
				{
					auto A1 = number_gran_for_delete.find(i->grans[j]->number);
					auto A2 = number_gran_for_delete.find(i->grans[k]->number);

					if (A1 == number_gran_for_delete.end() && A2 == number_gran_for_delete.end())
					{
						number_gran_for_delete.insert(i->grans[k]->number);
						i->grans[j]->cells.push_back(i->grans[k]->cells[0]);
					}
				}
			}
		}
		kkk++;
	}

	// ������ ���� ������� ��� ������ ����� � ������� �����
	for (auto& i : All_gran_)
	{
		// ���� ����� �� ���� �������, ��������� �
		if(number_gran_for_delete.find(i->number) == number_gran_for_delete.end())
		{
			auto G = new Gran();
			for (auto& j : i->yzels)
			{
				G->yzels.push_back(j);
			}
			for (auto& j : i->cells)
			{
				G->cells.push_back(j);
			}
			this->All_Gran.push_back(G);
		}
		delete i;
	}
	number_gran_for_delete.clear();
	All_gran_.clear();
	All_gran_.resize(0);

	// �������� ��� �����
	kkk = 1;
	for (auto& i : this->All_Gran)
	{
		i->number = kkk;
		kkk++;
	}

	// ������ ������ ����� � �����
	for (auto& i : this->All_Yzel)
	{
		i->grans.clear();
		i->grans.reserve(12);
	}

	// ����������� �� ���� ������ � ��������� ���� � �������
	for (auto& i : this->All_Gran)
	{
		for (auto& j : i->yzels)
		{
			j->grans.push_back(i);
		}
	}

	// ����������� �� ���� ������ � ��������� ����� � ��������
	for (auto& i : this->All_Gran)
	{
		for (auto& j : i->cells)
		{
			j->grans.push_back(i);
		}
	}

	/*if (this->Test_geometr() == false)
	{
		this->~Setka();
		std::exit(EXIT_FAILURE);
	}*/
	cout << "END: New_connect" << endl;

}

void Setka::New_append_surfaces()
{
	// ������� ������� ���� ����� (����� ������� ����� �� ����� ����� �����), 
	// ����� ���������� �� ������ � ����� ������

	for (auto& i : this->All_Yzel)
	{
		i->type = Type_yzel::Us;
	}

	for (auto& i : this->All_Luch)
	{
		if (i->type == "A_Luch" || i->type == "A2_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (j < this->geo->M0 + 1 + this->geo->M1 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_1;
				}
				else if (j < this->geo->M0 + 1 + this->geo->M1 + 1 + this->geo->M2 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_2;
				}
				else if (j < this->geo->M0 + 1 + this->geo->M1 
					+ 1 + this->geo->M2 + 1 + this->geo->M3 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_3;
				}
				else
				{
					i->Yzels[j]->type = Type_yzel::Zone_4;
				}
			}

			i->Yzels_opor[1]->type = Type_yzel::TS;
			i->Yzels_opor[2]->type = Type_yzel::HP;
			i->Yzels_opor[3]->type = Type_yzel::BS;
		}
		else if (i->type == "B_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (j < this->geo->M0 + 1 + this->geo->M1 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_1;
				}
				else if (j < this->geo->M0 + 1 + this->geo->M1 + 1 + this->geo->M2 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_2;
				}
				else if (j < this->geo->M0 + 1 + this->geo->M1
					+ 1 + this->geo->M2 + 1 + this->geo->M3 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_3;
				}
				else
				{
					i->Yzels[j]->type = Type_yzel::Zone_4;
				}
			}

			i->Yzels_opor[1]->type = Type_yzel::TS;
			i->Yzels_opor[2]->type = Type_yzel::HP;
		}
		else if (i->type == "C_Luch" || i->type == "C2_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (j < this->geo->M0 + 1 + this->geo->M1 + 1)
				{
					i->Yzels[j]->type = Type_yzel::Zone_1;
				}
				else
				{
					i->Yzels[j]->type = Type_yzel::Zone_2;
				}
			}

			i->Yzels_opor[1]->type = Type_yzel::TS;
		}
		else if (i->type == "E_Luch")
		{
			short int kl = 1;
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (i->Yzels[j] == i->Yzels_opor[1]) kl++;
				if (i->Yzels[j] == i->Yzels_opor[2]) kl++;

				if(kl == 1) i->Yzels[j]->type = Type_yzel::Zone_2;
				if(kl == 2) i->Yzels[j]->type = Type_yzel::Zone_3;
				if(kl == 3) i->Yzels[j]->type = Type_yzel::Zone_4;
			}

			i->Yzels_opor[2]->type = Type_yzel::Us;

			i->Yzels_opor[1]->type = Type_yzel::HP;
		}
		else if (i->type == "D_Luch")
		{
			short int kl = 1;
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				if (i->Yzels[j] == i->Yzels_opor[1]) kl++;
				if (i->Yzels[j] == i->Yzels_opor[2]) kl++;

				if (kl == 1) i->Yzels[j]->type = Type_yzel::Zone_2;
				if (kl == 2) i->Yzels[j]->type = Type_yzel::Zone_3;
				if (kl == 3) i->Yzels[j]->type = Type_yzel::Zone_4;
			}

			i->Yzels_opor[1]->type = Type_yzel::Us;
			i->Yzels_opor[2]->type = Type_yzel::Us;

			if (i->Yzels_opor[1]->coord[0][0] >= this->geo->L6 - 0.00001)
			{
				i->Yzels_opor[1]->type = Type_yzel::HP;
			}
		}
		else if (i->type == "H_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size(); j++)
			{
				i->Yzels[j]->type = Type_yzel::Zone_2;
			}
		}
		else if (i->type == "G_Luch")
		{
			for (size_t j = 0; j < i->Yzels.size() - 1; j++)
			{
				i->Yzels[j]->type = Type_yzel::Zone_2;
			}
		}
	}

	bool b1 = false;
	// ��� ����������� Gran_TS, Gran_HP, Gran_BS
	macros2(TS);
	macros2(HP);
	macros2(BS);

	for (auto& i : this->Gran_TS)
	{
		i->type2 = Type_Gran_surf::TS;
	}

	for (auto& i : this->Gran_BS)
	{
		i->type2 = Type_Gran_surf::BS;
	}

	for (auto& i : this->Gran_HP)
	{
		i->type2 = Type_Gran_surf::HP;
	}

	this->Renumerate();

	// ������� ���� �������
	if (true)
	{
		for (auto& Cell : this->All_Cell)
		{
			int z1 = 0;
			int z2 = 0;
			int z3 = 0;
			int z4 = 0;

			for (auto& yz : Cell->yzels)
			{
				if (yz->type == Type_yzel::Zone_1)
				{
					z1++;
				}
				if (yz->type == Type_yzel::Zone_2)
				{
					z2++;
				}
				if (yz->type == Type_yzel::Zone_3)
				{
					z3++;
				}
				if (yz->type == Type_yzel::Zone_4)
				{
					z4++;
				}
			}

			if (z1 > z2 && z1 > z3 && z1 > z4)
			{
				Cell->type = Type_cell::Zone_1;
			}
			else if (z2 > z1 && z2 > z3 && z2 > z4)
			{
				Cell->type = Type_cell::Zone_2;
			}
			else if (z3 > z1 && z3 > z2 && z3 > z4)
			{
				Cell->type = Type_cell::Zone_3;
			}
			else
			{
				Cell->type = Type_cell::Zone_4;
			}
		}
	}


	macros4(TS);
	macros4(HP);
	macros4(BS);

}

void Setka::Calculating_measure(unsigned short int st_time)
{
	// ��������� ���� ��� ��� ������������
	/*auto A = this->All_Cell[0];
	macros1(0, 0, 0, 0);
	macros1(1, 1, 0, 0);
	macros1(2, 1, 1, 0);
	macros1(3, 0, 1, 0);
	macros1(4, 0, 0, 1);
	macros1(5, 1, 0, 1);
	macros1(6, 1, 1, 1);
	macros1(7, 0, 1, 1);*/

	// ��������� �������� ���������� �������� ������
	for (auto& i : this->All_Cell)
	{
		i->Culc_center(st_time);
	}

	for (auto& i : this->All_Gran)
	{
		i->Culc_measure(st_time);
	}

	for (auto& i : this->All_Cell)
	{
		i->Culc_volume(st_time);
	}
}

void Setka::auto_set_luch_geo_parameter(int for_new)
{
	// �������������� ��������� �������� ����� � ������ ��������
	// for_new = 0 - ������ ��� ������ ������ ������� � ��������� ������ �� ���������
	// � ���� ������ �������� ���������� ����� �������
	// for_new = 1 - ��������� ��������, ���� ��� �� ������ ������ � ����������� ��� ����� ��� ����

	cout << "Start: Izmenenie geo parameters" << endl;
	// ��� ������� ��������� � �������� �������� ����� (luch.cpp) ����� �������� � ��������� � ��������
	bool izmen = true;
	int k = 0;
	int ii = 0;

	double procent = 50.0;
	if (for_new == 1)
	{
		procent = 1.0;
	}

	while(izmen == true)
	//for(int ii = 0; ii < 50; ii++)
	{
		ii++;

		if (ii > 80)
		{
			procent = 0.1;
		}

		if (for_new == 0) // ������������ ��������� ������������
		{
			if (ii > 80)
			{
				procent = 0.1;
			}
			else if (ii > 50)
			{
				procent = 1.0;
			}
			else if (ii > 20)
			{
				procent = 3.0;
			}
			else if (ii > 13)
			{
				procent = 10.0;
			}
			else if (ii > 3)
			{
				procent = 20.0;
			}
		}



		izmen = false;
		k = 0;

		for (auto& i : this->All_Luch)
		{

			if (i->type == "A_Luch" || i->type == "A2_Luch")
			{
				// da1
				auto a1 = i->get_yzel_near_opor(1, -1);
				auto a2 = i->Yzels_opor[1];

				auto b1 = i->get_yzel_near_opor(1, this->geo->M11);
				auto b2 = i->get_yzel_near_opor(1, this->geo->M11 + 1);

				double d1 = fabs(a1->func_R(0) - a2->func_R(0));
				double d2 = fabs(b1->func_R(0) - b2->func_R(0));

				macros3(da1, 5.0);

				// da2
				b1 = i->Yzels_opor[2];
				b2 = i->get_yzel_near_opor(2, -1);

				d2 = fabs(b1->func_R(0) - b2->func_R(0));

				macros3(da2, 5.0);

				//da3
				if (i->parameters.find("da3") == i->parameters.end() || i->parameters["da3"] > 0.05)
				{
					b1 = i->Yzels_opor[2];
					b2 = i->get_yzel_near_opor(2, 1);

					d2 = fabs(b1->func_R(0) - b2->func_R(0));

					macros3(da3, 5.0);
				}

				// da4
				if (i->parameters.find("da4") == i->parameters.end() || i->parameters["da4"] > 0.05)
				{
					b1 = i->Yzels_opor[3];
					b2 = i->get_yzel_near_opor(3, -1);

					d2 = fabs(b1->func_R(0) - b2->func_R(0));

					macros3(da4, 5.0);
				}

				//da5
				if (i->parameters.find("da5") == i->parameters.end() || i->parameters["da5"] > 0.05)
				{
					b1 = i->Yzels_opor[3];
					b2 = i->get_yzel_near_opor(3, 1);

					d2 = fabs(b1->func_R(0) - b2->func_R(0));

					macros3(da5, 5.0);
				}

			}

			if (i->type == "B_Luch")
			{
				// ba1
				auto a1 = i->get_yzel_near_opor(1, -1);
				auto a2 = i->Yzels_opor[1];

				auto b1 = i->get_yzel_near_opor(1, this->geo->M11);
				auto b2 = i->get_yzel_near_opor(1, this->geo->M11 + 1);

				double d1 = Yzel_distance(a1, a2, 0);
				double d2 = Yzel_distance_x(b1, b2, 0);

				macros3(ba1, 5.0);

				// ba2
				b1 = i->Yzels_opor[2];
				b2 = i->get_yzel_near_opor(2, -1);

				d2 = Yzel_distance(b1, b2, 0);

				macros3(ba2, 3.0);

				// ba3
				if (i->parameters.find("ba3") == i->parameters.end() || i->parameters["ba3"] > 0.05)
				{
					a1 = i->get_yzel_near_opor(1, -1);
					a2 = i->Yzels_opor[1];

					b1 = i->Yzels_opor[2];
					b2 = i->get_yzel_near_opor(2, 1);

					d1 = Yzel_distance(a1, a2, 0);
					d2 = Yzel_distance_x(b1, b2, 0);

					macros3(ba3, 1.0);
				}
				// ba4
				if (i->parameters.find("ba4") == i->parameters.end() || i->parameters["ba4"] > 0.05)
				{
					b1 = i->Yzels_opor[3];
					b2 = i->get_yzel_near_opor(3, -1);

					d2 = Yzel_distance(b1, b2, 0);

					macros3(ba4, 1.0);
				}
			}
		}

		// ��� B ����� �� TS
		for (int i = 0; i < this->B_Luch.size(); i++)
		{
			auto A = this->A_Luch[i].back();
			auto B = this->B_Luch[i][0];
			auto a1 = A->get_yzel_near_opor(3, 1);
			auto a2 = A->Yzels_opor[3];

			auto b1 = B->Yzels_opor[3];
			auto b2 = B->get_yzel_near_opor(3, 1);

			double d1 = Yzel_distance(a1, a2, 0);
			double d2 = Yzel_distance(b1, b2, 0);

			if ((100.0 - d2 * 100.0 / d1) > 1.0)
			{
				k++;
				izmen = true;
				for (auto& kk : this->B_Luch[i])
				{
					if (kk->parameters.find("ba5") != kk->parameters.end())
					{
						kk->parameters["ba5"] *= (1.0 + procent / 100.0);
					}
					else
					{
						kk->parameters["ba5"] = this->geo->ba5 * (1.0 + procent / 100.0);
					}
				}
			}
			else if((100.0 - d2 * 100.0 / d1) < -1.0)
			{
				k++;
				izmen = true;
				for (auto& kk : this->B_Luch[i])
				{
					if (kk->parameters.find("ba5") != kk->parameters.end())
					{
						kk->parameters["ba5"] *= (1.0 - procent / 100.0);
					}
					else
					{
						kk->parameters["ba5"] = this->geo->ba5 * (1.0 - procent / 100.0);
					}
				}
			}
		}

		// ��� E �����
		for (int i = 0; i < this->E_Luch.size(); i++)
		{
			auto A = this->E_Luch[i][0];
			auto B = this->B_Luch[i].back();
			
			// ea1
			auto a1 = A->get_yzel_near_opor(1, 1);
			auto a2 = A->Yzels_opor[1];

			auto b1 = B->Yzels_opor[2];
			auto b2 = B->get_yzel_near_opor(2, 1);

			double d2 = Yzel_distance(a1, a2, 0);
			double d1 = Yzel_distance(b1, b2, 0);

			if (this->E_Luch[i][0]->parameters.find("ea1") == this->E_Luch[i][0]->parameters.end() || this->E_Luch[i][0]->parameters["ea1"] > 0.01)
			{
				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea1") != kk->parameters.end())
						{
							kk->parameters["ea1"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["ea1"] = this->geo->ea1 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea1") != kk->parameters.end())
						{
							kk->parameters["ea1"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["ea1"] = this->geo->ea1 * (1.0 - procent / 100.0);
						}
					}
				}
			}

			// ea2

			if (this->E_Luch[i][0]->parameters.find("ea2") == this->E_Luch[i][0]->parameters.end() || this->E_Luch[i][0]->parameters["ea2"] > 0.01)
			{
				a1 = A->get_yzel_near_opor(2, -1);
				a2 = A->Yzels_opor[2];

				b1 = B->Yzels_opor[3];
				b2 = B->get_yzel_near_opor(3, -1);

				d2 = Yzel_distance(a1, a2, 0);
				d1 = Yzel_distance(b1, b2, 0);

				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea2") != kk->parameters.end())
						{
							kk->parameters["ea2"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["ea2"] = this->geo->ea2 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea2") != kk->parameters.end())
						{
							kk->parameters["ea2"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["ea2"] = this->geo->ea2 * (1.0 - procent / 100.0);
						}
					}
				}
			}

			// ea3

			if (this->E_Luch[i][0]->parameters.find("ea3") == this->E_Luch[i][0]->parameters.end() || this->E_Luch[i][0]->parameters["ea3"] > 0.01)
			{
				a1 = A->get_yzel_near_opor(2, 1);
				a2 = A->Yzels_opor[2];

				b1 = B->Yzels_opor[3];
				b2 = B->get_yzel_near_opor(3, 1);

				d2 = Yzel_distance(a1, a2, 0);
				d1 = Yzel_distance(b1, b2, 0);

				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea3") != kk->parameters.end())
						{
							kk->parameters["ea3"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["ea3"] = this->geo->ea3 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->E_Luch[i])
					{
						if (kk->parameters.find("ea3") != kk->parameters.end())
						{
							kk->parameters["ea3"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["ea3"] = this->geo->ea3 * (1.0 - procent / 100.0);
						}
					}
				}
			}
		}

		// ��� D �����
		for (int i = 0; i < this->D_Luch.size(); i++)
		{
			auto A = this->D_Luch[i][0];
			auto B = this->E_Luch[i].back();

			// md1
			auto a1 = A->get_yzel_near_opor(1, 1);
			auto a2 = A->Yzels_opor[1];

			auto b1 = B->Yzels_opor[1];
			auto b2 = B->get_yzel_near_opor(1, 1);

			double d2 = Yzel_distance(a1, a2, 0);
			double d1 = Yzel_distance(b1, b2, 0);

			if (this->D_Luch[i][0]->parameters.find("md1") == this->D_Luch[i][0]->parameters.end() || this->D_Luch[i][0]->parameters["md1"] > 0.01)
			{
				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md1") != kk->parameters.end())
						{
							kk->parameters["md1"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["md1"] = this->geo->md1 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md1") != kk->parameters.end())
						{
							kk->parameters["md1"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["md1"] = this->geo->md1 * (1.0 - procent / 100.0);
						}
					}
				}
			}

			// md2

			if (this->D_Luch[i][0]->parameters.find("md2") == this->D_Luch[i][0]->parameters.end() || this->D_Luch[i][0]->parameters["md2"] > 0.01)
			{
				a1 = A->get_yzel_near_opor(2, -1);
				a2 = A->Yzels_opor[2];

				b1 = B->Yzels_opor[2];
				b2 = B->get_yzel_near_opor(2, -1);

				d2 = Yzel_distance(a1, a2, 0);
				d1 = Yzel_distance(b1, b2, 0);

				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md2") != kk->parameters.end())
						{
							kk->parameters["md2"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["md2"] = this->geo->md2 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md2") != kk->parameters.end())
						{
							kk->parameters["md2"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["md2"] = this->geo->md2 * (1.0 - procent / 100.0);
						}
					}
				}
			}

			// md3

			if (this->D_Luch[i][0]->parameters.find("md3") == this->D_Luch[i][0]->parameters.end() || this->D_Luch[i][0]->parameters["md3"] > 0.01)
			{
				a1 = A->get_yzel_near_opor(2, 1);
				a2 = A->Yzels_opor[2];

				b1 = B->Yzels_opor[2];
				b2 = B->get_yzel_near_opor(2, 1);

				d2 = Yzel_distance(a1, a2, 0);
				d1 = Yzel_distance(b1, b2, 0);

				if ((100.0 - d2 * 100.0 / d1) > 1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md3") != kk->parameters.end())
						{
							kk->parameters["md3"] *= (1.0 + procent / 100.0);
						}
						else
						{
							kk->parameters["md3"] = this->geo->md3 * (1.0 + procent / 100.0);
						}
					}
				}
				else if ((100.0 - d2 * 100.0 / d1) < -1.0)
				{
					k++;
					izmen = true;
					for (auto& kk : this->D_Luch[i])
					{
						if (kk->parameters.find("md3") != kk->parameters.end())
						{
							kk->parameters["md3"] *= (1.0 - procent / 100.0);
						}
						else
						{
							kk->parameters["md3"] = this->geo->md3 * (1.0 - procent / 100.0);
						}
					}
				}
			}
		}

		// ��� G ����� 
		for (int ik = 0; ik < this->G_Luch.size(); ik++)
		{
				auto A = this->G_Luch[ik][0];
				auto B = this->B_Luch[ik].back();

				auto a1 = A->get_yzel_near_opor(2, -3);
				auto a2 = A->Yzels_opor[2];

				auto b1 = B->Yzels_opor[2];
				auto b2 = B->get_yzel_near_opor(2, -3);

				double d2 = Yzel_distance(a1, a2, 0);
				double d1 = Yzel_distance(b1, b2, 0);

				auto i = this->G_Luch[ik][0];

				if (i->parameters.find("dd2") == i->parameters.end() || i->parameters["dd2"] > 0.01)
				{
					if ((100.0 - d2 * 100.0 / d1) > 0.5)
					{
						k++;
						izmen = true;
						for (auto& kk : this->G_Luch[ik])
						{
							if (kk->parameters.find("dd2") != kk->parameters.end())
							{
								kk->parameters["dd2"] *= (1.0 + procent / 100.0);
							}
							else
							{
								kk->parameters["dd2"] = this->geo->dd2 * (1.0 + procent / 100.0);
							}
						}
					}
					else if ((100.0 - d2 * 100.0 / d1) < -0.5)
					{
						k++;
						izmen = true;
						for (auto& kk : this->G_Luch[ik])
						{
							if (kk->parameters.find("dd2") != kk->parameters.end())
							{
								kk->parameters["dd2"] *= (1.0 - procent / 100.0);
							}
							else
							{
								kk->parameters["dd2"] = this->geo->dd2 * (1.0 - procent / 100.0);
							}
						}
					}
				}
	

		}

		for (auto& i : this->All_Luch)
		{
			i->dvigenie(0);
		}

		//cout << "k = " << k << endl;
	}


	cout << "END: Izmenenie geo parameters" << endl;
}

void Setka::Winslow_method(void)
{
	int N, M;
	double R = 1.0;
	double** x;
	double** y;
	double** z;
	double** x2;
	double** y2;
	double xx, yy;
	double du, dv;
	double dxu, dxv, dyu, dyv, ddxuv, ddyuv, A, B, C, err;

	N = 16;
	M = 16;
	

	// ������ �������
	if (true)
	{
		x = new double* [N];
		for (int i = 0; i < N; i++) x[i] = new double[M];

		y = new double* [N];
		for (int i = 0; i < N; i++) y[i] = new double[M];

		x2 = new double* [N];
		for (int i = 0; i < N; i++) x2[i] = new double[M];

		y2 = new double* [N];
		for (int i = 0; i < N; i++) y2[i] = new double[M];

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				x[i][j] = 0.0;
				x2[i][j] = 0.0;
				y[i][j] = 0.0;
				y2[i][j] = 0.0;
			}
		}
	}

	// ������������� ��������
	if (true)
	{
		for (int i = 0; i < N; i++)
		{
			xx = R * cos(-const_pi/ 2.0 - const_pi/ 4.0 + (i)*const_pi/ 2.0 / (N - 1.0));
			yy = R * sin(-const_pi/ 2.0 - const_pi/ 4.0 + (i)*const_pi/ 2.0 / (N - 1.0));
			x[i][0] = xx;
			x2[i][0] = xx;
			y[i][0] = yy;
			y2[i][0] = yy;
		}

		for (int i = 0; i < N; i++)
		{
			xx = R * cos(const_pi/ 2.0 + const_pi/ 4.0 - (i)*const_pi/ 2.0 / (N - 1.0));
			yy = R * sin(const_pi/ 2.0 + const_pi/ 4.0 - (i)*const_pi/ 2.0 / (N - 1.0));
			x[i][M - 1] = xx;
			x2[i][M - 1] = xx;
			y[i][M - 1] = yy;
			y2[i][M - 1] = yy;
		}

		for (int j = 0; j < M; j++)
		{
			xx = R * cos(-const_pi/ 4.0 + (j)*const_pi/ 2.0 / (M - 1.0));
			yy = R * sin(-const_pi/ 4.0 + (j)*const_pi/ 2.0 / (M - 1.0));
			x[N - 1][j] = xx;
			x2[N - 1][j] = xx;
			y[N - 1][j] = yy;
			y2[N - 1][j] = yy;
		}

		for (int j = 0; j < M; j++)
		{
			xx = R * cos(const_pi+ const_pi/ 4.0 - (j)*const_pi/ 2.0 / (M - 1.0));
			yy = R * sin(const_pi+ const_pi/ 4.0 - (j)*const_pi/ 2.0 / (M - 1.0));
			x[0][j] = xx;
			x2[0][j] = xx;
			y[0][j] = yy;
			y2[0][j] = yy;
		}
	}

	du = 1.0 / N;
	dv = 1.0 / M;

	do
	{
		err = 0.0;
		for (int i = 1; i < N - 1; i++)
		{
			for (int j = 1; j < M - 1; j++)
			{
				dxu = (x[i + 1][j] - x[i - 1][j]) / (2.0 * du);
				dxv = (x[i][j + 1] - x[i][j - 1]) / (2.0 * dv);
				dyu = (y[i + 1][j] - y[i - 1][j]) / (2.0 * du);
				dyv = (y[i][j + 1] - y[i][j - 1]) / (2.0 * dv);
				ddxuv = (x[i + 1][j + 1] + x[i - 1][j - 1] - x[i + 1][j - 1] - x[i - 1][j + 1]) / (4.0 * du * dv);
				ddyuv = (y[i + 1][j + 1] + y[i - 1][j - 1] - y[i + 1][j - 1] - y[i - 1][j + 1]) / (4.0 * du * dv);
				A = dxv * dxv + dyv * dyv;
				B = 2.0 * (dxu * dxv + dyu * dyv) * ddxuv;
				C = dxu * dxu + dyu * dyu;
				/*if (i == 1 && j == 1)
				{
					cout << const_pi<< " " << dxu << " " << dxv << " " << A << " a " << B << " " << C << endl;
					cout << x[i + 1][j] << " " << x[i - 1][j] << " " << du << endl;
				}*/
				if (fabs(2.0 * A / (du * du) + 2.0 * C / (dv * dv)) > 0.000001)
				{
					x2[i][j] = (A * (x[i + 1][j] + x[i - 1][j]) / (du * du) - B + C * (x[i][j + 1] + x[i][j - 1]) / (dv * dv)) / (2.0 * A / (du * du) + 2.0 * C / (dv * dv));
					err = max(err, fabs(x2[i][j] - x[i][j]));
				}
				B = 2.0 * (dxu * dxv + dyu * dyv) * ddyuv;
				if (fabs(2.0 * A / (du * du) + 2.0 * C / (dv * dv)) > 0.000001)
				{
					y2[i][j] = (A * (y[i + 1][j] + y[i - 1][j]) / (du * du) - B + C * (y[i][j + 1] + y[i][j - 1]) / (dv * dv)) / (2.0 * A / (du * du) + 2.0 * C / (dv * dv));
					err = max(err, fabs(y2[i][j] - y[i][j]));
				}
			}
		}
		z = x2;
		x2 = x;
		x = z;

		z = y2;
		y2 = y;
		y = z;

		//cout << err << " " << x[1][1] << " " << x2[1][1] << endl;
		//return;
	} while (err > 0.000001);

	for (int i = 0; i < int(N/4); i++)
	{
		for (int j = 0; j < int(M/4); j++)
		{
			x[i][j] = (x[i][j] - x[N - 1 - i][M - 1 - j] + x[i][M - 1 - j] - x[N - 1 - i][j]) / 4.0;
			x[i][M - 1 - j] = x[i][j];
			x[N - 1 - i][j] = -x[i][j];
			x[N - 1 - i][M - 1 - j] = -x[i][j];

			y[i][j] = (y[i][j] - y[N - 1 - i][M - 1 - j] - y[i][M - 1 - j] + y[N - 1 - i][j]) / 4.0;
			y[i][M - 1 - j] = -y[i][j];
			y[N - 1 - i][j] = y[i][j];
			y[N - 1 - i][M - 1 - j] = -y[i][j];
		}
	}

	ofstream fout;
	fout.open("winslow.txt");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			//if( i%10 == 0 && j%10 == 0)  
			fout << x[i][j] << " " << y[i][j] << endl;
			/*else if( i == N - 1 && j % 10 == 0)  fout << x[i][j] << " " << y[i][j] << endl;
			else if(i % 10 == 0 && j == M - 1)  fout << x[i][j] << " " << y[i][j] << endl;
			else if(i == N - 1 && j == M - 1)  fout << x[i][j] << " " << y[i][j] << endl;*/
		}
	}
	fout.close();

}

void Setka::Tecplot_print_all_yzel_in_3D(string name)
{
	// name - ��� ��� �����
	int k = 1;
	ofstream fout;
	string name_f = "Tecplot_all_yzel_in_3D_" + name + ".txt";
	

	fout.open(name_f);



	fout << "TITLE = HP  VARIABLES = X, Y, Z"  << endl;

	for (auto& ii : this->Yzel_2D)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();


	name_f = "Tecplot_all_kryg_yzel_in_3D_" + name + ".txt";
	k = 1;
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->Krug_Yzel)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();


	name_f = "Tecplot_all_kryg_yzel_2_in_3D_" + name + ".txt";
	k = 1;
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->Krug_Yzel_2)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_yzel_with_condition()
{
	// name - ��� ��� �����
	int k = 1;
	ofstream fout;
	string name_f = "Tecplot_all_yzel_with_condition.txt";


	fout.open(name_f);



	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;
	fout << "ZONE T = HP" << endl;

	for (auto& i : this->All_Yzel)
	{
		if (i->type == Type_yzel::TS)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
	}

	fout.close();
}

void Setka::Tecplot_print_krug_yzel_in_3D(int num)
{
	// num - ����� ���� ������� - �������� ��� ���������
	int k = 1;
	vector <vector<Yzel*>> VVVV;
	ofstream fout;
	string name_f = "Tecplot_krug_yzel_in_3D_" + to_string(num) +  ".txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	VVVV = Krug_Yzel;
	if (num == 2) VVVV = Krug_Yzel_2;

	for (auto& ii : VVVV)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_opor_yzel_in_luchs_3D(string name)
{
	ofstream fout;
	string name_f = "Tecplot_opor_yzel_in_luchs_3D_" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->All_Luch)
	{
		if (ii->type != name) continue;
		for (auto& i : ii->Yzels_opor)
		//for (auto& i : ii->Yzels)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
	}

	fout.close();
}
 
void Setka::Tecplot_print_all_lush_in_3D(string name)
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_all_lush_in_3D_" + name + ".txt";

	auto A = *All_name_luch[name];
	int k = 0;
	for (auto& j : A[0])
	{
		k += (j->Yzels.size() - 1);
	}

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int kkk = 1;

	
	for(auto& i : A)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << 2 * k << ", E = " << k << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& j : i)
		{
			for (int m = 0; m < j->Yzels.size() - 1; m++)
			{
				fout << j->Yzels[m]->coord[0][0] << " " << j->Yzels[m]->coord[0][1] << " " << j->Yzels[m]->coord[0][2] << endl;
				fout << j->Yzels[m + 1]->coord[0][0] << " " << j->Yzels[m + 1]->coord[0][1] << " " << j->Yzels[m + 1]->coord[0][2] << endl;
			}
		}


		for (int m = 0; m < k; m++)
		{
			fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_lush_in_2D()
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_all_lush_in_2D.txt";

	int k = 0;

	for (auto& nn : this->name_luch)
	{
		auto A = *(this->All_name_luch[nn]);
		for (auto& j : A[0])
		{
			k += (j->Yzels.size() - 1);
		}
	}

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;

	int kkk = 1;


	for (int i = 0; i < this->geo->Nphi; i++)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << 2 * k << ", E = " << k << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& nn : this->name_luch)
		{
			auto A = *(this->All_name_luch[nn]);
			for (auto& j : A[i])
			{
				for (int m = 0; m < j->Yzels.size() - 1; m++)
				{
					fout << j->Yzels[m]->coord[0][0] << " " << sqrt(kv(j->Yzels[m]->coord[0][1]) + kv(j->Yzels[m]->coord[0][2])) << endl;
					fout << j->Yzels[m + 1]->coord[0][0] << " " << sqrt(kv(j->Yzels[m + 1]->coord[0][1]) + kv(j->Yzels[m + 1]->coord[0][2])) << endl;
				}
			}
		}


		for (int m = 0; m < k; m++)
		{
			fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_plane_lush(int plane)
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_print_plane_lush_" + to_string(plane) + ".txt";

	int k = 0;

	for (auto& nn : this->name_luch)
	{
		auto A = *(this->All_name_luch[nn]);
		for (auto& j : A[0])
		{
			k += (j->Yzels.size() - 1);
		}
	}

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;


	for (int i = plane; i <= plane; i++)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << 2 * k << ", E = " << k << ", F=FEPOINT, ET=LINESEG" << endl;
	
		for (auto& nn : this->name_luch)
		{
			auto A = *(this->All_name_luch[nn]);
			for (auto& j : A[i])
			{
				for (int m = 0; m < j->Yzels.size() - 1; m++)
				{
					fout << j->Yzels[m]->coord[0][0] << " " << sqrt(kv(j->Yzels[m]->coord[0][1]) + kv(j->Yzels[m]->coord[0][2])) << endl;
					fout << j->Yzels[m + 1]->coord[0][0] << " " << sqrt(kv(j->Yzels[m + 1]->coord[0][1]) + kv(j->Yzels[m + 1]->coord[0][2])) << endl;
				}
			}
		}


		for (int m = 0; m < k; m++)
		{
			fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_plane_surfase(int plane)
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_print_plane_surfase_" + to_string(plane) + ".txt";

	int kTS = 0;
	kTS += this->A_Luch[0].size();
	kTS += this->B_Luch[0].size();
	kTS += this->C_Luch[0].size();

	int kHP = 0;
	kHP += this->A_Luch[0].size();
	kHP += this->B_Luch[0].size();
	kHP += this->E_Luch[0].size();
	kHP += this->D_Luch[0].size();

	int kBS = 0;
	kBS += this->A_Luch[0].size();


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;

	for (int i = plane; i <= plane; i++)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << kTS + kHP + kBS << ", E = " << kTS + kHP + kBS - 3 << ", F=FEPOINT, ET=LINESEG" << endl;

		// TS
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->B_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->C_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		// HP
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[2]->coord[0][0];
			double y = i->Yzels_opor[2]->coord[0][1];
			double z = i->Yzels_opor[2]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->B_Luch[i])
		{
			double x = i->Yzels_opor[2]->coord[0][0];
			double y = i->Yzels_opor[2]->coord[0][1];
			double z = i->Yzels_opor[2]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->E_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->D_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		// BS
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[3]->coord[0][0];
			double y = i->Yzels_opor[3]->coord[0][1];
			double z = i->Yzels_opor[3]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}


		for (int m = 1; m < kTS; m++)
		{
			fout << m << " " << m + 1 << endl;
		}

		for (int m = 1; m < kHP; m++)
		{
			fout << kTS + m << " " << kTS + m + 1 << endl;
		}

		for (int m = 1; m < kBS; m++)
		{
			fout << kTS + kHP + m << " " << kTS + kHP + m + 1 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_All_surfase_in_2D()
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_print_All_surfase_in_2D.txt";

	int kTS = 0;
	kTS += this->A_Luch[0].size();
	kTS += this->B_Luch[0].size();
	kTS += this->C_Luch[0].size();

	int kHP = 0;
	kHP += this->A_Luch[0].size();
	kHP += this->B_Luch[0].size();
	kHP += this->E_Luch[0].size();
	kHP += this->D_Luch[0].size();

	int kBS = 0;
	kBS += this->A_Luch[0].size();


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;
	int kkk = 0;

	for (int i = 0; i < this->geo->Nphi; i++)
	{
		kkk++;
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << kTS + kHP + kBS << ", E = " << kTS + kHP + kBS - 3 << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;

		// TS
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->B_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->C_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		// HP
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[2]->coord[0][0];
			double y = i->Yzels_opor[2]->coord[0][1];
			double z = i->Yzels_opor[2]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->B_Luch[i])
		{
			double x = i->Yzels_opor[2]->coord[0][0];
			double y = i->Yzels_opor[2]->coord[0][1];
			double z = i->Yzels_opor[2]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->E_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		for (auto& i : this->D_Luch[i])
		{
			double x = i->Yzels_opor[1]->coord[0][0];
			double y = i->Yzels_opor[1]->coord[0][1];
			double z = i->Yzels_opor[1]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}

		// BS
		for (auto& i : this->A_Luch[i])
		{
			double x = i->Yzels_opor[3]->coord[0][0];
			double y = i->Yzels_opor[3]->coord[0][1];
			double z = i->Yzels_opor[3]->coord[0][2];
			fout << x << " " << sqrt(kv(y) + kv(z)) << endl;
		}


		for (int m = 1; m < kTS; m++)
		{
			fout << m << " " << m + 1 << endl;
		}

		for (int m = 1; m < kHP; m++)
		{
			fout << kTS + m << " " << kTS + m + 1 << endl;
		}

		for (int m = 1; m < kBS; m++)
		{
			fout << kTS + kHP + m << " " << kTS + kHP + m + 1 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_cell_in_3D()
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_all_cell_in_3D_.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int kkk = 1;


	for (auto& i : this->All_Cell)
	{
		if (i->type != Type_cell::Zone_1) continue;
		if (kkk > 1000) continue;

		fout << "ZONE T=HP, N = " << i->yzels.size() << ", E = " << i->yzels.size() - 1 << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& j : i->yzels)
		{
			fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
		}

		for (int m = 0; m < i->yzels.size() - 1; m++)
		{
			fout << m + 1 << " " << m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_cut_plane_parameters(const Eigen::Vector3d& A,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2)
{
	Eigen::Vector3d B, C;
	Eigen::Vector3d B1, B2;

	Eigen::Vector3d n;
	n = v1.cross(v2);

	bool a1, a2;

	// ������ ��������� ���������



	for (auto& cell : this->All_Cell)
	{
		a1 = false;
		a2 = false;
		// ������� ���� ������ ��������� �� ��������� ������
		for (auto& yz : cell->yzels)
		{
			B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			C = B - A;
			if (n.dot(C) >= 0.0)
			{
				a1 = true;
			}
			else
			{
				a2 = true;
			}
		}

		if (a1 == false || a2 == false) continue; // ������ �� ������������ ����������

		// ����� ���� ����� ��� ������� � �����������
		for (auto& gran : cell->grans)
		{
			int ii;
			for (int i = 0; i < gran->yzels.size(); i++)
			{
				ii = i + 1;
				if (ii >= gran->yzels.size()) ii = 0;
				auto Y1 = gran->yzels[i];
				auto Y2 = gran->yzels[ii];
				a1 = false;
				a2 = false;

				// �������� ��������� �� ��� ���� �� ������ ������� �� ���������
				B1 << Y1->coord[0][0], Y1->coord[0][1], Y1->coord[0][2];
				C = B1 - A;
				if (n.dot(C) >= 0.0)
				{
					a1 = true;
				}
				else
				{
					a2 = true;
				}

				B2 << Y1->coord[0][0], Y1->coord[0][1], Y1->coord[0][2];
				C = B2 - A;
				if (n.dot(C) >= 0.0)
				{
					a1 = true;
				}
				else
				{
					a2 = true;
				}

				if (a1 == false || a2 == false) continue;

			}
		}
	}
}

void Setka::Tecplot_print_cell_plane_parameters()
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_print_cell_plane_parameters.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y";
	for (auto& nn : this->phys_param->param_names)
	{
		fout << ", " << nn;
	}
	fout << ", my_zone, Mach";
	fout << endl;

	int kkk = 1;


	for (auto& i : this->Cell_2D[0])
	{
		fout << i->center[0][0] << " " << norm2(0.0, i->center[0][1], i->center[0][2]);
		for (auto& nn : this->phys_param->param_names)
		{
			if (i->parameters[0].find(nn) != i->parameters[0].end())
			{
				if (nn != "Q")
				{
					fout << " " << i->parameters[0][nn];
				}
				else
				{
					fout << " " << i->parameters[0]["Q"]/ i->parameters[0]["rho"];
				}
			}
			else
			{
				fout << " " << 0.0;
			}
		}

		if (i->type == Type_cell::Zone_1)
		{
			fout << " " << 1.0;
		}
		else if (i->type == Type_cell::Zone_2) 
		{
			fout << " " << 2.0;
		}
		else if (i->type == Type_cell::Zone_3)
		{
			fout << " " << 3.0;
		}
		else if (i->type == Type_cell::Zone_4)
		{
			fout << " " << 4.0;
		}
		else
		{
			fout << " " << 0.0;
		}

		fout << " " << norm2(i->parameters[0]["Vx"], i->parameters[0]["Vy"], i->parameters[0]["Vz"])/
			sqrt(this->phys_param->gamma * i->parameters[0]["p"] / i->parameters[0]["rho"]);

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_gran_in_cell()
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_print_all_gran_in_3D.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (int i = 0; i < 1000; i++) // this->All_Gran.size()
	{
		auto C = this->All_Cell[i];

		fout << "ZONE T=HP, N = " << C->grans.size() * 4 << ", E = " << C->grans.size() << ", F=FEPOINT, ET=quadrilateral, SOLUTIONTIME = " << i + 1 << endl;

		for (auto& A : C->grans)
		{
			for (auto& j : A->yzels)
			{
				fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
			}
		}

		for (int k = 0; k < C->grans.size(); k++)
		{
			fout << 4 * k + 1 << " " << 4 * k + 2 << " " << 4 * k + 3 << " " << 4 * k + 4 << endl;
		}
	}
}

void Setka::Tecplot_print_all_gran_in_surface(string name_surf)
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_print_all_gran_in_surface_" + name_surf + ".txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	vector<Gran*>* C = nullptr;

	if (name_surf == "TS")
	{
		C = &this->Gran_TS;
	}
	else if (name_surf == "HP")
	{
		C = &this->Gran_HP;
	}
	else if (name_surf == "BS")
	{
		C = &this->Gran_BS;
	}
	else
	{
		cout << "Error  9785656643" << endl;
	}


	fout << "ZONE T=HP, N = " << C->size() * 4 << ", E = " << C->size() << ", F=FEPOINT, ET=quadrilateral" <<  endl;

	for (auto& A : *C)
	{
		for (auto& j : A->yzels)
		{
			fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
		}
	}

	for (int k = 0; k < C->size(); k++)
	{
		fout << 4 * k + 1 << " " << 4 * k + 2 << " " << 4 * k + 3 << " " << 4 * k + 4 << endl;
	}
	
}


void Setka::Tecplot_print_gran_with_condition(short int iii)
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_print_gran_with_condition_" + to_string(iii) + "_.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int k = 0; 
	int kk = 0;

	for (auto& i : this->MK_Grans[iii])
	{
		k++;
		kk += i->yzels.size();
	}

	fout << "ZONE T=HP, N = " << kk << ", E = " << k << ", F=FEPOINT, ET=quadrilateral" << endl;

	for (auto& C : this->MK_Grans[iii])
	{
		for (auto& j : C->yzels)
		{
			fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
		}
	}
	
	k = 1;
	for (auto& C : this->MK_Grans[iii])
	{
		for (int i = k; i < k + C->yzels.size(); i++)
		{
			fout << i << " ";
		}
		fout << endl;
		k += C->yzels.size();
	}
}

void Setka::Tecplot_print_1D(Interpol* Int1, const Eigen::Vector3d& Origin, 
	const Eigen::Vector3d& vec, string name, const double& leng)
{
	ofstream fout;
	string name_f = "Tecplot_Tecplot_print_1D" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = R, X, Y, Z";

	for (auto& nam : Int1->param_names)
	{
		fout << ", " << nam;
	}
	fout << ", BB_8Pi";
	fout << endl;

	unsigned int N = 5000;
	Eigen::Vector3d C;
	std::unordered_map<string, double> parameters;
	Cell_handle next_cell;
	Cell_handle prev_cell = Cell_handle();
	bool fine_int;

	for (size_t i = 0; i < N; i++)
	{
		C = Origin + i * vec / N * leng;
		fine_int = Int1->Get_param(C(0), C(1), C(2), parameters, prev_cell, next_cell);
		if (fine_int == false) continue;
		prev_cell = next_cell;

		fout << 1.0 * i / N * leng << " " << C(0) << " " << C(1) << " " << C(2);
		for (auto& nam : Int1->param_names)
		{
			if (nam != "Q")
			{
				fout << " " << parameters[nam];
			}
			else
			{
				fout << " " << parameters["Q"]/ parameters["rho"];
			}
		}
		fout << " " << kvv(parameters["Bx"], parameters["By"], parameters["Bz"]) / (8.0 * const_pi);
		fout << endl;

	}

	fout.close();

}

void Setka::Tecplot_print_2D(Interpol* Int1, const double& a, const double& b, const double& c, const double& d, string name)
{
	// ������� ������� � ���������
	cout << "Start: Tecplot_print_2D " << name << endl;
	std::array<double, 3> normal;
	normal[0] = a;
	normal[1] = b;
	normal[2] = c;

	double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
	if (length > 0)
	{
		normal[0] /= length;
		normal[1] /= length;
		normal[2] /= length;
	}

	// �������� ������������ ����������� ��� ���������� (��������, ��� OX � ���������)
	std::array<double, 3> reference_dir;
	if (std::abs(normal[0]) > 0.9)
	{ // ���� ������� ������ � OX, �������� OY
		reference_dir = { 0.0, 1.0, 0.0 };
	}
	else
	{
		reference_dir = { 1.0, 0.0, 0.0 };
	}

	// ������� ������ � ���������, ���������������� �������
	std::array<double, 3>  tangent_dir = {
		reference_dir[1] * normal[2] - reference_dir[2] * normal[1],
		reference_dir[2] * normal[0] - reference_dir[0] * normal[2],
		reference_dir[0] * normal[1] - reference_dir[1] * normal[0]
	};



	std::vector< std::array<short int, 2> > rebro;
	rebro.resize(12);
	rebro[0][0] = 1; rebro[0][1] = 2;
	rebro[1][0] = 2; rebro[1][1] = 3;
	rebro[2][0] = 3; rebro[2][1] = 4;
	rebro[3][0] = 4; rebro[3][1] = 1;
	rebro[4][0] = 5; rebro[4][1] = 6;
	rebro[5][0] = 6; rebro[5][1] = 7;
	rebro[6][0] = 7; rebro[6][1] = 8;
	rebro[7][0] = 8; rebro[7][1] = 5;
	rebro[8][0] = 1; rebro[8][1] = 5;
	rebro[9][0] = 2; rebro[9][1] = 6;
	rebro[10][0] = 3; rebro[10][1] = 7;
	rebro[11][0] = 4; rebro[11][1] = 8;

	Yzel* A1;
	Yzel* A2;

	std::vector < std::vector< std::array<double, 3> >> all_setka;

	bool aa;

	cout << "Tecplot_print_2D: srez" << endl;

	// ������� ���� �����
	for (const auto& cell : this->All_Cell)
	{
		std::vector< std::array<double, 3> > all_point;
		std::array<double, 3> P1;
		std::array<double, 3> P2;
		std::array<double, 3> outIntersection;

		for (size_t i = 0; i < 12; i++) // ����������� �� �����
		{
			A1 = cell->yzels[rebro[i][0] - 1];
			A2 = cell->yzels[rebro[i][1] - 1];

			P1 = { A1->coord[0][0], A1->coord[0][1], A1->coord[0][2] };
			P2 = { A2->coord[0][0], A2->coord[0][1], A2->coord[0][2] };

			//cout << "A" << endl;
			aa = findIntersection(P1, P2, a, b, c, d, outIntersection);
			//cout << "B" << endl;
			if (aa == true)
			{
				all_point.push_back(outIntersection);
			}
			else
			{
				continue;
			}

			if (fabs(a * outIntersection[0] + b * outIntersection[1] + c * outIntersection[2] + d) > 0.001)
			{
				cout << "Error  798646575467845456" << endl;
				cout << outIntersection[0] << endl;
				cout << outIntersection[1] << endl;
				cout << outIntersection[2] << endl;
				cout << fabs(a * outIntersection[0] + b * outIntersection[1] + c * outIntersection[2] + d) << endl;
				exit(-1);
			}

			//cout << "C" << endl;
		}

		//cout << "________" << endl;

		if (all_point.size() == 0) continue;
		
		if (all_point.size() < 3)
		{
			cout << "Error 9867531090    " << all_point.size() << endl;
			continue;
		}


		// std::vector< std::array<double, 3> > all_point;
		// ������ ��� �������� ��� ��������� �����, ������� ���� �������������� �� �����

		std::array<double, 3> centroid = { 0.0, 0.0, 0.0 };
		for (const auto& p : all_point) {
			centroid[0] += p[0];
			centroid[1] += p[1];
			centroid[2] += p[2];
		}
		centroid[0] /= all_point.size();
		centroid[1] /= all_point.size();
		centroid[2] /= all_point.size();

		/*for (const auto& i : all_point)
		{
			cout << i[0] << " " << i[1] << " " << i[2] << endl;
		}
		cout << "___" << endl;*/


		// ��������� ����� �� ���� ������������ tangent_dir
		std::sort(all_point.begin(), all_point.end(), [centroid, &tangent_dir]
		(std::array<double, 3> aa, std::array<double, 3> bb)
			{
				std::array<double, 3> vec_a = { aa[0] - centroid[0], aa[1] - centroid[1], aa[2] - centroid[2] };
				std::array<double, 3> vec_b = { bb[0] - centroid[0], bb[1] - centroid[1], bb[2] - centroid[2] };;

				// ���� ����� vec_a � tangent_dir
				double dot_a = vec_a[0] * tangent_dir[0] + vec_a[1] * tangent_dir[1] + vec_a[2] * tangent_dir[2];
				double cross_a =
					tangent_dir[1] * vec_a[2] - tangent_dir[2] * vec_a[1] -
					tangent_dir[0] * vec_a[2] + tangent_dir[2] * vec_a[0] +
					tangent_dir[0] * vec_a[1] - tangent_dir[1] * vec_a[0];

				double angle_a = std::atan2(cross_a, dot_a);

				// ���� ����� vec_b � tangent_dir
				double dot_b = vec_b[0] * tangent_dir[0] + vec_b[1] * tangent_dir[1] + vec_b[2] * tangent_dir[2];
				double cross_b =
					tangent_dir[1] * vec_b[2] - tangent_dir[2] * vec_b[1] -
					tangent_dir[0] * vec_b[2] + tangent_dir[2] * vec_b[0] +
					tangent_dir[0] * vec_b[1] - tangent_dir[1] * vec_b[0];

				double angle_b = std::atan2(cross_b, dot_b);

				return angle_a < angle_b;
			});

		/*for (const auto& i : all_point)
		{
			cout << i[0] << " " << i[1] << " " << i[2] << endl;
		}
		cout << endl;

		exit(-1);*/


		all_setka.push_back(all_point);
	}

	unsigned int NN = 0;
	unsigned int NN2 = 0;

	for (const auto& i : all_setka)
	{
		NN += i.size() + 1;
		NN2 += i.size();
	}

	
	cout << "Tecplot_print_2D: print" << endl;

	// ������ ���� �����
	ofstream fout;
	string name_f = "Tecplot_setka_srez_" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP" << endl;
	fout << "VARIABLES = X, Y, Z" << endl;
	fout << "ZONE T=HP, NODES = " << NN2 << ", ELEMENTS = " << NN2 << ", F = FEPOINT, ET = LINESEG" << endl;

	Eigen::Vector3d C;
	for (const auto& i : all_setka)
	{
		for (const auto& j : i)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];
			fout << C(0) << " " << C(1) << " " << C(2) << endl;
		}
	}


	size_t all_k1 = 1;
	for (const auto& i : all_setka)
	{
		size_t k1 = i.size();
		for (size_t ii = 0; ii < k1; ii++)
		{
			size_t k2 = ii + 1;
			if (k2 >= k1) k2 = 0;
			fout << all_k1 + ii << " " << all_k1 + k2 << endl;
		}

		all_k1 = all_k1 + k1;
	}

	fout.close();

	// ������ ���� ����������
	name_f = "Tecplot_Tecplot_print_2D_" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP" << endl;
	fout << "VARIABLES = X, Y, Z";

	for (auto& nam : Int1->param_names)
	{
		fout << ", " << nam;
	}
	fout << ", Mach, BB_8pi";
	fout << endl;

	fout << "ZONE T=HP, ";

	fout << "NODES = " << NN << ", ELEMENTS = " << NN2 << ", F = FEPOINT, ET = TRIANGLE" << endl;

	Eigen::Vector3d C2;
	std::unordered_map<string, double> parameters;
	Cell_handle next_cell;
	Cell_handle prev_cell = Cell_handle();
	bool fine_int;

	for (const auto& i : all_setka)
	{
		std::vector< std::array<double, 3> > kj(i);
		C2 << 0.0, 0.0, 0.0;

		for (const auto& j : i)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];
			C2 += C;
		}
		C2 /= (1.0 * i.size());

		kj.push_back({ C2[0], C2[1], C2[2] });

		bool visible = true;

		for (const auto& j : kj)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];
			fine_int = Int1->Get_param(C(0), C(1), C(2), parameters, prev_cell, next_cell);
			if (fine_int == false)
			{
				C = C * 0.999;
				fine_int = Int1->Get_param(C(0), C(1), C(2), parameters, prev_cell, next_cell);
			}

			if (fine_int == true) next_cell = prev_cell;
			if (fine_int == false)
			{
				visible = false;
				break;
			}
		}

		for (const auto& j : kj)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];
			fine_int = Int1->Get_param(C(0), C(1), C(2), parameters, prev_cell, next_cell);
			if (fine_int == false)
			{
				C = C * 0.999;
				fine_int = Int1->Get_param(C(0), C(1), C(2), parameters, prev_cell, next_cell);
			}

			if (fine_int == true) next_cell = prev_cell;

			if (visible == false)
			{
				fout << 0.0 << " " << 0.0 << " " << 0.0;
			}
			else
			{
				fout << C(0) << " " << C(1) << " " << C(2);
			}


			for (auto& nam : Int1->param_names)
			{
				if (fine_int == true && visible == true)
				{
					if (nam != "Q")
					{
						fout << " " << parameters[nam];
					}
					else
					{
						fout << " " << parameters["Q"] / parameters["rho"];
					}
				}
				else
				{
					fout << " " << 0.0;
				}
			}

			fout << " " << sqrt(parameters["rho"]) * norm2(parameters["Vx"], parameters["Vy"], parameters["Vz"])/
				sqrt(this->phys_param->gamma * parameters["p"]) << " " 
				<< norm2(parameters["Bx"], parameters["By"], parameters["Bz"]) / (8.0 * const_pi);


			fout << endl;
		}

		kj.clear();
	}

	size_t all_k = 1;
	for (const auto& i : all_setka)
	{
		size_t k1 = i.size();
		for (size_t ii = 0; ii < k1; ii++)
		{
			size_t k2 = ii + 1;
			if (k2 >= k1) k2 = 0;
			fout << all_k + k1 << " " << all_k + ii << " " << all_k + k2 << endl;
		}

		all_k = all_k + k1 + 1;
	}

	fout.close();
	cout << "End: Tecplot_print_2D " << name << endl;

}

void Setka::Tecplot_print_plane_interpolation(Eigen::Vector3d A, Eigen::Vector3d v1, 
	Eigen::Vector3d v2, int l1, int r1, int l2, int r2)
{
	cout << "Start: Tecplot_print_plane_interpolation" << endl;
	ofstream fout;
	string name_f = "Tecplot_print_plane_interpolation.txt";
	Eigen::Vector3d P;
	Cell* previos = nullptr;
	Cell* cel = nullptr;
	unordered_map<string, double> par;

	auto names = this->phys_param->param_names;

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y";

	for (auto& nam : names)
	{
		fout << ", " << nam;
	}

	fout << endl;


	for (int i = l1; i <= r1; i++)
	{
		for (int j = l2; j <= r2; j++)
		{
			//cout << "i - " << i << "    j - " << j << endl;
			par.clear();
			P = A + i * v1 + j * v2;

			if (P.norm() < this->geo->R0 - 0.0000001) continue;
			if(P(0) < this->geo->L7 - 0.000001) continue;
			if(norm2(0.0, P(1), P(2)) > this->geo->R5 + 0.000001) continue;
			if(P.norm() > this->geo->R5 + 0.000001 && P(0) > 0.00000001) continue;


			cel = this->Find_cell_point(P(0), P(1), P(2), 0, previos);
			previos = cel;

			if (cel != nullptr)
			{
				//cel->Get_RBF_interpolation(P(0), P(1), P(2), par);
				cel->Get_IDW_interpolation(P(0), P(1), P(2), par, this->phys_param);

				fout << i * v1.norm() << " " << j * v2.norm();
				for (auto& nam : names)
				{
					fout << " " << par[nam];
				}
				fout << endl;
			}
			else
			{
				//cout << "NO = " << P(0) << " " << P(1) << " " << P(2) << endl;
			}
		}
	}

	cout << "End: Tecplot_print_plane_interpolation" << endl;
}