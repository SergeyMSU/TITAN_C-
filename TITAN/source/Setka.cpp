#include "Setka.h"
#include <algorithm>

#define  macros1(n, x, y, z) A->yzels[n]->coord[0][0] = x; \
	A->yzels[n]->coord[0][1] = y; \
	A->yzels[n]->coord[0][2] = z

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



using namespace std;

Setka::Setka()
{
	this->Surf1 = nullptr;
	this->geo = new Geo_param();
	this->phys_param = new Phys_param();
	Luch::geo = this->geo;
	this->geo->Nphi = 60;     // &INIT&
	// ! Это число зависит от сетки триангуляции круга (сколько там вращений по углу она подразумевает)

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

	this->New_initial();                     // Начальное создание узлов и ячеек
	this->New_connect();                     // Начальное создание граней и связывание их с ячейками
	// Также добавляет соседей для каждой ячеки
	this->New_append_surfaces();

	this->Renumerate();
	if (this->Test_geometr() == false)
	{
		this->~Setka();
		std::exit(EXIT_FAILURE);
	}
}

Setka::~Setka()
{
	for (auto& i : this->All_Yzel)
	{
		delete i;
	}
	this->All_Yzel.clear();

	for (auto& i : this->All_Cell)
	{
		delete i;
	}
	this->All_Cell.clear();

	for (auto& i : this->All_Gran)
	{
		delete i;
	}
	this->All_Gran.clear();

	for (auto& i : this->All_Luch)
	{
		delete i;
	}
	this->All_Luch.clear();
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

	cout << "Test 1:" << endl;
	for (auto& i : this->All_Gran)
	{
		if (i->cells.size() != 1 && i->cells.size() != 2)
		{
			cout << "Failure: i->cells.size()  = " << i->cells.size() << endl;
			return false;
		}
	}
	cout << "Success" << endl;


	cout << "Test 2:" << endl;
	for (auto& i : this->All_Cell)
	{
		if (i->grans.size() != 6)
		{
			cout << "Failure: i->grans.size()  = " << i->grans.size() << endl;
			return false;
		}

		if (i->yzels.size() != 8)
		{
			cout << "Failure: i->yzels.size()  = " << i->yzels.size() << endl;
			return false;
		}
	}
	cout << "Success" << endl;

	cout << "Test 3:" << endl;
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
				return false;
			}
		}

		if (i->type == "C_Luch" || i->type == "C2_Luch")
		{
			if (i->Yzels_opor.size() != nnn2)
			{
				cout << "Failure: i->Yzels_opor.size()  = " << i->Yzels_opor.size() << "   " << nnn2 <<
					"   " << i->type << endl;
				return false;
			}
		}
	}
	cout << "Success" << endl << endl;

	return true;
}

void Setka::Set_luch_parametr()
{
	// Добавляет лучам необходимые параметры 
	// (например, углы the и phi для радиальных лучей A, B, C типов) - это ускорит расчёты, так как 
	// они остаются постоянными для луча, и можно постоянно не пересчитывать значения
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
	double R_BS; // положение внешней ударной волны (для движение B E D лучей)
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
	
		// Двигаем BS для B E D лучей
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

	// Считываем файл 2Д сетки
	ifstream ffin("SDK1_2D_Setka.bin", ios::binary | ios::in);
	if (!ffin)
	{
		cout << "Net takogo fajla (fajl 2D setki)" << endl;
		exit(-1);
	}

	ffin.read((char*)&n, sizeof n);  // Версия файла
	cout << "Versiya fajla 2D setki: " << n << endl;

	// Считываем геометрические параметры
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

	// Меняем геометрические параметры, которые нужны
	this->geo->dd3 = 19.0;
	this->geo->dd4 = 15.0;
	this->geo->dd5 = 19.0;
	this->geo->dd6 = 14.0;


	// Считываем узлы
	int sdvig_yz = 0;
	if (true)
	{
		ffin.read((char*)&n, sizeof n); // Число узлов
		for (int i = 0; i < n; i++)
		{
			ffin.read((char*)&x, sizeof x);
			ffin.read((char*)&y, sizeof y);
			auto yzel = new Yzel(x, y, 0.0);  // Создаём узел
			this->All_Yzel.push_back(yzel);
			Yz_.push_back(yzel);
			sdvig_yz++;
		}

		Yzel_2D.push_back(Yz_);  // В этом случае происходит копирование элементов. Списки остаются независимыми, т.е. можно менять Yz_, это не повлияет на Yzel_2D
		Yz_.clear();
	}

	// Считываем лучи
	if (true)
	{
		vector<Luch*> Luch_;
		//for (auto& [name, vec_ptr] : this->All_name_luch)
		for(auto& name : this->name_luch)
		{
			auto vec_ptr = this->All_name_luch[name];
			cout << "Zopolnyaem  " << name << endl;
			ffin.read((char*)&n, sizeof n); // Число лучей
			for (i = 0; i < n; i++)
			{
				auto Luch_0 = new Luch();  // Создаём луч
				All_Luch.push_back(Luch_0);
				Luch_0->type = name;
				Luch_.push_back(Luch_0);
				ffin.read((char*)&n2, sizeof n2); // Число узлов
				for (j = 0; j < n2; j++)
				{
					ffin.read((char*)&n3, sizeof n3); // Номер узла
					if (Yzel_2D[0].size() < static_cast<size_t>(n3))
					{
						cout << "86740576536 PROBLEM   " << Yzel_2D[0].size() << " " << n3;
					}
					Luch_0->Yzels.push_back(Yzel_2D[0][n3 - 1]);
				}

				ffin.read((char*)&n2, sizeof n2); // Число опорных узлов
				//cout << name << "  Opor yzel = " << n2 << endl;
				for (j = 0; j < n2; j++)
				{
					ffin.read((char*)&n3, sizeof n3); // Номер узла
					Luch_0->Yzels_opor.push_back(Yzel_2D[0][n3 - 1]);
				}
			}
			vec_ptr->push_back(Luch_);
			cout << "Dobavili luchey:  " << Luch_.size() << " " << n << endl;
			Luch_.clear();
		}
	}

	// Теперь надо сделать сетку в 3Д, увеличивая число узлов и лучей в ней

	// Увеличиваем число узлов
	if (true)
	{
		for (int iii1 = 1; iii1 < this->geo->Nphi; iii1++)
		{
			Yz_.clear();
			for (auto& i : Yzel_2D[0])
			{
				auto yzel = new Yzel(i->coord[0][0], i->coord[0][1], i->coord[0][2]);  // Создаём узел
				Yz_.push_back(yzel);
				this->All_Yzel.push_back(yzel);
			}
			this->Yzel_2D.push_back(Yz_);

			// Вращаем координаты (поворячиваем плоскость)
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

	// Увеличиваем число лучей
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
					auto Luch_0 = new Luch();  // Создаём луч
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

	// Считываем ячейки
	if (true)
	{
		this->Cell_2D.resize(this->geo->Nphi);

		ffin.read((char*)&n, sizeof n); // Число ячеек
		for (int i = 0; i < this->geo->Nphi; i++)
		{
			this->Cell_2D[i].resize(n);
		}
		this->All_Cell.reserve(this->geo->Nphi* n * 2);

		for (int i = 0; i < n; i++) // Считываем каждую ячейку
		{
			ffin.read((char*)&n2, sizeof n2);
			auto cell_ = new Cell();  // Создаём ячейку
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
			for (int i = 0; i < n; i++) // Считываем каждую ячейку
			{
				n2 = this->Cell_2D[0][i]->yzels.size();
				auto cell_ = new Cell();  // Создаём ячейку
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

	// Теперь надо в каждую ячейку добавить узлы с верхнего слоя узлов
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

	// Считываем файл круга
	ifstream fin("SDK1_krug_setka.bin", ios::binary | ios::in);
	if (!fin)
	{
		cout << "Net takogo fajla (fajl setki v krugu)" << endl;
		exit(-1);
	}

	fin.read((char*)&n, sizeof n);  // Версия файла
	cout << "Versiya fajla setki v krugu: " << n << endl;

	// Считываем узлы (сохраняем их в первый слой головных узлов)
	if (true) 
	{
		this->Krug_Yzel.resize(this->A_Luch[0][0]->Yzels.size());
		this->Krug_Yzel_2.resize(this->C_Luch[0][0]->Yzels.size());

		fin.read((char*)&n, sizeof n); // Число узлов

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
		// Считаем радиус круга

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

		// Двигаем координаты круга 1 на сферу
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

		// Вращаем координаты для совпадения с точками второй сетки
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

		// Двигаем координаты круга 2 на сферу
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

		// Вращаем координаты для совпадения с точками второй сетки
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

	// Перенумеровываем все узлы
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}
	}

	// надо удалить лишние узлы (убрать дубликаты)
	// Следующий цикл ищет дубликаты
	// это можно было бы не делать, если бы точки в кругу в файле лежали в правильном порядке (сначала граничные)
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

	// добавляем узлы в общий массив узлов
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

	// Перенумеровываем все узлы
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}
	}

	// надо сделать лучи из точек на кругах (A2 и C2 типов)  // ---------------------------------------------------------------------------------------------
	if (true)
	{
		this->A2_Luch.reserve(this->Krug_Yzel[0].size() - 60);
		this->C2_Luch.reserve(this->Krug_Yzel_2[0].size() - 60); // Нужно вычесть точки по контуру, так как
			// они уже попали в лучи A

		for (int j = 0; j < this->Krug_Yzel[0].size(); j++) 
		{
			int nf = this->Krug_Yzel[0][j]->number;

			bool hasN = std::any_of(  // Содержит ли вектор элемент, удовлетворяющий условию
				this->A_Luch.begin(),
				this->A_Luch.end(),
				[nf](auto x) { return x[0]->Yzels[0]->number == nf; }
			);

			if (hasN) continue;

			auto Luch_0 = new Luch();  // Создаём луч
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

			bool hasN = std::any_of(  // Содержит ли вектор элемент, удовлетворяющий условию
				this->C_Luch.begin(),
				this->C_Luch.end(),
				[nf](auto x) { return x.back()->Yzels[0]->number == nf; }
			);

			if (hasN) continue;

			auto Luch_0 = new Luch();  // Создаём луч
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


	// Теперь надо двигать узлы в соответствии с функциями движения
	if (true) 
	{
		// SET_PARAMETER   &INIT&
		this->geo->R2 = 20.0;
		this->geo->R3 = 30.0;
		this->geo->R4 = 140.0;
		this->geo->R5 = 400.0;
		this->geo->L6 = -60.0;
		this->geo->L7 = -160.0;
		
		this->Set_luch_parametr();
		// Выставляем опорные точки у А2-лучей в правильном порядке
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

		// Надо бы отсортировать вектор лучей, чтобы двигать их в правильном порядке
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

		// Добавим в E-лучи опорную тотчку из C-луча
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
		//for (auto& nn : this->name_luch) // Что-бы перебирать лучи в правильном порядке
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

	// Считываем ячейки (для точек в кругах)
	if (true)
	{
		this->Cell_layer_head.resize(this->Krug_Yzel.size() - 1);
		this->Cell_layer_tail.resize(this->Krug_Yzel_2.size() - 1);

		int nn, ny;
		fin.read((char*)&n, sizeof n); // Число ячеек

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
			fin.read((char*)&nn, sizeof nn); // Число узлов в ячейке
			for (int j = 0; j < nn; j++)
			{
				fin.read((char*)&ny, sizeof ny); // Номер узла

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

		// Меняем местами узлы, чтобы они шли в правильном порядке
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

	// На настоящий момент ячейки созданы, теперь необходимо их связать (добавить соседей)!
	// Сделать грани, рёбра, связать их
	// Сделать определение ячейки по точке

	fin.close();
	ffin.close();
	cout << "---END New_initial---" << endl;
}

void Setka::New_connect()
{
	// План следующий: сначала создаём для каждой ячейки 6 граней
	// Потом находим повторы, удаляем лишние и связываем грани с ячейками
	// Для ускорения удаления и т.д. сначала всё делаем в локально созданном списке, потом перенесём всё в вектор

	cout << "START: New_connect" << endl;
	vector<Gran*> All_gran_;
	std::unordered_set<int> number_gran_for_delete;

	this->Renumerate();

	// Создаём все грани
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

	// Удаляем дубликаты граней

	// Пробегаемся по всем граням и связываем узлы с гранями
	for (auto& i : All_gran_)
	{
		for (auto& j : i->yzels)
		{
			j->grans.push_back(i);
		}
	}

	// нумеруем все грани (для того, чтобы запомнить, какие надо удалить)
	int kkk = 1;
	for (auto& i : All_gran_)
	{
		i->number = kkk;
		kkk++;
	}

	//пробегаемся по всем узлам и вычисляем грани, которые нужно удалить
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

	// теперь надо удалить все лишние грани и создать новые
	for (auto& i : All_gran_)
	{
		// Если грань не надо удалять, сохраняем её
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

	// нумеруем все грани
	kkk = 1;
	for (auto& i : this->All_Gran)
	{
		i->number = kkk;
		kkk++;
	}

	// чистим старые грани в узлах
	for (auto& i : this->All_Yzel)
	{
		i->grans.clear();
		i->grans.reserve(12);
	}

	// Пробегаемся по всем граням и связываем узлы с гранями
	for (auto& i : this->All_Gran)
	{
		for (auto& j : i->yzels)
		{
			j->grans.push_back(i);
		}
	}

	// Пробегаемся по всем граням и связываем грани с ячейками
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
	// Сначала зададим типы узлам (через опорные точки на лучах проще всего), 
	// потом пробежимся по граням и найдём нужные

	for (auto& i : this->All_Yzel)
	{
		i->type = Type_yzel::Us;
	}

	for (auto& i : this->All_Luch)
	{
		if (i->type == "A_Luch" || i->type == "A2_Luch")
		{
			i->Yzels_opor[1]->type = Type_yzel::TS;
			i->Yzels_opor[2]->type = Type_yzel::HP;
			i->Yzels_opor[3]->type = Type_yzel::BS;
		}
		else if (i->type == "B_Luch")
		{
			i->Yzels_opor[1]->type = Type_yzel::TS;
			i->Yzels_opor[2]->type = Type_yzel::HP;
		}
		else if (i->type == "C_Luch" || i->type == "C2_Luch")
		{
			i->Yzels_opor[1]->type = Type_yzel::TS;
		}
		else if (i->type == "E_Luch")
		{
			i->Yzels_opor[1]->type = Type_yzel::HP;
		}
		else if (i->type == "D_Luch")
		{
			if (i->Yzels_opor[1]->coord[0][0] >= this->geo->L6 - 0.00001)
			{
				i->Yzels_opor[1]->type = Type_yzel::HP;
			}
		}
	}

	bool b1 = false;
	macros2(TS);
	macros2(HP);
	macros2(BS);
}

void Setka::Calculating_measure(unsigned short int st_time)
{
	// Следующий блок был для тестирования
	/*auto A = this->All_Cell[0];
	macros1(0, 0, 0, 0);
	macros1(1, 1, 0, 0);
	macros1(2, 1, 1, 0);
	macros1(3, 0, 1, 0);
	macros1(4, 0, 0, 1);
	macros1(5, 1, 0, 1);
	macros1(6, 1, 1, 1);
	macros1(7, 0, 1, 1);*/

	// Следующий порядкок вычисления является важным
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
	// for_new = 0 - значит это первый запуск функции и параметры двлеки от идеальных
	// в этом случае движение изначально будет большое
	// for_new = 1 - небольшое движение, если это не первый запуск и поверхности уже стоят где надо

	cout << "Start: Izmenenie geo parameters" << endl;
	// Эту функцию совместно с функцией движения сетки (luch.cpp) можно улучшать и дополнять в процессе
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

		if (for_new == 0) // динамическое изменение коэффициента
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

		// для B лучей за TS
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

		// для E лучей
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

		// для D лучей
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
	

	// Создаём массивы
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

	// Инициализация массивов
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
	// name - это имя сетки
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
	// name - это имя сетки
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
	// num - какой круг выводим - головной или хвостовой
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
	// name - это имя лучей
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
	// name - это имя лучей
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
	// name - это имя лучей
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
	// name - это имя лучей
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
	// name - это имя лучей
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
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_all_cell_in_3D_.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int kkk = 1;


	for (auto& i : this->All_Cell)
	{

		if (kkk > 1000)
		{
			break;
		}

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

void Setka::Tecplot_print_all_gran_in_cell()
{
	// name - это имя лучей
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
	// name - это имя лучей
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


void Setka::Tecplot_print_gran_with_condition()
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_print_gran_with_condition.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int k = 0; 
	int kk = 0;

	for (auto& i : this->All_Gran)
	{
		if(i->cells.size() == 1)
		{ 
			k++;
			kk += i->yzels.size();
		}

		if (i->cells.size() <= 0 || i->cells.size() > 2)
		{
			cout << "Error 0909090016" << endl;
		}
	}

	fout << "ZONE T=HP, N = " << kk << ", E = " << k << ", F=FEPOINT, ET=quadrilateral" << endl;

	for (auto& C : this->All_Gran)
	{
		if (C->cells.size() != 1)
		{
			continue;
		}

		for (auto& j : C->yzels)
		{
			fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
		}
	}
	
	k = 1;
	for (auto& C : this->All_Gran)
	{
		if (C->cells.size() != 1)
		{
			continue;
		}

		for (int i = k; i < k + C->yzels.size(); i++)
		{
			fout << i << " ";
		}
		fout << endl;
		k += C->yzels.size();
	}
}