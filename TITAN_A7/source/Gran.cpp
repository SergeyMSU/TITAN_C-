#include "Gran.h"

short int Gran::Get_method()
{
	// 0 - Лакс
	// 1 - HLL
	// 2 - HLLC
	// 3 - HLLD

	//if (this->type2 != Type_Gran_surf::Us) return 3;
	//return 0;

	//if (norm2(0.0, this->center[0][1], this->center[0][2]) >= 330.0) return 0;

	if (this->type2 == Type_Gran_surf::Us)
	{

		if (this->cells[0]->type == Type_cell::Zone_1)
		{
			return 2;
		}

		bool b1 = false;
		bool b2 = false;
		bool b3 = false;
		bool b4 = false;
		for (auto& i : this->yzels)
		{
			if (i->type == Type_yzel::HP)
			{
				return 0;                           // Лакс вдоль гелиопаузы
			}

			if (i->dist_from_HP == 1) b1 = true;
			if (i->dist_from_HP == 2) b2 = true;
			if (i->dist_from_HP == 3) b3 = true;
			if (i->dist_from_HP == 4) b4 = true;

			if (b1 == true || b2 == true || b3 == true || b4 == true)            // Сейчас вообще почти три ряда ячеек считаем Лаксом
			{
				return 0;                            // Лакс во втором ряду от гелиопаузы
			}
		}

		// Временно убираем HLLD
		/*if (this->cells[0]->type == Type_cell::Zone_2 && this->center[0][0] < -70.0)
		{
			return 0;
		}*/

	}

	return 3;
}


void Print_AMR(short int nH, vector<Gran*>& Gran_for_print)
{
	std::vector<std::pair<Point, size_t>> points_1; // точки и их номера для построения триангуляции
	std::vector < std::vector <Int_point*>> ALL_Cells_1(Gran_for_print.size());     // Точки в которых хранятся параметры

	std::vector<AMR_cell*> cells_amr;
	std::array<double, 3> center;

	// Выделяем память под массив указателей
	Delaunay** Delone_array = new Delaunay * [Gran_for_print.size()];
	std::vector <double> koeff(Gran_for_print.size());

	double SumS = 0.0;
	// Инициализируем все указатели nullptr
	for (size_t i = 0; i < Gran_for_print.size(); ++i) 
	{
		SumS += Gran_for_print[i]->area[0];
		Delone_array[i] = nullptr;
	}

	for (size_t i = 0; i < Gran_for_print.size(); ++i)
	{
		koeff[i] = Gran_for_print[i]->area[0] / SumS;
	}

	//Delaunay* Delone_1;


	double Vx, Vy, Vz;
	double x, y, z;
	int NN = 2000;
	double VzL = -10.0;
	double VzR = 10.0;
	double dVz = (VzR - VzL) / NN;
	int Nx = 2000;
	double VxL = -10.0;
	double VxR = 10.0;
	double dVx = (VxR - VxL) / Nx;
	int Ny = 2000;
	double VyL = -10.0;
	double VyR = 10.0;
	double dVy = (VyR - VyL) / Ny;

	boost::multi_array<double, 2> fff(boost::extents[Nx][Ny]);
	for (size_t i = 0; i < Nx; ++i) 
	{
		for (size_t j = 0; j < Ny; ++j) 
		{
			fff[i][j] = 0.0;
		}
	}

	boost::multi_array<double, 1> f1d(boost::extents[Nx]);
	for (size_t i = 0; i < Nx; ++i)
	{
			f1d[i] = 0.0;
	}

	

	//cout << "A1" << endl;
	for (size_t j = 0; j < Gran_for_print.size(); ++j)
	{
		points_1.clear();
		unsigned int i = 0;

		for (const auto& amr : Gran_for_print[j]->AMR[nH - 1])
		{
			cells_amr.clear();
			if (amr == nullptr)
			{
				cout << "Error ertert34tr34t45erfer" << endl;
				return;
			}
			//cout << "B1" << endl;
			amr->Get_all_cells(cells_amr);
			//cout << "B2" << endl;

			for (const auto& cell : cells_amr)
			{
				double u1 = -2.54327;
				double n = 1.0;
				double c = 1.0;

				cell->Get_Center(amr, center);
				amr->Get_real_koordinate(center[0], center[1], center[2], Vx, Vy, Vz);
				auto A = new Int_point(Vx, Vy, Vz);
				//A->parameters["f"] = cell->f;                                     // Так строим саму функцию распределения
				//A->parameters["f"] = cell->f * fabs(center[0]);                     // Так строим поток функции
				//A->parameters["f"] = n / (sqrt_pi * c) * exp(-kv(Vx - u1) / kv(c)); // Так строим поток функции
				A->parameters["f"] = n / pow3(sqrt_pi * c) * exp(-(kv(Vx - u1) + kv(Vy) + kv(Vz))/ kv(c)); // Так строим поток функции


				points_1.push_back({ {Vx, Vy, Vz}, i });
				ALL_Cells_1[j].push_back(A);
				i++;
			}

		}

		//cout << "A2" << endl;

		// Делаем триангуляцию
		Delone_array[j] = new Delaunay(points_1.begin(), points_1.end());
		// Триангулировали функцию распределения, теперь надо её проинтегрировать и вывести в файл
	}

	//cout << "A3" << endl;

	// Инициализируем ограничивающую коробку "пустой"
	CGAL::Bbox_3 bbox;

	// Проходим по всем конечным вершинам триангуляции
	for (size_t j = 0; j < Gran_for_print.size(); ++j)
	{
		for (Delaunay::Finite_vertices_iterator it = Delone_array[j]->finite_vertices_begin();
			it != Delone_array[j]->finite_vertices_end(); ++it) {
			// Получаем точку вершины и добавляем ее bbox в общий bbox
			bbox = bbox + it->point().bbox();
		}
	}
	VzL = bbox.zmin() / 2;
	VzR = bbox.zmax() / 2;
	dVz = (VzR - VzL) / NN;
	VxL = bbox.xmin() / 2;
	VxR = bbox.xmax() / 2;
	dVx = (VxR - VxL) / Nx;
	VyL = bbox.ymin() / 2;
	VyR = bbox.ymax() / 2;
	dVy = (VyR - VyL) / Ny;

	//cout << "A4" << endl;

	// Печатаем 2Д карту
	if (true)
	{
		#pragma omp parallel for schedule(dynamic)
		for (int k1 = 0; k1 < Nx; k1++)
		{
			std::vector < Cell_handle> prev_cell(Gran_for_print.size());
			std::vector < Cell_handle> next_cell(Gran_for_print.size());
			for (short int i = 0; i < Gran_for_print.size(); i++) prev_cell[i] = Cell_handle();

			std::unordered_map<string, double> parameters;
			double Vx = VxL + (VxR - VxL) * (k1 + 0.5) / Nx;
			double SS = 0.0;
			for (int k2 = 0; k2 < Ny; k2++)
			{
				double Vy = VyL + (VyR - VyL) * (k2 + 0.5) / Ny;
				double S = 0.0;

				for (int k = 0; k < 100; k++)
				{
					Vz = VzL + (VzR - VzL) * (k + 0.5) / 100;
					for (size_t j = 0; j < Gran_for_print.size(); ++j)
					{
						if (Get_param_amr(Vx, Vy, Vz, parameters, ALL_Cells_1[j], Delone_array[j], prev_cell[j], next_cell[j]) == true)
						{
							prev_cell[j] = next_cell[j];
							S += parameters["f"] * dVz * koeff[j];
						}
					}
				}
				
				if (S > 0.000001)
				{
					S = 0.0;
						for (int k = 0; k < NN; k++)
						{
							double Vz = VzL + (VzR - VzL) * (k + 0.5) / NN;
							for (size_t j = 0; j < Gran_for_print.size(); ++j)
							{
								if (Get_param_amr(Vx, Vy, Vz, parameters, ALL_Cells_1[j], Delone_array[j], prev_cell[j], next_cell[j]) == true)
								{
									prev_cell[j] = next_cell[j];
									S += parameters["f"] * dVz * koeff[j];
								}
							}
						}
					
				}
				

				SS += S * dVy;
				fff[k1][k2] += S;
			}
			f1d[k1] = SS;
		}
		
	}

	//cout << "A5" << endl;


	//delete Delone_1;
	// Освобождаем память
	for (size_t i = 0; i < Gran_for_print.size(); ++i)
	{
		delete Delone_array[i];
		Delone_array[i] = nullptr;
	}

	cells_amr.clear();
	for (size_t j = 0; j < Gran_for_print.size(); ++j)
	{
		for (auto& cel : ALL_Cells_1[j])
		{
			delete cel;
			cel = nullptr;
		}
		ALL_Cells_1[j].clear();
	}
	points_1.clear();

	// ====== СОЗДАНИЕ ФАЙЛА И ПЕЧАТЬ МАССИВА ======

	// Создаем имя файла на основе номера грани и nH
	cout << "Print AMR: " << Gran_for_print[0]->number << " " << Gran_for_print[0]->center[0][0] << endl;
	cout << "normal: " << Gran_for_print[0]->normal[0][0] << " " << Gran_for_print[0]->normal[0][1] << " " << Gran_for_print[0]->normal[0][2] << endl;
	string filename = "AMR_result_" + to_string(Gran_for_print[0]->number) + "_H" + to_string(nH) + ".txt";
	string filename2 = "1d_AMR_result_" + to_string(Gran_for_print[0]->number) + "_H" + to_string(nH) + ".txt";

	// Открываем файл для записи
	ofstream outfile;
	ofstream outfile2;
	outfile.open(filename);
	outfile2.open(filename2);

	if (!outfile.is_open()) {
		cout << "Error: Cannot open file " << filename << " for writing" << endl;
		return;
	}

	// Заголовок файла в формате Tecplot
	outfile << "TITLE = AMR Distribution" << endl;
	outfile << "VARIABLES = Vx, Vy, f" << endl;
	outfile << "ZONE I = " << Nx << ", J = " << Ny << ", F = POINT" << endl;

	outfile2 << "TITLE = AMR Distribution" << endl;
	outfile2 << "VARIABLES = Vx, f, f_maxwell" << endl;

	// Печатаем массив с координатами
	for (int k2 = 0; k2 < Ny; k2++)  // сначала по Y
	{
		for (int k1 = 0; k1 < Nx; k1++)  // потом по X
		{
			Vx = VxL + (VxR - VxL) * (k1 + 0.5) / Nx;
			Vy = VyL + (VyR - VyL) * (k2 + 0.5) / Ny;
			outfile << Vx << " " << Vy << " " << fff[k1][k2] << endl;
		}
	}

	for (int k1 = 0; k1 < Nx; k1++)  // потом по X
	{
		Vx = VxL + (VxR - VxL) * (k1 + 0.5) / Nx;
		double u1 = -2.54327;
		double u2 = 0.0;
		double u3 = 0.0;
		double v1, v2, v3;
		Gran_for_print[0]->AMR[3][0]->Get_lokal_koordinate(u1, u2, u3, v1, v2, v3);
		double n = 1.0;
		double c = 1.0;
		outfile2 << Vx << " " << fabs(Vx) * f1d[k1] << " " <<
			// n / (sqrt_pi * c) * exp(-kv(Vx - u1) / kv(c))  << endl;
			fabs(Vx) * n / (sqrt_pi * c) * exp((-kv(Vx - v1)) / kv(c)) << endl;
	}

	// Закрываем файл
	outfile.close();
	outfile2.close();

	// Открываем файл для записи информации
	outfile.open("info_AMR_print.txt", std::ios::app);

	outfile << "AMR data saved to file:  " << filename << endl;
	outfile << "Gran centr: " << Gran_for_print[0]->center[0][0] << " " << Gran_for_print[0]->center[0][1] << " " << Gran_for_print[0]->center[0][1] << endl;
}

void Gran::Read_AMR(short int ni, short int nH, Phys_param* ph_param, bool need_refine)
{
	// nH - начинается с 1 !!!!!!!!!!!!!!!!!!!!!!!!!

	if (this->AMR.size() < nH)
	{
		cout << "Error 9yfhgrydttrjik743" << endl;
		exit(-1);
	}

	if (this->AMR[nH - 1][ni] == nullptr)
	{
		this->AMR[nH - 1][ni] = new AMR_f();
		this->AMR[nH - 1][ni]->AMR_self = this->AMR[nH - 1][ni];
	}

	if (this->AMR[nH - 1][ni]->cells.size() > 0)
	{
		cout << "Error 5678567yrewertfewgr  " << endl;
		cout << nH << ni << endl;
		cout << (int)(this->type) << endl;
		cout << (int)(this->type2) << endl;
		exit(-1);
	}

	string name_f = "func_grans_AMR_" + to_string(ni) + "_H" +
		to_string(nH) + "_" + to_string(this->number) + ".bin";
	if (file_exists(ph_param->AMR_folder + "/" + name_f) && this->type == Type_Gran::Us)
	{
		// В этом случае просто считываем AMR - сетку
		this->AMR[nH - 1][ni]->Read(ph_param->AMR_folder + "/" + name_f);

		//if (nH == 4) this->AMR[nH - 1][ni]->procent_signif = 0.1;  // Эксперементально понижаем существенный процент

		// Вроде нормально настроил вторую зону

		if (need_refine == true)
		{
			unsigned int N = this->AMR[nH - 1][ni]->Size();
			unsigned int Nmax = 300000;

			this->AMR[nH - 1][ni]->procent_signif = 0.8;
			this->AMR[nH - 1][ni]->procent_devide = 4.0; // 2.0

			if (nH == 8)
			{
				this->AMR[nH - 1][ni]->procent_signif = 4.0;  // 0.3
				this->AMR[nH - 1][ni]->procent_devide = 7.0; // 2.0
			}

			if (nH == 9)
			{
				this->AMR[nH - 1][ni]->procent_signif = 2.0;  // 0.3
				this->AMR[nH - 1][ni]->procent_devide = 7.0; // 2.0
			}

			if (nH == 4)
			{
				this->AMR[nH - 1][ni]->procent_signif = 0.2;  // 0.3
			}

			if (nH == 3)
			{
				this->AMR[nH - 1][ni]->procent_signif = 0.45;  // 0.3
			}

			if (nH == 2)
			{
				this->AMR[nH - 1][ni]->procent_signif = 1.0; // 0.4
				this->AMR[nH - 1][ni]->procent_devide = 4.0; // 2.0
			}

			if (nH == 1)
			{
				this->AMR[nH - 1][ni]->procent_signif = 0.8;
				this->AMR[nH - 1][ni]->procent_devide = 2.0; // 2.0
			}


			//if (nH == 4) Nmax = 6000;  // Эксперементально
			//if (nH == 3) Nmax = 4000;  // Эксперементально

			if (N < Nmax)
			{
				this->AMR[nH - 1][ni]->Refine(nH);
			}

			/*if (kvv(this->center[0][1], 0.0, this->center[0][2]) < 20.0)
			{
				cout << "x = " << this->center[0][0] << "   nH = " << nH << "     Do = " << 
					N << "    Posle = " << this->AMR[nH - 1][ni]->Size() << endl;
			}*/

			/*else
			{
				cout << "Warning:  AMR = " << N << endl;
				cout << "Info: nH >  " << Nmax << "   =  " << nH << "   gran -> " << this->center[0][0] << " " << this->center[0][1] << " " <<
					this->center[0][2] << endl;
				cout << "gran number = " << this->number << endl;
			}*/
		}

		//this->AMR[nH - 1][ni]->Clean_low();

	}
	else
	{
		if (this->type == Type_Gran::Us)
		{
			if (nH == 1 || nH == 2)
			{
				this->AMR[nH - 1][ni]->AMR_resize(0.0, 50.0, -50.0, 50.0,
					-50.0, 50.0, 3, 6, 6);
			}
			else if (nH == 5 || nH == 6)
			{
				this->AMR[nH - 1][ni]->AMR_resize(0.0, 100.0, -100.0, 100.0,
					-100.0, 100.0, 3, 6, 6);
			}
			else
			{
				this->AMR[nH - 1][ni]->AMR_resize(0.0, 20.0, -20.0, 20.0,
					-20.0, 20.0, 3, 6, 6);
				//if (nH == 4) this->AMR[nH - 1][ni]->procent_signif = 0.1; // Эксперементально понижаем существенный процент
			}
		}
		else
		{
			this->AMR[nH - 1][ni]->AMR_resize(0.0, 20.0, -20.0, 20.0,
				-20.0, 20.0, 1, 1, 1);
		}
	}

	// На всякий случай задаём нормаль
	if (ni == 0)
	{
		this->AMR[nH - 1][ni]->Vn[0] = this->normal[0][0];
		this->AMR[nH - 1][ni]->Vn[1] = this->normal[0][1];
		this->AMR[nH - 1][ni]->Vn[2] = this->normal[0][2];
	}
	else
	{
		this->AMR[nH - 1][ni]->Vn[0] = -this->normal[0][0];
		this->AMR[nH - 1][ni]->Vn[1] = -this->normal[0][1];
		this->AMR[nH - 1][ni]->Vn[2] = -this->normal[0][2];
	}



	this->AMR[nH - 1][ni]->Set_bazis();

	// Заполняем параметры на AMR
	this->AMR[nH - 1][ni]->parameters["n"] = 0.0;
	this->AMR[nH - 1][ni]->parameters["nn"] = 0.0;
	this->AMR[nH - 1][ni]->parameters["Smu"] = 0.0;

}

bool Gran::Have_zone_number(short int z)
{
	for (const auto& i : this->MK_type)
	{
		if (i == z) return true;
	}

	return false;
}

void Gran::Culc_measure(unsigned short int st_time)
{
	double xc, yc, zc;
	int Ny = this->yzels.size();
	xc = 0.0;
	yc = 0.0;
	zc = 0.0;

	// Вычисляем центр грани
	// он понадобиться для проверки нормали
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

	// Старый вариант, где грань = 4 треугольника
	if (false)
	{
		// Вычисляем площадь грани
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
	}
	else  // Грань = 2 треугольника
	{
		// Вычисляем площадь грани
		double S = 0.0;
		S += triangleArea3D(this->yzels[0]->coord[st_time][0], this->yzels[0]->coord[st_time][1],
			this->yzels[0]->coord[st_time][2],
			this->yzels[1]->coord[st_time][0], this->yzels[1]->coord[st_time][1],
			this->yzels[1]->coord[st_time][2],
			this->yzels[2]->coord[st_time][0], this->yzels[2]->coord[st_time][1],
			this->yzels[2]->coord[st_time][2]);

		S += triangleArea3D(this->yzels[0]->coord[st_time][0], this->yzels[0]->coord[st_time][1],
			this->yzels[0]->coord[st_time][2],
			this->yzels[3]->coord[st_time][0], this->yzels[3]->coord[st_time][1],
			this->yzels[3]->coord[st_time][2],
			this->yzels[2]->coord[st_time][0], this->yzels[2]->coord[st_time][1],
			this->yzels[2]->coord[st_time][2]);

		this->area[st_time] = S;
	}


	// Вычисляем нормаль грани
	if (Ny != 4) // для других случаев следующий блок может работать не правильно, 
		// нужно отдельно их рассматривать
	{
		cout << "Error  0906764104" << endl;
		exit(-1);
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

	// Нужно чтобы нормаль смотрела от первой ячейки ко второй
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

void Gran::Get_Random_pozition(Eigen::Vector3d& poz, Sensor* Sens)
{
	vector <double> sqv(2);

	// Вычисляем площадь грани
	double S = 0.0;
	sqv[0] = triangleArea3D(this->yzels[0]->coord[0][0],
		this->yzels[0]->coord[0][1], this->yzels[0]->coord[0][2],
		this->yzels[1]->coord[0][0],
		this->yzels[1]->coord[0][1], this->yzels[1]->coord[0][2],
		this->yzels[2]->coord[0][0],
		this->yzels[2]->coord[0][1], this->yzels[2]->coord[0][2]);

	sqv[1] = triangleArea3D(this->yzels[0]->coord[0][0],
		this->yzels[0]->coord[0][1], this->yzels[0]->coord[0][2],
		this->yzels[3]->coord[0][0],
		this->yzels[3]->coord[0][1], this->yzels[3]->coord[0][2],
		this->yzels[2]->coord[0][0],
		this->yzels[2]->coord[0][1], this->yzels[2]->coord[0][2]);


	S = sqv[0] + sqv[1];



	
	double SS = 0.0;
	double ksi = Sens->MakeRandom();
	Eigen::Vector3d A, B, C;

	for (int i = 0; i < 2; i++)
	{
		if (ksi <= SS + sqv[i] / S || i == 1)
		{
			if (i == 0)
			{
				A << this->yzels[0]->coord[0][0],
					this->yzels[0]->coord[0][1],
					this->yzels[0]->coord[0][2];
				B << this->yzels[1]->coord[0][0],
					this->yzels[1]->coord[0][1],
					this->yzels[1]->coord[0][2];
				C << this->yzels[2]->coord[0][0],
					this->yzels[2]->coord[0][1],
					this->yzels[2]->coord[0][2];
			}
			else
			{
				A << this->yzels[0]->coord[0][0],
					this->yzels[0]->coord[0][1],
					this->yzels[0]->coord[0][2];
				B << this->yzels[3]->coord[0][0],
					this->yzels[3]->coord[0][1],
					this->yzels[3]->coord[0][2];
				C << this->yzels[2]->coord[0][0],
					this->yzels[2]->coord[0][1],
					this->yzels[2]->coord[0][2];
			}

			double r1 = Sens->MakeRandom();
			double r2 = Sens->MakeRandom();
			poz = (1.0 - sqrt(r1)) * A + (sqrt(r1) * (1.0 - r2)) * B +
				(r2 * sqrt(r1)) * C;
			return;
			break;
		}
		else
		{
			SS += sqv[i] / S;
		}
	}

	

}

void Gran::Set_Gran_Geo_for_MK(void)
{
	double x_min = this->yzels[0]->coord[0][0];
	double x_max = this->yzels[0]->coord[0][0];
	double y_min = this->yzels[0]->coord[0][1];
	double y_max = this->yzels[0]->coord[0][1];
	double z_min = this->yzels[0]->coord[0][2];
	double z_max = this->yzels[0]->coord[0][2];

	for (const auto& yz : this->yzels)
	{
		x_min = min(x_min, yz->coord[0][0]);
		y_min = min(y_min, yz->coord[0][1]);
		z_min = min(z_min, yz->coord[0][2]);

		x_max = max(x_max, yz->coord[0][0]);
		y_max = max(y_max, yz->coord[0][1]);
		z_max = max(z_max, yz->coord[0][2]);
	}

	this->geo_parameters["x_min"] = x_min;
	this->geo_parameters["y_min"] = y_min;
	this->geo_parameters["z_min"] = z_min;

	this->geo_parameters["x_max"] = x_max;
	this->geo_parameters["y_max"] = y_max;
	this->geo_parameters["z_max"] = z_max;
}

bool Gran::Luch_iz_cross_approx(const Eigen::Vector3d& R, const Eigen::Vector3d& V)
{
	double tx1, tx2, ty1, ty2, tz1, tz2;
	bool px, py, pz;

	if (this->geo_parameters.find("x_min") == this->geo_parameters.end())
	{
		cout << "Error 9768547854" << endl;
		exit(-1);
	}

	if (this->geo_parameters.find("x_max") == this->geo_parameters.end())
	{
		cout << "Error 9768547852" << endl;
		exit(-1);
	}
		
	px = false;
	if (R[0] > this->geo_parameters["x_min"] &&
		R[0] < this->geo_parameters["x_max"])
	{
		px = true;
	}

	py = false;
	if (R[1] > this->geo_parameters["y_min"] &&
		R[1] < this->geo_parameters["y_max"])
	{
		py = true;
	}

	pz = false;
	if (R[2] > this->geo_parameters["z_min"] &&
		R[2] < this->geo_parameters["z_max"])
	{
		pz = true;
	}



	if (fabs(V[0]) > 1e-6)
	{
		tx1 = (this->geo_parameters["x_min"] - R[0]) / V[0];
		tx2 = (this->geo_parameters["x_max"] - R[0]) / V[0];
		if (tx1 < 0.0 && tx2 < 0.0)
		{
			return false;
		}
	}
	else
	{
		if(px == false)
		{
			return false;
		}
		else
		{
			tx1 = -1e10;
			tx2 = 1e10;
		}
	}

	if (fabs(V[1]) > 1e-6)
	{
		ty1 = (this->geo_parameters["y_min"] - R[1]) / V[1];
		ty2 = (this->geo_parameters["y_max"] - R[1]) / V[1];
		if (ty1 < 0.0 && ty2 < 0.0)
		{
			return false;
		}
	}
	else
	{
		if (py == false)
		{
			return false;
		}
		else
		{
			ty1 = -1e10;
			ty2 = 1e10;
		}
	}

	if (fabs(V[2]) > 1e-6)
	{
		tz1 = (this->geo_parameters["z_min"] - R[2]) / V[2];
		tz2 = (this->geo_parameters["z_max"] - R[2]) / V[2];
		if (tz1 < 0.0 && tz2 < 0.0)
		{
			return false;
		}
	}
	else
	{
		if (pz == false)
		{
			return false;
		}
		else
		{
			tz1 = -1e10;
			tz2 = 1e10;
		}
	}

	double t_enter = max3(min(tx1, tx2), min(ty1, ty2), min(tz1, tz2));
	double t_exit = min3(max(tx1, tx2), max(ty1, ty2), max(tz1, tz2));

	if (t_enter <= t_exit && t_exit >= 0.0)
	{
		return true;
	}

	return false;
}


/**
 * Проверяет пересечение луча с гранью (алгоритм Мёллера — Трумбора).
 *
 * @param orig Начальная точка луча.
 * @param Vel скорость
 * @param t Возвращает время до пересечения.
 * @return true, если луч пересекает треугольник, иначе false.
 */
bool Gran::Luch_crossing(const Eigen::Vector3d& orig, const Eigen::Vector3d& Vel, double& time)
{
	if (this->yzels.size() != 4)
	{
		cout << "Error 7543234305   " << this->yzels.size() << endl;
		exit(-1);
		// Этот алгоритм строго для четырёх-угольных граней, иначе нужен другой 
	}

	auto yz1 = this->yzels[0];
	auto yz2 = this->yzels[1];
	auto yz3 = this->yzels[2];

	Eigen::Vector3d v0, v1, v2;
	v0 << yz1->coord[0][0], yz1->coord[0][1], yz1->coord[0][2];
	v1 << yz2->coord[0][0], yz2->coord[0][1], yz2->coord[0][2];
	v2 << yz3->coord[0][0], yz3->coord[0][1], yz3->coord[0][2];

	Eigen::Vector3d dir;
	dir = Vel;
	double norm = Vel.norm();
	dir /= norm;
	double time1 = 0.0;

	bool b1 = this->rayTriangleIntersect(orig, dir, v0, v1, v2, time1);
	time1 /= norm;


	yz1 = this->yzels[0];
	yz2 = this->yzels[2];
	yz3 = this->yzels[3];
	v0 << yz1->coord[0][0], yz1->coord[0][1], yz1->coord[0][2];
	v1 << yz2->coord[0][0], yz2->coord[0][1], yz2->coord[0][2];
	v2 << yz3->coord[0][0], yz3->coord[0][1], yz3->coord[0][2];
	double time2 = 0.0;

	bool b2 = this->rayTriangleIntersect(orig, dir, v0, v1, v2, time2);
	time2 /= norm;

	if (b1 == false && b2 == false) return false;

	if (b1 == true && b2 == true)
	{
		time = min(time1, time2);
		return true;
	}

	if (b1 == true)
	{
		time = time1;
		return true;
	}

	if (b2 == true)
	{
		time = time2;
		return true;
	}

	return true;
}

/**
 * Проверяет пересечение луча с треугольником (алгоритм Мёллера — Трумбора).
 *
 * @param orig Начальная точка луча (Vector3f).
 * @param dir Направление луча (должно быть нормализовано, Vector3f).
 * @param v0, v1, v2 Вершины треугольника (Vector3f).
 * @param t Возвращает параметр пересечения (расстояние от orig до точки пересечения).
 * @return true, если луч пересекает треугольник, иначе false.
 */
bool Gran::rayTriangleIntersect(
	const Eigen::Vector3d& orig, const Eigen::Vector3d& dir,
	const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
	double& t)
{
	const double EPSILON = 1e-6;

	Eigen::Vector3d edge1 = v1 - v0;
	Eigen::Vector3d edge2 = v2 - v0;
	Eigen::Vector3d pvec = dir.cross(edge2);  // Векторное произведение D ? E2

	double det = edge1.dot(pvec);  // Определитель (для проверки параллельности)

	// Если луч параллелен плоскости треугольника (или почти параллелен)
	if (fabs(det) < EPSILON)
		return false;

	double inv_det = 1.0 / det;

	// Вектор от вершины треугольника до начала луча
	Eigen::Vector3d tvec = orig - v0;

	// Вычисляем барицентрическую координату u
	double u = tvec.dot(pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	// Вектор для вычисления v
	Eigen::Vector3d qvec = tvec.cross(edge1);

	// Вычисляем барицентрическую координату v
	double v = dir.dot(qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	// Вычисляем параметр t (расстояние до пересечения)
	t = edge2.dot(qvec) * inv_det;

	return t > 1e-6;
}

Gran::Gran()
{
	this->yzels.reserve(4);
	this->cells.reserve(2);
	this->cells_TVD.reserve(2);
	this->area[0] = this->area[1] = 0.0;
}

double Gran::culc_velosity(short int now1, const double& time)
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

	ddot = ddot / this->yzels.size()/time;
	return ddot;
}

double Gran::func_R(unsigned short int i_time)
{
	return norm2(this->center[i_time][0], this->center[i_time][1], this->center[i_time][2]);
}