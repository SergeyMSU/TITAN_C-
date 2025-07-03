#include "Setka.h"


void Setka::Set_MK_Zone(void)
{
	cout << "Start Set_MK_Zone" << endl;
	this->Renumerate();
	Eigen::Vector3d Centr;

	this->Cell_Center->MK_zone_r = 1;
	this->Cell_Center->MK_zone_phi = 0;

	// Задаём зону для каждой ячейки
	for (auto& cell : this->All_Cell)
	{
		Centr[0] = cell->center[0][0];
		Centr[1] = cell->center[0][1];
		Centr[2] = cell->center[0][2];

		if (cell->type == Type_cell::Zone_1)
		{
			cell->MK_zone_r = 1;
			cell->MK_zone = 1;
			cell->MK_zone_phi = 0;
		}
		else if (cell->type == Type_cell::Zone_2)
		{
			cell->MK_zone_r = 2;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
				cell->MK_zone = 2;
			}
			else
			{
				cell->MK_zone_phi = 2;
				cell->MK_zone = 3;
			}
		}
		else if (cell->type == Type_cell::Zone_3)
		{
			cell->MK_zone_r = 3;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
				cell->MK_zone = 4;
			}
			else
			{
				cell->MK_zone_phi = 2;
				cell->MK_zone = 5;
			}
		}
		else if (cell->type == Type_cell::Zone_4)
		{
			cell->MK_zone_r = 4;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
				cell->MK_zone = 6;
			}
			else
			{
				cell->MK_zone_phi = 2;
				cell->MK_zone = 7;
			}
		}
		else
		{
			cout << "Error 86743207564" << endl;
			exit(-1);
		}
	}

	this->MK_Grans.resize(7);
	this->MK_Potoks.resize(7);
	for (short int i = 0; i < 7; ++i) this->MK_Potoks[i] = 0.0;

	// 1 зона
	for (auto& gr : this->Gran_TS)
	{
		this->MK_Grans[0].push_back(gr);
		gr->MK_type.push_back(1);
	}

	// 2 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (Centr[0] < -0.000001) continue;

		if (gr->type2 == Type_Gran_surf::HP)
		{
			this->MK_Grans[1].push_back(gr);
			gr->MK_type.push_back(2);
		}


		if (gr->type2 == Type_Gran_surf::TS)
		{
			this->MK_Grans[1].push_back(gr);
			gr->MK_type.push_back(2);
		}
		

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2)
			{
				this->MK_Grans[1].push_back(gr);
				gr->MK_type.push_back(2);
			}
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_2)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[1].push_back(gr);
					gr->MK_type.push_back(2);
				}
			}
		}
	}

	// 3 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (Centr[0] > 0.000001) continue;

		if (gr->type2 == Type_Gran_surf::HP)
		{
			this->MK_Grans[2].push_back(gr);
			gr->MK_type.push_back(3);
		}

		if (gr->type2 == Type_Gran_surf::TS)
		{
			this->MK_Grans[2].push_back(gr);
			gr->MK_type.push_back(3);
		}

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2)
			{
				this->MK_Grans[2].push_back(gr);
				gr->MK_type.push_back(3);
			}
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_3) 
			{
				this->MK_Grans[2].push_back(gr);
				gr->MK_type.push_back(3);
			}

			
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_2)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[2].push_back(gr);
					gr->MK_type.push_back(3);
				}
			}
		}
	}

	// 4 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::HP)
		{
			if (Centr[0] > 0)
			{
				this->MK_Grans[3].push_back(gr);
				gr->MK_type.push_back(4);
			}
		}

		if (gr->type2 == Type_Gran_surf::BS)
		{
			this->MK_Grans[3].push_back(gr);
			gr->MK_type.push_back(4);
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[3].push_back(gr);
					gr->MK_type.push_back(4);
				}
			}

		}
	}

	// 5 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::HP)
		{
			if (Centr[0] < 0.0)
			{
				this->MK_Grans[4].push_back(gr);
				gr->MK_type.push_back(5);
			}
		}

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3)
			{
				this->MK_Grans[4].push_back(gr);
				gr->MK_type.push_back(5);
			}
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				this->MK_Grans[4].push_back(gr);
				gr->MK_type.push_back(5);
			}


			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[4].push_back(gr);
					gr->MK_type.push_back(5);
				}
			}

			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] < 0)
				{
					this->MK_Grans[4].push_back(gr);
					gr->MK_type.push_back(5);
				}
			}

		}
	}

	// 6 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::BS)
		{
			this->MK_Grans[5].push_back(gr);
			gr->MK_type.push_back(6);
		}

		if (gr->type == Type_Gran::Outer_Hard)
		{
			if (Centr[0] > 0)
			{
				this->MK_Grans[5].push_back(gr);
				gr->MK_type.push_back(6);
			}
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_4 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[5].push_back(gr);
					gr->MK_type.push_back(6);
				}
			}

		}
	}

	// 7 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type == Type_Gran::Outer_Hard )
		{
			if (Centr[0] < 0)
			{
				this->MK_Grans[6].push_back(gr);
				gr->MK_type.push_back(7);
			}
		}

		if (gr->type != Type_Gran::Us && Centr[0] < 0)
		{
			if (gr->cells[0]->type == Type_cell::Zone_4)
			{
				this->MK_Grans[6].push_back(gr);
				gr->MK_type.push_back(7);
			}
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] < 0)
				{
					this->MK_Grans[6].push_back(gr);
					gr->MK_type.push_back(7);
				}
			}

			if (gr->cells[0]->type == Type_cell::Zone_4 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[6].push_back(gr);
					gr->MK_type.push_back(7);
				}
			}

		}


	}

	for (size_t jj = 0; jj < 7; jj++)
	{
		cout << "MK_grans:  " << jj << "   size = " << this->MK_Grans[jj].size() << endl;
	}

	// Проверяем, что в массивах нет повторов
	if (false)
	{
		for (size_t jj = 0; jj < 7; jj++)
		{
			for (auto& i : this->All_Gran)
			{
				i->work1 = false;
			}

			for (auto& i : this->MK_Grans[jj])
			{
				if (i->work1 == false)
				{
					i->work1 = true;
				}
				else
				{
					cout << "ERROR 6435856408" << endl;
					cout << jj << endl;
				}
			}
		}
			

	}

	cout << "END Set_MK_Zone" << endl;
}

void Setka::MK_prepare(short int zone_MK)
{
	cout << "Start MK_prepare   zone_MK = " << zone_MK << endl;
	// zone_MK должно начинаться с единицы
	if (zone_MK == 0)
	{
		cout << "Error 2341238507" << endl;
		exit(-1);
	}

	this->Renumerate();

	// Блок загрузки датчиков случайных чисел
	if (true)
	{
		ifstream fin2;
		fin2.open("rnd_my.txt");
		if (fin2.is_open() == false)
		{
			cout << "ERROR open  rnd_my.txt " << endl;
			exit(-100);
		}
		double d, a1, b1, c;
		for (int i = 0; i < 1021; i++)
		{
			fin2 >> a1 >> b1 >> c;
			auto s = new Sensor(a1, b1, c);
			this->Sensors.push_back(s);
		}
		fin2.close();
	}

	if (this->MK_Grans.size() < zone_MK || this->MK_Grans[zone_MK - 1].size() == 0)
	{
		cout << "Error 2341963846" << endl;
		exit(-1);
	}

	// Готовим/загружаем AMR сетку для граней
	if (true)
	{
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			if (gr->AMR.size() < this->phys_param->num_H)
			{
				gr->AMR.resize(this->phys_param->num_H);
			}

			for(short int ii = 0; ii <= 1; ii++)\
			{
				for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
				{
					if (gr->AMR[iH - 1][ii] == nullptr)
					{
						string name_f = "func_grans_AMR_" + to_string(ii) + "_H" + 
							to_string(iH) + "_" + to_string(gr->number) + ".bin";
						if (file_exists("data_AMR/" + name_f))
						{
							gr->AMR[iH - 1][ii] = new AMR_f();
							gr->AMR[iH - 1][ii]->AMR_self = gr->AMR[iH - 1][ii];
							gr->AMR[iH - 1][ii]->Read("data_AMR/" + name_f);
						}
						else
						{
							//cout << "ih = " << iH << " " << ii << endl;
							gr->AMR[iH - 1][ii] = new AMR_f(0.0, 20.0, -20.0, 20.0,
								-20.0, 20.0, 3, 6, 6);
							gr->AMR[iH - 1][ii]->AMR_self = gr->AMR[iH - 1][ii];

							if (ii == 0)
							{
								gr->AMR[iH - 1][ii]->Vn[0] = gr->normal[0][0];
								gr->AMR[iH - 1][ii]->Vn[1] = gr->normal[0][1];
								gr->AMR[iH - 1][ii]->Vn[2] = gr->normal[0][2];
							}
							else
							{
								gr->AMR[iH - 1][ii]->Vn[0] = -gr->normal[0][0];
								gr->AMR[iH - 1][ii]->Vn[1] = -gr->normal[0][1];
								gr->AMR[iH - 1][ii]->Vn[2] = -gr->normal[0][2];
							}
							gr->AMR[iH - 1][ii]->Set_bazis();
						}

					}
				}
			}
		}
	}


	// Заполняем граничные условия для AMR сетки (на границе расчётной области)

	if (true)
	{
		for (auto& gr : this->All_boundary_Gran)
		{
			if (gr->type == Type_Gran::Outer_Hard || gr->type == Type_Gran::Outer_Soft)
			{
				if (gr->AMR.size() != 0)
				{
					// Здесь надо задать граничные условия для четвёртого сорта 
					// а также помельчить сетку, если необходимо
					gr->AMR[3][1]->Fill_maxwel_inf(this->phys_param->Velosity_inf);    // <<1>> - внутрь, <<3>> - 4 сорт водорода
					//gr->AMR[3][1]->Print_all_center_Tecplot(gr->AMR[3][1]);
					//gr->AMR[3][1]->Print_1D_Tecplot(gr->AMR[3][1], this->phys_param->Velosity_inf);
					//cout << "DONE  " << endl;
					//exit(-1);
				}
			}
		}
	}

	// Теперь для каждой функции распределения вычисляем поток через неё
	if (true)
	{
		double S = 0.0;
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			gr->Culc_measure(0); // Вычисляем площадь грани (на всякий случай ещё раз)
			// Нужно вычислять поток только у входящей части функции распределения
			short int ni = 0;
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 1;
			}
			gr->MK_Potok = 0.0;

			for (auto& ai : gr->AMR)
			{
					ai[ni]->Culk_SpotokV(gr->area[0]);
					S += ai[ni]->SpotokV;
					gr->MK_Potok += ai[ni]->SpotokV;
			}

		}
		this->MK_Potoks[zone_MK - 1] = S; // Входящий поток через всю границу зоны
	}

	// Счмтаем необходимые геометрические параметры для МК
	if (true)
	{
		for (auto& i : this->All_Cell)
		{
			i->Set_Cell_Geo_for_MK();
		}

		for (auto& i : this->All_Gran)
		{
			i->Set_Gran_Geo_for_MK();
		}
	}

	cout << "END MK_prepare   zone_MK = " << zone_MK << endl;
}

void Setka::MK_delete(short int zone_MK)
{
	cout << "Start MK_delete" << endl;
	// Блок удаления датчиков случайных чисел
	if (true)
	{
		for (auto& i : this->Sensors)
		{
			delete i;
		}
		this->Sensors.clear();
	}

	// Готовим/загружаем AMR сетку для граней
	if (true)
	{
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			for (short int ii = 0; ii <= 1; ii++)\
			{
				for (short int iH = 1; iH <= gr->AMR.size(); iH++)
				{
					string name_f = "func_grans_AMR_" + to_string(ii) + "_H" +
						to_string(iH) + "_" + to_string(gr->number) + ".bin";
					gr->AMR[iH - 1][ii]->Save("data_AMR/" + name_f);
					//cout << "Delete " << ii << " " << iH << " " << 
					//	gr->number << endl;
					gr->AMR[iH - 1][ii]->Delete();
					//cout << "Delete2" << endl;
					delete gr->AMR[iH - 1][ii];
				}
			}
			gr->AMR.clear();
			gr->AMR.resize(0);
		}

		this->MK_Grans.clear();
	}

	cout << "END MK_delete" << endl;
}


void Setka::MK_go(short int zone_MK)
{
	cout << "Start MK_go " << zone_MK << endl;
	int N_on_gran = 10;   // Сколько запускаем частиц на грань в среднем
	double mu_expect = 0.0;
	mu_expect = this->MK_Potoks[zone_MK - 1] / 
		(1.0 * N_on_gran * this->MK_Grans[zone_MK - 1].size());

	unsigned int k1 = 0;
	for (auto& gr : this->MK_Grans[zone_MK - 1])
	{
		k1++;
		if (k1 % 100 == 0)
		{
			cout << "Gran = " << k1 << "    Iz: " << this->MK_Grans[zone_MK - 1].size() << endl;
		}
		// Выбираем конкретный номер датчика случайных чисел
		unsigned int sens_num = 0;
		short int ni = 0; // Номер "входящей" функции распределения
		if (gr->cells[0]->MK_zone == zone_MK)
		{
			ni = 1;
		}
		double full_gran_potok = gr->MK_Potok;

		if (full_gran_potok < 0.000001 * MK_Potoks[zone_MK - 1]/ this->MK_Grans[zone_MK - 1].size())
		{
			continue;
		}

		// Разыгрываем каждый сорт отдельно
		for (short int nh_ = 0; nh_ < this->phys_param->num_H; ++nh_)
		{
			auto& func = gr->AMR[nh_][ni];

			if (func->SpotokV < 0.000001 * full_gran_potok)
			{
				// Функуция распределния нулевая, можно не разыгрывать
				continue;
			}

			// Расчитываем число запускаемых частиц
			unsigned int N_particle = max(static_cast<int>(func->SpotokV / mu_expect) + 1, 
				min(N_on_gran, 10));
			double mu = func->SpotokV / N_particle; // Вес каждой частицы

			// Запускаем каждую частицу
			for (unsigned int num = 0; num < N_particle; ++num)
			{
				MK_particle P = MK_particle();
				if (ni == 0)
				{
					P.cel = gr->cells[1];               // Ячейка в которой находится частица
				}
				else
				{
					P.cel = gr->cells[0];
				}
				
				P.mu = mu;                           // Вес частицы
				P.sort = nh_ + 1;                    // Сорт частицы

				Eigen::Vector3d poz;

				this->Sensors[sens_num]->MakeRandom();
				this->Sensors[sens_num]->MakeRandom();
				this->Sensors[sens_num]->MakeRandom();

				// Находим положение точки на грани
				gr->Get_Random_pozition(poz, this->Sensors[sens_num]);
				P.coord = poz;

				// Находим скорость частицы
				func->Get_random_velosity(func, gr->area[0], poz, this->Sensors[sens_num]);
				P.Vel = poz;

				// Некоторые проверки разыгрынной скорости частицы
				if (this->regim_otladki)
				{
					if (ni == 0)
					{
						if (scalarProductFast(poz(0), poz(1), poz(2),
							gr->normal[0][0], gr->normal[0][1], gr->normal[0][2]) < 0.0)
						{
							cout << "Error  8765656431" << endl;
							exit(-1);
						}
					}
					else
					{
						if (scalarProductFast(poz(0), poz(1), poz(2),
							gr->normal[0][0], gr->normal[0][1], gr->normal[0][2]) > 0.0)
						{
							cout << "Error  7411100090" << endl;
							whach(poz(0));
							whach(poz(1));
							whach(poz(2));
							exit(-1);
						}
					}
				}


				// Проверяем будет ли точка находиться в данной ячейке или нет
				Cell* previos = P.cel;
				double dt = previos->geo_parameters["l_size"] / P.Vel.norm() / 1000.0;
				Cell* ppp = this->Find_cell_point(P.coord[0] + P.Vel[0] * dt,
					P.coord[1] + P.Vel[1] * dt,
					P.coord[2] + P.Vel[2] * dt,
					0, previos);

				if (ppp == nullptr)
				{
					cout << "Error 9756567412" << endl;
					cout << P.coord[0] + P.Vel[0] * dt << " " <<
						P.coord[1] + P.Vel[1] * dt << " " <<
						P.coord[2] + P.Vel[2] * dt << endl;
					cout << P.coord[0] << " " <<
						P.coord[1] << " " <<
						P.coord[2] << endl;
					P.cel->Tecplot_print_cell();
					exit(-1);
				}

				if (ppp != P.cel)
				{
					P.cel = ppp;
					P.coord[0] += P.Vel[0] * dt;
					P.coord[1] += P.Vel[1] * dt;
					P.coord[2] += P.Vel[2] * dt;
				}


				cout << "_______________________________" << endl;
				cout << "Zapusk test  " << P.sort << "    " << gr->number << endl;
				whach(P.coord[0]);
				whach(P.coord[1]);
				whach(P.coord[2]);
				whach(P.Vel[0]);
				whach(P.Vel[1]);
				whach(P.Vel[2]);
				cout << "mu = " << P.mu << endl;
				cout << "_______________________________" << endl;
				this->MK_fly_immit(P, zone_MK); // Запускаем частицу в полёт   // !! Не написана
				cout << "END" << endl;
				exit(-1);
			}
		}
	}

	exit(-1);
}


void Setka::MK_fly_immit(MK_particle& P, short int zone_MK)
{
	/*cout << "______Start_MK_fly_immit___________" << endl;
	whach(P.coord[0]);
	whach(P.coord[1]);
	whach(P.coord[2]);
	cout << "_______________________________" << endl;*/

	cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;


	double time = 0.0;            // время нахождения частицы в ячейке
	Gran* gran = nullptr;         // Через какую грань ячейка выйдер из ячейки


	// Находим время до выхода частицы из ячейки, а также через какую грань будет выход
	bool b1 = false;
	unsigned short int k1 = 0;
	while (b1 == false)
	{
		k1++;
		//cout << "A " << endl;
		b1 = this->Time_to_vilet(P, time, gran);
		if (b1 == true && gran == nullptr)
		{
			cout << "Error 7786341271" << endl;
			exit(-1);
		}
		//cout << "B " << b1 << endl;
		if (b1 == false)
		{
			// Немного двигаем точку
			P.coord += 1e-7 * P.Vel;

			//cout << "C " << endl;
			P.cel =  Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, P.cel);
			// Здесь надо проверить, что во время микро-движения точка не
			// вышла в другую ячейку или за пределы расчётной области
		}

		if (k1 > 100)
		{
			cout << "Error 8614098634" << endl;
			exit(-1);
		}
	}


	/*cout << "D " << time << endl;
	whach(P.coord[0]);
	whach(P.coord[1]);
	whach(P.coord[2]);
	cout << "_______________________________" << endl;*/

	if (gran == nullptr)
	{
		cout << "Error 6438609412" << endl;
		exit(-1);
	}

	// Здесь время до выхода из ячейки определено time
	// Также определено через какую грань это произойдёт  gran
	
	// далее блок основной программы в ячейке





	// Находим следующую ячейку

	P.coord += 1.000001 * time * P.Vel;

	if (gran->Have_zone_number(zone_MK))
	{
		// В этом случае долетели до границы, записываем что надо и выключаем частицу

		return;
	}

	Cell* Cell_next = P.cel->Get_Sosed(gran);
	// точно находим следующую ячейку
	P.cel = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cell_next);
	if (P.cel == nullptr)
	{
		cout << "Error  8545342078" << endl;
		whach(P.coord[0]);
		whach(P.coord[1]);
		whach(P.coord[2]);
		whach(P.Vel[0]);
		whach(P.Vel[1]);
		whach(P.Vel[2]);
		exit(-1);
	}

	return this->MK_fly_immit(P, zone_MK);
}


bool Setka::Time_to_vilet(const MK_particle& P, double& time, Gran*& gran)
{
	Cell* C = P.cel;
	Eigen::Vector3d R, V;

	R = P.coord;
	V = P.Vel;

	Gran* gran_min = nullptr;
	double time_min = 1e10;
	double time1;

	bool b1 = false;

	for (const auto& gr : C->grans)
	{
		if (gr == nullptr)
		{
			cout << "Error 9453286475" << endl;
			exit(-1);
		}

		if (gr->Luch_iz_cross_approx(R, V) == true)
		{
			if (gr->Luch_crossing(R, V, time1) == true)
			{
				//if (time1 > 1e-06) continue;

				if (time_min > time1)
				{
					time_min = time1;
					gran_min = gr;
					b1 = true;
				}
			}
		}
	}

	if (b1 == false)
	{
		// Пересечение ни с одной гранью не произошло!
		return false;
	}
	else
	{
		gran = gran_min;
		time = time_min;

		if (gran == nullptr)
		{
			cout << "Error 1112389046" << endl;
			exit(-1);
		}

		return true;
	}
}