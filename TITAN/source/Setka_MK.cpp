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
	}

	// 2 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (Centr[0] < -0.000001) continue;

		if (gr->type2 == Type_Gran_surf::HP) this->MK_Grans[1].push_back(gr);


		if (gr->type2 == Type_Gran_surf::TS) this->MK_Grans[1].push_back(gr);
		

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2) this->MK_Grans[1].push_back(gr);
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_2)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[1].push_back(gr);
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

		if (gr->type2 == Type_Gran_surf::HP) this->MK_Grans[2].push_back(gr);

		if (gr->type2 == Type_Gran_surf::TS) this->MK_Grans[2].push_back(gr);

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2) this->MK_Grans[2].push_back(gr);
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_3) 
			{
				this->MK_Grans[2].push_back(gr);
			}

			
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_2)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[2].push_back(gr);
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
			if(Centr[0] > 0) this->MK_Grans[3].push_back(gr);
		}

		if (gr->type2 == Type_Gran_surf::BS) this->MK_Grans[3].push_back(gr);

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[3].push_back(gr);
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
			if (Centr[0] < 0.0) this->MK_Grans[4].push_back(gr);
		}

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3) this->MK_Grans[4].push_back(gr);
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_3) this->MK_Grans[4].push_back(gr);


			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[4].push_back(gr);
				}
			}

			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if(gr->cells[0]->center[0][0] < 0) this->MK_Grans[4].push_back(gr);
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
		}

		if (gr->type == Type_Gran::Outer_Hard)
		{
			if(Centr[0] > 0) this->MK_Grans[5].push_back(gr);
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_4 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[5].push_back(gr);
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
			if (Centr[0] < 0) this->MK_Grans[6].push_back(gr);
		}

		if (gr->type != Type_Gran::Us && Centr[0] < 0)
		{
			if(gr->cells[0]->type == Type_cell::Zone_4) this->MK_Grans[6].push_back(gr);
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] < 0) this->MK_Grans[6].push_back(gr);
			}

			if (gr->cells[0]->type == Type_cell::Zone_4 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[6].push_back(gr);
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
	int N_on_gran = 1000;   // Сколько запускаем частиц на грань в среднем
	double mu_expect = 0.0;
	mu_expect = this->MK_Potoks[zone_MK - 1] / 
		(1.0 * N_on_gran * this->MK_Grans[zone_MK - 1].size());

	for (auto& gr : this->MK_Grans[zone_MK - 1])
	{
		// Выбираем конкретный номер датчика случайных чисел
		unsigned int sens_num = 0;
		short int ni = 0; // Номер "входящей" функции распределения
		if (gr->cells[0]->MK_zone == zone_MK)
		{
			ni = 1;
		}
		double full_gran_potok = gr->MK_Potok;

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
			unsigned int N_particle = max(static_cast<int>(func->SpotokV / mu_expect) + 1, 10);
			double mu = func->SpotokV / N_particle; // Вес каждой частицы

			// Запускаем каждую частицу
			for (unsigned int num = 0; num < N_particle; ++num)
			{
				MK_particle P = MK_particle();
				P.mu = mu;

				Eigen::Vector3d poz;

				// Находим положение точки на грани
				gr->Get_Random_pozition(poz, this->Sensors[sens_num]);
				P.coord = poz;

				// Разыгрываем скорость точки
			}
		}
	}
}