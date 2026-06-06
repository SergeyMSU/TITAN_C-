
#include "Setka.h"
void scatter_1(const Eigen::Vector3d& Vel, double g, double ksi1, double ksi2, Eigen::Vector3d& Vel2);
void perturbVectorKeepingMagnitude_1(double& Vx, double& Vy, double& Vz, const double& noiseScale = 0.01);
double randomNoise_1(double scale);
void perturbVector_1(double& Vx, double& Vy, double& Vz, double noiseScale = 0.01);
bool findSphereIntersectionTime_1(
	const Eigen::Vector3d& X,  // Положение частицы
	const Eigen::Vector3d& V,  // Скорость частицы
	const double& R,        // Радиус сферы
	double& time);
void isotropic_direction(double xi1, double xi2, Eigen::Vector3d& Vel2);


void Setka::MK_go_dust()
{
	auto start = std::chrono::high_resolution_clock::now();
	int N_package = 12000000;   // Сколько запускаем пакетов
	unsigned int k1 = 0;

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
		double a1, b1, c;
		for (int i = 0; i < 1021; i++)
		{
			fin2 >> a1 >> b1 >> c;
			auto s = new Sensor(a1, b1, c);
			this->Sensors.push_back(s);
		}
		fin2.close();
	}

	this->Renumerate();
	// Считаем необходимые геометрические параметры для МК
	// И добавляем переменные в ячейки (заполняем нулями)
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

	// Надо считать температуру пыли из файла
	if (true)
	{
		std::string filename = "dust_paremeter_2-rho10.bin";
		std::ifstream file(filename, std::ios::binary);
		if (!file.is_open())
		{
			std::cerr << "Error gjiuerhgyh7845ygfudhger " << filename << std::endl;
			exit(-1);
		}

		for (auto& A : this->All_Cell)
		{
			double a1, a2;

			file.read(reinterpret_cast<char*>(&a1), sizeof(a1));
			file.read(reinterpret_cast<char*>(&a2), sizeof(a2));
			A->parameters[0]["Tdust"] = A->parameters[1]["Tdust"] = a2;
		}
	}


	for (auto& A : this->All_Cell)
	{
		A->parameters[0]["E_abs"] = 0.0;
	}
	this->DDD->E_esc = 0.0;


	double LL = 4.0 * const_pi * this->DDD->Itot * 1E24;  // Квадрат радиуса звезды
	double Energ = LL / N_package;

	cout << "LL = " << LL << "   Energ = " << Energ << endl;

#pragma omp parallel for schedule(dynamic)
	for (int pc = 1; pc <= N_package; ++pc)  
	{
		#pragma omp critical (first) 
		{
			k1++;
			if (k1 % 50000 == 0)
			{
				cout << "pc = " << k1 << endl;
			}
		}

		// Выбираем конкретный номер датчика случайных чисел
		unsigned int sens_num1 = 2 * omp_get_thread_num();
		unsigned int sens_num2 = 2 * omp_get_thread_num() + 1;

		MK_particle P = MK_particle();
		P.mu = Energ;
		P.lambda = this->DDD->sample_F_kurucz(this->Sensors[sens_num1]->MakeRandom()); // Вычислили длину волны пакета

		// Находим направление движения пакета
		Eigen::Vector3d poz;
		double costhe = 1 - 2.0 * this->Sensors[sens_num1]->MakeRandom();
		double phi = 2.0 * const_pi * this->Sensors[sens_num1]->MakeRandom();
		double sin_theta = std::sqrt(1.0 - costhe * costhe);
		double cos_phi = std::cos(phi);
		double sin_phi = std::sin(phi);
		poz[0] = sin_theta * cos_phi;
		poz[1] = sin_theta * sin_phi;
		poz[2] = costhe;

		//poz[0] = 1.0; poz[1] = 0.0; poz[2] = 0.0; // Направление в апвинд

		P.AddVel(poz);  // Скорость
		P.Addcoord(poz * 25.0); // Положение

		Cell* previos = nullptr;
		Cell* ppp = this->Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, previos);
		P.cel = ppp;               // Ячейка в которой находится частица
		

		P.KSI = -log(1.0 - this->Sensors[sens_num1]->MakeRandom());
		P.I_do = 0.0;

		//cout << k1 << " " << P.lambda << " " << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
		this->MK_fly_dust(P, this->Sensors[sens_num2]); // Запускаем пакет в полёт 
	}


	// Нормировка энергии
	double ff = LL / this->DDD->E_esc;
	cout << "ff = " << ff << endl;
	for (auto& A : this->All_Cell)
	{
		A->parameters[0]["E_abs"] *= ff;
	}

	// Рассчитываем температуру пыли
	// Вычисляем температуру пыли в каждой точке
	for (auto& A : this->All_Cell)
	{
		double Tdust = 0.0;
		if (A->parameters[0]["rhodust"] > 0.00000001)
		{
			double aa = A->parameters[0]["E_abs"] / (4.0 * const_pi) / A->volume[0] / A->parameters[0]["rhodust"] / 3.71063E26; // плотность * объём
			Tdust = DDD->getTemperatureFromL(aa);
		}

		A->parameters[0]["Tdust"] = A->parameters[1]["Tdust"] = Tdust;
	}

	// Надо сохранить температуру и плотность пыли в файл
	if (true)
	{
		std::string filename = "dust_paremeter_1-rho10.bin";
		std::ofstream file(filename, std::ios::binary);
		if (!file.is_open())
		{
			std::cerr << "Error gjiuerhgyh7845ygfudhger " << filename << std::endl;
			exit(-1);
		}

		for (auto& A : this->All_Cell)
		{
			double a1 = A->parameters[0]["rhodust"];
			double a2 = A->parameters[0]["Tdust"];
			file.write(reinterpret_cast<const char*>(&a1), sizeof(a1));
			file.write(reinterpret_cast<const char*>(&a2), sizeof(a2));
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << "MK all time: " << duration.count() / 1000.0 / 60.0 << " minutes" << std::endl;
}


void Setka::MK_fly_dust(MK_particle& P, Sensor* Sens)
{
	Eigen::Vector3d coord_init;
	Eigen::Vector3d Vel_init;

	coord_init << P.coord[0], P.coord[1], P.coord[2];
	Vel_init << P.Vel[0], P.Vel[1], P.Vel[2];

	//cout << "Start " << endl;
	//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;


	unsigned int k_cikl = 0;
	bool vtoroy_shans = false;
	bool vtoroy_shans2 = false;


	std::array<Cell_handle, 6> prev_cell;
	std::array<Cell_handle, 6> next_cell;
	for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();

	Cell* Cell_main = nullptr;   // Ячейка из основной сетки, где находится атом-частица
	// Её нужно знать, так как f_pui и интеграллы от неё хранятся именно в основной сетке
	Cell* Cell_main_prev = nullptr;

	// Главный цикл по ячейкам
	// Выйти из него можно только если частица достигнет конца области
	while (true)
	{
		Eigen::Vector3d coord_do;
		coord_do[0] = P.coord[0];
		coord_do[1] = P.coord[1];
		coord_do[2] = P.coord[2];

		//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;

		k_cikl++;
		if (k_cikl > 10000)
		{
			cout << "Error 8675498765" << endl;
			cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
			cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
			cout << P.sort << " " << P.KSI << " " << P.I_do << endl;
			//exit(-1);
		}

		Cell* Cell_do = P.cel;  // На всякий случай сохраним стартовую ячейку, вдруг надо будет вернуться
		if (Cell_do == nullptr)
		{
			cout << "Error 9865749586" << endl;
			exit(-1);
		}

		double time = 0.0;            // время нахождения частицы в ячейке
		Gran* gran = nullptr;         // Через какую грань ячейка выйдет из ячейки

		//cout << "A1 " << endl;
		//cout << "P.coord = " << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
		// Находим время до выхода частицы из ячейки, а также через какую грань будет выход
		bool b1 = false;
		unsigned short int k1 = 0;

		// Цикл на случай, если точка по какой-то причине не выходит из ячейки
		while (b1 == false)
		{
			k1++;
			if (P.cel == nullptr)
			{
				cout << "Error 12569834678" << endl;
				exit(-1);
			}
			b1 = this->Time_to_vilet(P, time, gran);

			if (b1 == true && gran == nullptr)
			{
				cout << "Error 7786341271" << endl;
				exit(-1);
			}

			if (b1 == false)
			{
				if (P.cel == nullptr)
				{
					cout << "Error 8675463895" << endl;
					exit(-1);
				}
				Eigen::Vector3d Cell_centerr;
				Cell_centerr << P.cel->center[0][0], P.cel->center[0][1],
					P.cel->center[0][2];

				// Подвинем немного точку к центру ячейки

				for (short int i = 0; i < 3; i++)
				{
					/*if (k1 < 2)
					{
						P.coord[i] += 1e-6 * P.Vel[i];
					}
					else if (k1 < 4)
					{
						P.coord[i] += 1e-5 * P.Vel[i];
					}*/
					if (k1 < 10)
					{
						P.coord[i] = P.coord[i] + (Cell_centerr[i] - P.coord[i]) / 800.0;
					}
					else
					{
						P.coord[i] = Cell_centerr[i];
					}
				}
				// Немного двигаем точку
				//P.coord += 1e-6 * P.Vel;

				auto cepp_prev = P.cel;
				P.cel = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, cepp_prev);
				// Здесь надо проверить, что во время микро-движения точка не
				// вышла в другую ячейку или за пределы расчётной области
				if (P.cel == nullptr)
				{
					//cout << "TUT  1875408695" << endl;
					P.cel = Cell_do;
				}
			}

			if (k1 > 11)
			{
				if (vtoroy_shans2 == false)
				{
					vtoroy_shans2 = true;
					P.cel = Cell_do;
					P.coord[0] = coord_do[0];
					P.coord[1] = coord_do[1];
					P.coord[2] = coord_do[2];
					perturbVectorKeepingMagnitude_1(P.Vel[0], P.Vel[1], P.Vel[2], 0.1 * norm2(P.Vel[0], P.Vel[1], P.Vel[2]));



					continue;
				}

				cout << "Poteryal D" << endl;
				//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
				//cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
				//cout << coord_do[0] << " " << coord_do[1] << " " << coord_do[2] << endl;
				//cout << P.cel->center[0][0] << " " << P.cel->center[0][1] << " " << P.cel->center[0][2] << endl;

				return;
				//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
				//whach(P.cel->number);

				//P.cel->Tecplot_print_cell();
				//exit(-1);
			}
		}

		//cout << "B " << endl;
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

		vtoroy_shans2 = false;

		// Здесь время до выхода из ячейки определено time
		// Также определено через какую грань это произойдёт  gran

		// далее блок основной программы в ячейке
		// ****************************************************************************

		// Получаем параметры в ячейке ----------------------------

		double rhodust = P.cel->parameters[0]["rhodust"];    // плотность пыли
		double Tdust = P.cel->parameters[0]["Tdust"];        // температура пыли

		double I = P.I_do;
		double l = sqrt(kvv(time * P.Vel[0], time * P.Vel[1], time * P.Vel[2]));    // Расстояние, которое атом потенциально пролетает внутри ячейки
		double Vel_norm = sqrt(kvv(P.Vel[0], P.Vel[1], P.Vel[2]));                  // Модуль скорости атома

		double Kabs = this->DDD->interpolate_K_abs(P.lambda);
		double Ksca = this->DDD->interpolate_K_sca(P.lambda);
		double nu_ex = (Kabs + Ksca) * rhodust;

		double sig = 0.0;

		if (nu_ex >= 0.000000001)
		{
			// Иначе если частота процессов нулевая, то в этой ячейке не произошло никакое событие
			sig = Vel_norm / nu_ex;
			I += l / sig * 3.35078E-7;  // Множитель из-за размера * характерную плотность
		}

		if (vtoroy_shans == false)
		{
			if (I < P.KSI)
			{
				P.I_do = I;  // В этом случае перезарядки в ячейке не произошло

				P.cel->mut.lock();
				P.cel->parameters[0]["E_abs"] += P.mu * (Kabs * rhodust) * l * 3.35078E-7;
				P.cel->mut.unlock();
			}
			else
			{
				double ksi = (P.KSI - P.I_do) * sig / 3.35078E-7;
				double t_ex = ksi / Vel_norm;
				if (t_ex > time * 1.0001)
				{
					cout << "Error e84ut7yhe8fh[oejfr  " << t_ex << " " << time << endl;
				}
				P.I_do = 0.0;
				for (short int i = 0; i < 3; i++) P.coord[i] += t_ex * P.Vel[i];
				Cell* Cnow = P.cel;
				Cell* CC = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cnow);

				if (P.cel != CC)
				{
					// Если перезарядка произошла за пределами текущей ячейки
					for (short int i = 0; i < 3; i++)
					{
						P.coord[i] -= t_ex / 1000.0 * P.Vel[i];
					}

					Cnow = CC;
					CC = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cnow);
					if (P.cel != CC)
					{
						Eigen::Vector3d Cell_center;
						Cell_center << P.cel->center[0][0], P.cel->center[0][1],
							P.cel->center[0][2];
						unsigned short int kklk = 0;
					dchj12:
						kklk++;
						if (kklk > 20)
						{
							cout << "Poteryal C" << endl;
							return;
						}

						// Подвинем немного точку к центру ячейки
						for (short int i = 0; i < 3; i++)
						{
							if (kklk < 18)
							{
								P.coord[i] += (Cell_center[i] - P.coord[i]) / 100.0;
							}
							else
							{
								P.coord[i] = Cell_center[i];
							}
						}

						Cnow = CC;
						CC = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cnow);

						if (P.cel != CC)
						{
							goto dchj12;
							cout << "Error 8674539765" << endl;
							whach(CC->number);
							whach(P.cel->number);
							whach(P.coord[0]);
							whach(P.coord[1]);
							whach(P.coord[2]);
							whach(P.Vel[0]);
							whach(P.Vel[1]);
							whach(P.Vel[2]);
							whach(t_ex);
							whach(time);
							exit(-1);
						}
					}
				}

				//double uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * const_pi * sqrtpi_);
				//double uz_E = Velosity_3(u, cp);

				P.cel->mut.lock();
				P.cel->parameters[0]["E_abs"] += P.mu * (this->DDD->interpolate_K_abs(P.lambda) * rhodust) * (t_ex * Vel_norm) * 3.35078E-7;
				P.cel->mut.unlock();

				// -------------------------------------------------------------------
				// теперь нужно определить процесс, который произошёл
				double ksi_ = Sens->MakeRandom();

				double omega = Ksca / (Ksca + Kabs);

				if (ksi_ < omega)
				{
					//return;  // Пока просто вырубим эти пакеты
					// В этом случае произошло рассеяние пакета
					double g = this->DDD->interpolate_g(P.lambda);
					//if (g < 0.9) return; // Не вырубаем пакеты, которые почти не рассеялись
					double ksi1 = Sens->MakeRandom();
					double ksi2 = Sens->MakeRandom();
					Eigen::Vector3d Vel, Vel2;
					Vel[0] = P.Vel[0];
					Vel[1] = P.Vel[1];
					Vel[2] = P.Vel[2];
					scatter_1(Vel, g, ksi1, ksi2, Vel2);  // Разыгрываем скорость пакета при рассеянии
					P.Vel[0] = Vel2[0];
					P.Vel[1] = Vel2[1];
					P.Vel[2] = Vel2[2];
				}
				else
				{
					// В этом случае произошло поглощение пакета
					//return;  // Пока просто вырубим эти пакеты
					// Здесь надо поменять частоту покеты и выбрать изотропное направление
					double ksi1 = Sens->MakeRandom();
					double ksi2 = Sens->MakeRandom();
					double ksi3 = Sens->MakeRandom();
					P.lambda = this->DDD->sample_frequency_sca(ksi1, Tdust);

					Eigen::Vector3d Vel2;
					isotropic_direction(ksi2, ksi3, Vel2);
					P.Vel[0] = Vel2[0];
					P.Vel[1] = Vel2[1];
					P.Vel[2] = Vel2[2];

					//this->DDD->mut.lock();
					//this->DDD->E_esc += P.mu;
					//this->DDD->mut.unlock();
					//return;  // Пока просто вырубим эти пакеты
				}

				P.KSI = -log(1.0 - Sens->MakeRandom());
				vtoroy_shans = false;
				continue;
			}
		}

		// ****************************************************************************
		// Находим следующую ячейку
		for (short int i = 0; i < 3; i++)
		{
			P.coord[i] += 1.000001 * time * P.Vel[i];
		}

		if (norm2(P.coord[0], P.coord[1], P.coord[2]) < 1.01 * this->geo->R0)
		{
			// Частица попала во внутреннюю сферу, надо, чтобы они пролетели мимо неё
			Eigen::Vector3d X(P.coord[0], P.coord[1], P.coord[2]);
			Eigen::Vector3d V(P.Vel[0], P.Vel[1], P.Vel[2]);
			double time_;
			if (findSphereIntersectionTime_1(X, V, 1.01 * this->geo->R0,
				time_) == true)
			{
				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += time_ * P.Vel[i];
				}
			}
			else
			{
				double norm_ = norm2(P.Vel[0], P.Vel[1], P.Vel[2]);
				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += (2.02 * this->geo->R0) * P.Vel[i] / norm_;
				}
			}
		}

		if (gran->type != Type_Gran::Us)
		//if (gran->type != Type_Gran::Us || norm2(P.coord[0], P.coord[1], P.coord[2]) > 55.635)  // Вырубаем за пределами внешнего ударного слоя
		{
		a1:
			// В этом случае долетели до границы, записываем что надо и выключаем частицу

			gran->mut.lock(); // Мьютекс для записи в гранб
			gran->N_particle++;
			gran->mut.unlock();

			this->DDD->mut.lock();
			this->DDD->E_esc += P.mu;
			this->DDD->mut.unlock();
			return;
		}

		// Ручная проверка вылета за пределы области
		if (P.coord[0] < this->geo->L7 || (P.coord[0] < 0.0 && norm2(0.0, P.coord[1], P.coord[2]) > this->geo->R5) ||
			(P.coord[0] >= 0.0 && norm2(P.coord[0], P.coord[1], P.coord[2]) > this->geo->R5))
		{
			this->DDD->mut.lock();
			this->DDD->E_esc += P.mu;
			this->DDD->mut.unlock();
			return;
		}

		Cell* Cell_next = P.cel->Get_Sosed(gran);
		short unsigned int kkk2 = 0;
	vv1:
		kkk2++;
		// точно находим следующую ячейку
		Cell* Cell_next_ = Cell_next;
		P.cel = Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, Cell_next_);

		// Кажется в случае проблем надо просто подтянуть ячейку к центру грани

		if (P.cel == nullptr)
		{
			for (auto& gr : Cell_do->grans)
			{
				gran = gr;
				goto a1;
			}
		}

		if (P.cel == nullptr)
		{
			// В этом случае точка часто попадает в угол ячеки или на грань

			if (kkk2 < 3)
			{
				Eigen::Vector3d normal;
				normal << gran->normal[0][0], gran->normal[0][1], gran->normal[0][2];
				if (gran->cells[0] != Cell_do)
				{
					normal *= -1.0;
				}
				double l_ = norm2(gran->yzels[0]->coord[0][0] - gran->yzels[1]->coord[0][0],
					gran->yzels[0]->coord[0][1] - gran->yzels[1]->coord[0][1],
					gran->yzels[0]->coord[0][2] - gran->yzels[1]->coord[0][2]);
				l_ = min(l_, norm2(gran->yzels[0]->coord[0][0] - gran->yzels[2]->coord[0][0],
					gran->yzels[0]->coord[0][1] - gran->yzels[2]->coord[0][1],
					gran->yzels[0]->coord[0][2] - gran->yzels[2]->coord[0][2]));
				l_ = min(l_, norm2(gran->yzels[0]->coord[0][0] - gran->yzels[3]->coord[0][0],
					gran->yzels[0]->coord[0][1] - gran->yzels[3]->coord[0][1],
					gran->yzels[0]->coord[0][2] - gran->yzels[3]->coord[0][2]));

				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += (gran->center[0][i] - P.coord[i]) / 200.0 + normal[i] * l_ / 200.0;
				}
				goto vv1;
			}
			if (kkk2 < 4)
			{
				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += 0.001 * time * P.Vel[i];
				}
				goto vv1;
			}
			else if (kkk2 < 6)
			{
				perturbVector_1(P.coord[0], P.coord[1], P.coord[2], Cell_do->geo_parameters["l_size"] / 100);
				goto vv1;
			}
			else if (kkk2 < 9)
			{
				for (short int i = 0; i < 3; i++)
				{
					P.coord[i] += (-Cell_do->center[0][i] + P.coord[i]) / 100.0;
				}
				goto vv1;
			}

			cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
			cout << "Poteryal B" << endl;
			this->DDD->mut.lock();
			this->DDD->E_esc += P.mu;
			this->DDD->mut.unlock();
			return;
		}


		vtoroy_shans = false;

		continue;
	}
}


/**
 * Розыгрыш нового направления фотона при рассеянии по фазовой функции Хеньи-Гринштейна.
 *
 * @param Vel   исходная скорость (единичный вектор) до рассеяния
 * @param g     параметр асимметрии (-1..1); при |g| < 1e-8 считается изотропным
 * @param ksi1  случайное число из [0,1) для розыгрыша cos?
 * @param ksi2  случайное число из [0,1) для розыгрыша азимута ?
 * @param Vel2  выходной вектор: новое направление (нормированное)
 */
void scatter_1(const Eigen::Vector3d& Vel, double g, double ksi1, double ksi2, Eigen::Vector3d& Vel2)
{
	// 1. Розыгрыш cos? по распределению Хеньи-Гринштейна
	double cos_theta;
	const double eps = 1e-8;
	if (std::fabs(g) < eps) {
		// Изотропное рассеяние
		cos_theta = 2.0 * ksi1 - 1.0;
	}
	else {
		double g2 = g * g;
		double denom = 1.0 - g + 2.0 * g * ksi1;
		double term = (1.0 - g2) / denom;
		cos_theta = (1.0 + g2 - term * term) / (2.0 * g);
		// Защита от ошибок округления
		if (cos_theta > 1.0) cos_theta = 1.0;
		if (cos_theta < -1.0) cos_theta = -1.0;
	}

	double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
	double phi = 2.0 * const_pi * ksi2;
	double cos_phi = std::cos(phi);
	double sin_phi = std::sin(phi);

	// 2. Локальный вектор рассеяния (в системе, где исходное направление — ось Z)
	Eigen::Vector3d local_dir(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);

	// 3. Построение ортонормированного базиса с Z' = Vel
	Eigen::Vector3d Z = Vel.normalized();
	Eigen::Vector3d X, Y;

	// Выбираем вспомогательный вектор, не коллинеарный Z
	if (std::fabs(Z.x()) < 0.9) {
		X = Z.cross(Eigen::Vector3d::UnitX()).normalized();
	}
	else {
		X = Z.cross(Eigen::Vector3d::UnitY()).normalized();
	}
	Y = Z.cross(X).normalized();  // Y ортогонален и X, и Z

	// 4. Преобразование в глобальные координаты
	Vel2 = X * local_dir.x() + Y * local_dir.y() + Z * local_dir.z();
	Vel2.normalize();  // финальная нормировка на случай накопления ошибок
}

void isotropic_direction(double xi1, double xi2, Eigen::Vector3d& Vel2)
{
	double cos_theta = 2.0 * xi1 - 1.0;
	double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
	double phi = 2.0 * const_pi * xi2;
	Vel2[0] = sin_theta * std::cos(phi);
	Vel2[0] = sin_theta * std::sin(phi);
	Vel2[0] = cos_theta;
	return;
}



// Функция, которая отклоняет вектор, сохраняя его длину
void perturbVectorKeepingMagnitude_1(double& Vx, double& Vy, double& Vz, const double& noiseScale) 
{
	// 1. Вычисляем исходную длину
	const double magnitude = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
	if (magnitude == 0.0) return;  // нулевой вектор нельзя отклонить

	// 2. Нормализуем вектор (делаем единичным)
	const double invMag = 1.0 / magnitude;
	Vx *= invMag;
	Vy *= invMag;
	Vz *= invMag;

	// 3. Генерируем случайный перпендикулярный вектор (шум)
	double noiseX = randomNoise_1(noiseScale);
	double noiseY = randomNoise_1(noiseScale);
	double noiseZ = randomNoise_1(noiseScale);

	// 4. Делаем шум строго перпендикулярным исходному вектору (чтобы не менять длину)
	const double dot = Vx * noiseX + Vy * noiseY + Vz * noiseZ;
	noiseX -= dot * Vx;
	noiseY -= dot * Vy;
	noiseZ -= dot * Vz;

	// 5. Добавляем шум и нормализуем
	Vx += noiseX;
	Vy += noiseY;
	Vz += noiseZ;

	const double newMag = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
	const double correction = magnitude / newMag;
	Vx *= correction;
	Vy *= correction;
	Vz *= correction;
}

// Генерация случайного числа в диапазоне [-scale, +scale]
double randomNoise_1(double scale) 
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(-scale, scale);
	return dist(gen);
}

bool findSphereIntersectionTime_1(
	const Eigen::Vector3d& X,  // Положение частицы
	const Eigen::Vector3d& V,  // Скорость частицы
	const double& R,        // Радиус сферы
	double& time)
{
	const double a = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
	const double b = 2.0 * (X[0] * V[0] + X[1] * V[1] + X[2] * V[2]);
	const double c = (X[0] * X[0] + X[1] * X[1] + X[2] * X[2]) - R * R;

	const double D = b * b - 4 * a * c;

	if (D < 0)
	{
		return false;  // Нет пересечений
	}

	const double sqrtD = std::sqrt(D);
	const double t1 = (-b - sqrtD) / (2 * a);
	const double t2 = (-b + sqrtD) / (2 * a);

	// Находим минимальное положительное время
	if (t1 >= 0)
	{
		time = t1;
		return true;
	}
	if (t2 >= 0)
	{
		time = t2;
		return true;
	}

	return false;  // Оба времени отрицательные (пересечение было в прошлом)
}

// Функция, которая отклоняет вектор
void perturbVector_1(double& Vx, double& Vy, double& Vz, double noiseScale) 
{

	// 3. Генерируем случайный перпендикулярный вектор (шум)
	double noiseX = randomNoise_1(noiseScale);
	double noiseY = randomNoise_1(noiseScale);
	double noiseZ = randomNoise_1(noiseScale);

	// 5. Добавляем шум и нормализуем
	Vx += noiseX;
	Vy += noiseY;
	Vz += noiseZ;
}
