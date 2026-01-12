#include "Setka.h"


bool findSphereIntersectionTime(
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


// Генерация случайного числа в диапазоне [-scale, +scale]
double randomNoise(double scale) {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(-scale, scale);
	return dist(gen);
}

// Функция, которая отклоняет вектор, сохраняя его длину
void perturbVectorKeepingMagnitude(double& Vx, double& Vy, double& Vz, const double& noiseScale = 0.01) {
	// 1. Вычисляем исходную длину
	const double magnitude = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
	if (magnitude == 0.0) return;  // нулевой вектор нельзя отклонить

	// 2. Нормализуем вектор (делаем единичным)
	const double invMag = 1.0 / magnitude;
	Vx *= invMag;
	Vy *= invMag;
	Vz *= invMag;

	// 3. Генерируем случайный перпендикулярный вектор (шум)
	double noiseX = randomNoise(noiseScale);
	double noiseY = randomNoise(noiseScale);
	double noiseZ = randomNoise(noiseScale);

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

// Функция, которая отклоняет вектор
void perturbVector(double& Vx, double& Vy, double& Vz, double noiseScale = 0.01) {

	// 3. Генерируем случайный перпендикулярный вектор (шум)
	double noiseX = randomNoise(noiseScale);
	double noiseY = randomNoise(noiseScale);
	double noiseZ = randomNoise(noiseScale);

	// 5. Добавляем шум и нормализуем
	Vx += noiseX;
	Vy += noiseY;
	Vz += noiseZ;
}

void Setka::Print_fH(short int zoneMK, Type_Gran_surf type, const double ex, const double ey, const double ez, const double dphi)
{
	if (this->MK_Grans.size() == 0)
	{
		this->Set_MK_Zone();
	}

	vector<Gran*> Gran_for_print;

	Gran* gr2 = nullptr;
	double s1 = -1.0;
	Eigen::Vector3d e;
	e << ex, ey, ez;

	for (const auto& gr1 : this->MK_Grans[zoneMK - 1])
	{
		if (gr1->type2 != type) continue;

		Eigen::Vector3d v1;
		v1 << gr1->center[0][0], gr1->center[0][1], gr1->center[0][2];
		double d = e.dot(v1) / v1.norm() / e.norm();   // Это косинус угла между гранями
		if (acos(d) <= dphi)
		{
			//Gran_for_print.push_back(gr1);

			//s1 = d;
			//gr2 = gr1;
		}

		if (d > s1)
		{
			s1 = d;
			gr2 = gr1;
		}
	}

	Gran_for_print.push_back(gr2);

	if (Gran_for_print.size() == 0)
	{
		cout << "Ne nashli grans" << endl;
		return;
	}

	for (auto& gr2 : Gran_for_print)
	{
		gr2->Culc_measure(0); // Вычисляем площадь грани (на всякий случай ещё раз)
		gr2->MK_Potok = 0.0;
		// Выделяем место под AMR, сколько сортов водорода, столько и места
		if (gr2->AMR.size() < this->phys_param->num_H)
		{
			gr2->AMR.resize(this->phys_param->num_H);
			for (size_t i = 0; i < this->phys_param->num_H; i++)
			{
				gr2->AMR[i][0] = nullptr;
				gr2->AMR[i][1] = nullptr;
			}
		}

		short int ni = 1;  // Определяем выходящую функцию распределения
		if (gr2->cells[0]->MK_zone == zoneMK) ni = 0;

		short int ni2 = 0; // Определяем входящую функцию распределения
		if (gr2->cells[0]->MK_zone == zoneMK) ni2 = 1;

		// Загружаем все функции распределения на грани
		for (short int ii = 0; ii <= 1; ii++)
		{
			for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
			{
				gr2->Read_AMR(ii, iH, this->phys_param, false);
			}
		}
	}

	// Печатаем функции распределения на грани
	for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
	{
		cout << "Start print AMR:  iH = " << iH << "  grans: " << Gran_for_print.size() << endl;
		Print_AMR(iH, Gran_for_print);
	}


	for (auto& gr2 : Gran_for_print)
	{
		// Удаляем все функции распределения на грани
		for (short int ii = 0; ii <= 1; ii++)
		{
			for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
			{
				gr2->AMR[iH - 1][ii]->Delete();
				delete gr2->AMR[iH - 1][ii];
				gr2->AMR[iH - 1][ii] = nullptr;
			}
		}
		gr2->AMR.clear();
	}


}

void Setka::Print_f_proect_in_gran(short int nn)
{
	// nn == 1  TS
	// nn == 2  HP
	// nn == 3  BS
	Eigen::Vector3d e;
	e << 1.0, 0.0, 0.0;
	Gran* gr2 = nullptr;
	double s1 = -1.0;

	vector<Gran*>* List_Gran;

	if (nn == 1) List_Gran = &this->Gran_TS;
	if (nn == 2) List_Gran = &this->Gran_HP;
	if (nn == 3) List_Gran = &this->Gran_BS;

	for (const auto& gr1 : *List_Gran)
	{
		Eigen::Vector3d v1;
		v1 << gr1->center[0][0], gr1->center[0][1], gr1->center[0][2];
		double d = e.dot(v1) / v1.norm() / e.norm();   // Это косинус угла между гранями
		if (d > s1)
		{
			s1 = d;
			gr2 = gr1;
		}
	}

	cout << "Print_f_proect_in_gran: Gran center: " << gr2->center[0][0] << " " << gr2->center[0][1] << " " << gr2->center[0][2] << endl;

	Cell* A;
	Cell* B;

	A = gr2->cells[0];
	B = gr2->cells[1];

	A->Init_mas_pogl(this->phys_param->pogl_n, this->phys_param->num_H);
	A->read_mas_pogl_FromFile(this->phys_param);

	B->Init_mas_pogl(this->phys_param->pogl_n, this->phys_param->num_H);
	B->read_mas_pogl_FromFile(this->phys_param);

	ofstream fout;
	string name_f;
	if (nn == 1) name_f = "f_proekts_on_TS.txt";
	if (nn == 2) name_f = "f_proekts_on_HP.txt";
	if (nn == 3) name_f = "f_proekts_on_BS.txt";

	fout.open(name_f);
	//fout << "TITLE = HP  VARIABLES = u, f1L, f1R, f1L_fluid, f1Lmoment, f2L, f2R, f2L_fluid, f2Lmoment, f3L, f3R, f3L_fluid, f3Lmoment, f4L, f4R, f4L_fluid, f4Lmoment, ff_inf" << endl;
	fout << "TITLE = HP  VARIABLES = u, " << endl;
	for (int i = 0; i < this->phys_param->num_H; i++)
	{
		fout << "f_" + to_string(i + 1) << "_L, f_" + to_string(i + 1) << "_R, ";
	}

	fout << "ff_inf" << endl;


	double dv = (this->phys_param->pogl_R - this->phys_param->pogl_L) / this->phys_param->pogl_n;


	for (int j = 0; j < this->phys_param->pogl_n; j++)
	{
		double VV = this->phys_param->pogl_L + dv * (j + 0.5);
		fout << VV << " ";

		for (int i = 0; i < this->phys_param->num_H; i++)
		{
			double u1 = A->parameters[0]["Vx_H" + to_string(i + 1)];
			double u2 = A->parameters[0]["Vy_H" + to_string(i + 1)];
			double u3 = A->parameters[0]["Vz_H" + to_string(i + 1)];
			double n = A->parameters[0]["rho_H" + to_string(i + 1)];
			double p = A->parameters[0]["p_H" + to_string(i + 1)];
			double c = sqrt(2.0 * p / n);

			double u1_MK = A->parameters[0]["MK_Vx_H" + to_string(i + 1)];
			double u2_MK = A->parameters[0]["MK_Vy_H" + to_string(i + 1)];
			double u3_MK = A->parameters[0]["MK_Vz_H" + to_string(i + 1)];
			double n_MK = A->parameters[0]["MK_n_H" + to_string(i + 1)];
			double p_MK = A->parameters[0]["MK_T_H" + to_string(i + 1)];
			double c_MK = sqrt(p_MK);

			if (c_MK < 0.00000000001) c_MK = 1.0;

			fout << A->mas_pogl(i, j) * this->phys_param->par_n_H_LISM / dv << " "
				<< B->mas_pogl(i, j) * this->phys_param->par_n_H_LISM / dv << " ";// <<
				//n / (sqrt_pi * c) * exp(-kv(VV - u1) / kv(c)) * this->phys_param->par_n_H_LISM << " " <<
				//n_MK / (sqrt_pi * c_MK) * exp(-kv(VV - u1_MK) / kv(c_MK)) * this->phys_param->par_n_H_LISM << " ";
		}

		fout << 3.0 / (sqrt_pi * 1.0) * exp((-kv(VV - this->phys_param->Velosity_inf)) / kv(1.0));

		fout << endl;
	}

	fout.close();

	A->Delete_mas_pogl();
	B->Delete_mas_pogl();
}

void Setka::Print_f_proect_in_cell(const double& x, const double& y, const double& z)
{
	Eigen::Vector3d e;
	e << 1.0, 0.0, 0.0;
	Gran* gr2 = nullptr;
	double s1 = -1.0;

	Cell* A = nullptr;
	Cell* previos = nullptr;
	
	A = this->Find_cell_point(x, y, z, 0, previos);

	if (A == nullptr) return;

	cout << "Print_f_proect_in_cell: Cell center: " << A->center[0][0] << " " << A->center[0][1] << " " << A->center[0][2] << endl;

	A->Init_mas_pogl(this->phys_param->pogl_n, this->phys_param->num_H);
	A->read_mas_pogl_FromFile(this->phys_param);

	ofstream fout;
	string name_f = to_string(A->number) + "__" + to_string(x) + "_f_proekts_in_Cell.txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = u, f1, f1_fluid, f1_moment, f2, f2_fluid, f2_moment, f3, f3_fluid, f3_moment, f4, f4_fluid, f4_moment" << endl;
	double dv = (this->phys_param->pogl_R - this->phys_param->pogl_L) / this->phys_param->pogl_n;

	for (int j = 0; j < this->phys_param->pogl_n; j++)
	{
		double VV = this->phys_param->pogl_L + dv * (j + 0.5);
		fout << VV << " "; 

			for (int i = 0; i < this->phys_param->num_H; i++)
			{
				double u1 = A->parameters[0]["Vx_H" + to_string(i + 1)];
				double u2 = A->parameters[0]["Vy_H" + to_string(i + 1)];
				double u3 = A->parameters[0]["Vz_H" + to_string(i + 1)];
				double n = A->parameters[0]["rho_H" + to_string(i + 1)];
				double p = A->parameters[0]["p_H" + to_string(i + 1)];
				double c = sqrt(2.0 * p / n);

				double u1_MK = A->parameters[0]["MK_Vx_H" + to_string(i + 1)];
				double u2_MK = A->parameters[0]["MK_Vy_H" + to_string(i + 1)];
				double u3_MK = A->parameters[0]["MK_Vz_H" + to_string(i + 1)];
				double n_MK = A->parameters[0]["MK_n_H" + to_string(i + 1)];
				double p_MK = A->parameters[0]["MK_T_H" + to_string(i + 1)];
				double c_MK = sqrt(p_MK);
				fout << A->mas_pogl(i, j) * this->phys_param->par_n_H_LISM / dv << " " <<
					n / (sqrt_pi * c) * exp(-kv(VV - u1) / kv(c)) * this->phys_param->par_n_H_LISM << " " <<
					n_MK / (sqrt_pi * c_MK) * exp(-kv(VV - u1_MK) / kv(c_MK)) * this->phys_param->par_n_H_LISM << " ";
			}

		fout << endl;
	}

	fout.close();

	A->Delete_mas_pogl();
}

void Setka::Set_MK_Zone(void)
{
	cout << "Start Set_MK_Zone" << endl;
	this->Renumerate();
	Eigen::Vector3d Centr;

	this->Cell_Center->MK_zone_r = 1;
	this->Cell_Center->MK_zone_phi = 0;

	bool new_bound = true;   // Нужно ли подвинуть внешнюю гарницу ближе?

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
				if (new_bound & Centr[0] > this->phys_param->R_MK_Max)
				{
					cell->MK_zone = 8;         // Фиктивная зона (нужна для того, чтобы раздилить зону 6)
				}
				else
				{
					cell->MK_zone = 6;
				}
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
	this->MK_Potoks_on_sort.resize(7);
	for (short int i = 0; i < 7; ++i)
	{
		this->MK_Potoks_on_sort[i].resize(this->phys_param->num_H, 0.0);
	}



	for (short int i = 0; i < 7; ++i) this->MK_Potoks[i] = 0.0;
	this->MK_zone_4.resize(7, 4);
	this->MK_zone_4 <<  false, false, false, false,
						false, false, false, false,
						false, false, false, false,
						false, false, false, false,
						false, false, false, false,
						false, false, false, false,
						false, false, false, false;
	for (auto& cell : this->All_Cell)
	{
		if (cell->MK_zone == 8) continue;  // Фиктивная зона

		int zone = this->determ_zone(cell, 0);

		// Некоторые проверки для правильного определения зоны
		if (cell->MK_zone == 3 && zone == 1)
		{
			cout << "INFORM teegedcvfv" << endl;
			cout << cell->center[0][0] << " " << cell->center[0][1] << " " << cell->center[0][2] << endl;
		}


		this->MK_zone_4(cell->MK_zone - 1, zone - 1) = true;
	}

	this->H_komponent_in_zone.resize(4, this->phys_param->num_H);
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < this->phys_param->num_H; j++)
		{
			this->H_komponent_in_zone(i, j) = false;
		}
	}

	for (int i = 0; i < this->phys_param->hydrogen_arise_1.rows(); ++i) 
	{
		for (int j = 0; j < this->phys_param->hydrogen_arise_1.cols(); ++j) 
		{
			int8_t value = this->phys_param->hydrogen_arise_1(i, j);
			this->H_komponent_in_zone(0, value - 1) = true;
		}
	}

	for (int i = 0; i < this->phys_param->hydrogen_arise_2.rows(); ++i)
	{
		for (int j = 0; j < this->phys_param->hydrogen_arise_2.cols(); ++j)
		{
			int8_t value = this->phys_param->hydrogen_arise_2(i, j);
			this->H_komponent_in_zone(1, value - 1) = true;
		}
	}

	for (int i = 0; i < this->phys_param->hydrogen_arise_3.rows(); ++i)
	{
		for (int j = 0; j < this->phys_param->hydrogen_arise_3.cols(); ++j)
		{
			int8_t value = this->phys_param->hydrogen_arise_3(i, j);
			this->H_komponent_in_zone(2, value - 1) = true;
		}
	}

	for (int i = 0; i < this->phys_param->hydrogen_arise_4.rows(); ++i)
	{
		for (int j = 0; j < this->phys_param->hydrogen_arise_4.cols(); ++j)
		{
			int8_t value = this->phys_param->hydrogen_arise_4(i, j);
			this->H_komponent_in_zone(3, value - 1) = true;
		}
	}

	this->MK_zone_H.resize(7, this->phys_param->num_H);
	for (size_t i = 0; i < 7; i++)
	{
		for (size_t j = 0; j < this->phys_param->num_H; j++)
		{
			this->MK_zone_H(i, j) = false;
		}
	}

	for (size_t i = 0; i < 7; i++)
	{
		for (size_t j = 0; j < this->phys_param->num_H; j++)
		{
			bool is_H = false;
			for (size_t zone_ = 0; zone_ < 4; zone_++)
			{
				if (this->MK_zone_4(i, zone_) == true)
				{
					if (this->H_komponent_in_zone(zone_, j) == true)
					{
						is_H = true;
					}
				}
			}
			this->MK_zone_H(i, j) = is_H;
		}
	}


	// Посмотрим что получилось
	if (true)
	{
		cout << "See MK_zone_4" << endl;
		for (size_t i = 0; i < 7; i++)
		{
			for (size_t j = 0; j < 4; j++)
			{
				cout << this->MK_zone_4(i, j) << " ";
			}
			cout << endl;
		}
		cout << endl;

		cout << "See H_komponent_in_zone" << endl;

		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < this->phys_param->num_H; j++)
			{
				cout << this->H_komponent_in_zone(i, j) << " ";
			}
			cout << endl;
		}
		cout << endl;

		cout << "See MK_zone_H" << endl;

		for (size_t i = 0; i < 7; i++)
		{
			for (size_t j = 0; j < this->phys_param->num_H; j++)
			{
				cout << this->MK_zone_H(i, j) << " ";
			}
			cout << endl;
		}
		cout << endl;
	}


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
	// старый вариант до внешней границы
	if (!new_bound)
	{
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
	}
	else
	{
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

			if (gr->type == Type_Gran::Outer_Hard && gr->cells[0]->center[0][0] <= this->phys_param->R_MK_Max)
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

			// Это те самые добавленные грани, которые позволяют сократить область моделирования
			if (gr->cells.size() == 2)
			{
				if (gr->cells[0]->type == Type_cell::Zone_4 &&
					gr->cells[1]->type == Type_cell::Zone_4)
				{
					if ( (gr->cells[0]->center[0][0] - this->phys_param->R_MK_Max) * 
						(gr->cells[1]->center[0][0] - this->phys_param->R_MK_Max) < 0.0)
					{
						this->MK_Grans[5].push_back(gr);
						gr->MK_type.push_back(6);
						gr->type = Type_Gran::Outer_Hard;
						if (gr->cells[0]->MK_zone == 8)
						{
							gr->normal[0][0] = -gr->normal[0][0];
							gr->normal[0][1] = -gr->normal[0][1];
							gr->normal[0][2] = -gr->normal[0][2];
							auto dfd = gr->cells[1];
							gr->cells[1] = gr->cells[0];
							gr->cells[0] = dfd;
							
							//gr->cells.clear();
							//gr->cells.resize(1);
							//gr->cells[0] = dfd;

							// Здесь нужно полностью сымитировать граничную грань, поэтому надо сделать ей только одного соседа.
							// При этом газовую динамику считать больше нельзя, нужно заново строить сетку.
							//cout << "Num ef = " << gr->number << endl;
							//exit(-1);
						}
					}
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
		unsigned int N_ = 0;
		for (auto& gr : this->MK_Grans[jj])
		{
			if (gr->type == Type_Gran::Us) N_++;
		}
		cout << "Only active:  " << N_ << endl;

	}

	// Проверяем, что в массивах нет повторов
	if (true)
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
					cout << i->center[0][0] << " " << i->center[0][1] << " " << i->center[0][2] << endl;
					exit(-1);
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
		double a1, b1, c;
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
		cout << "Error 2341963846 " << this->MK_Grans.size() << endl;
		exit(-1);
	}

	// Готовим/загружаем AMR сетку для граней
	if (true)
	{
		cout << "Start: Zagruzka AMR" << endl;
		unsigned int N1[20];
		unsigned int N_vxod[20];
		for (short int i = 0; i < 20; i++)
		{
			N1[i] = 0;
			N_vxod[i] = 0;
		}

		unsigned short int NNall = 0;
		double S = 0.0;
		vector<double> SS(this->phys_param->num_H);
		for (short int iH = 0; iH < this->phys_param->num_H; iH++)
		{
			SS[iH] = 0.0;
		}

		unsigned int k1 = 0;

#pragma omp parallel for schedule(dynamic)
		for(size_t ijk = 0; ijk < this->MK_Grans[zone_MK - 1].size(); ijk++)
		{
			#pragma omp critical (first) 
			{
				k1++;
				if (k1 % 500 == 0)
				{
					cout << "Gran = " << k1 << "    Iz: " << this->MK_Grans[zone_MK - 1].size() << endl;
				}
			}

			auto gr = this->MK_Grans[zone_MK - 1][ijk];
			//if (gr->number != 1537) continue;                                                          // !DELETE
			//cout << "Center = " << gr->center[0][0] << " " <<
			//	gr->center[0][1] << " " << gr->center[0][2] << endl;                                  // !DELETE
			//gr->Culc_measure(0); // Вычисляем площадь грани (на всякий случай ещё раз)
			gr->MK_Potok = 0.0;

			// Выделяем место под AMR, сколько сортов водорода, столько и места
			if (gr->AMR.size() < this->phys_param->num_H)
			{
				gr->AMR.resize(this->phys_param->num_H);
				gr->MK_Potok_on_sort.resize(this->phys_param->num_H);
				for (size_t i = 0; i < this->phys_param->num_H; i++)
				{
					gr->AMR[i][0] = nullptr;
					gr->AMR[i][1] = nullptr;
					gr->MK_Potok_on_sort[i] = 0.0;
				}
			}

			short int ni = 1;  // Определяем выходящую функцию распределения
			if (gr->cells[0]->MK_zone == zone_MK) ni = 0;

			short int ni2 = 0; // Определяем входящую функцию распределения
			if (gr->cells[0]->MK_zone == zone_MK) ni2 = 1;


			//cout << "ni1, ni2 = " << ni << " " << ni2 << endl;                                                            // !DELETE
			//cout << "sosed = " << gr->cells.size() << endl;                                                            // !DELETE
			//cout << "Type_Gran = " << int(gr->type) << endl;                                                            // !DELETE


			for(short int ii = 0; ii <= 1; ii++)
			{
				for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
				{
					// Выходящую функцию загружаем, обнуляем и сохраняем
					if (false)//(ni == ii)
					{
						if (gr->type == Type_Gran::Us)
						{
							if (gr->AMR[iH - 1][ii] != nullptr)
							{
								cout << "ERROR ergiegkjoeigjergeg" << endl;
								exit(-1);
							}
							gr->Read_AMR(ii, iH, this->phys_param, this->phys_param->refine_AMR);
							gr->AMR[iH - 1][ni]->Fill_null();

							#pragma omp critical (vixod) 
							{
								//N_vixod[iH - 1] += gr->AMR[iH - 1][ni]->Size();
								//N1[iH - 1]++;
							}

							string name_f = "func_grans_AMR_" + to_string(ii) + "_H" +
								to_string(iH) + "_" + to_string(gr->number) + ".bin";

							if (this->phys_param->save_AMR == true)
							{
								gr->AMR[iH - 1][ii]->Save(this->phys_param->AMR_folder + "/" + name_f);
							}

							gr->AMR[iH - 1][ii]->Delete();
							delete gr->AMR[iH - 1][ii];
							gr->AMR[iH - 1][ii] = nullptr;
						}
					}

					if (ni2 == ii)
					{
						if (gr->AMR[iH - 1][ii] != nullptr)
						{
							cout << "ERROR ergeyh45t34tergdgrdsrg" << endl;
							cout << "iH = " << iH << "   ii" << ii << endl;
							exit(-1);
						}
						gr->Read_AMR(ii, iH, this->phys_param, false);

						/*if (iH == 4 && gr->type2 == Type_Gran_surf::BS && kvv(gr->center[0][1], gr->center[0][2], 0.0) < 200.0)
						{
							gr->AMR[iH - 1][ii]->Re_Partially_free_space();
							gr->AMR[iH - 1][ii]->Analyze_memory_usage(true);
							exit(-1);
						}*/

						if (gr->type == Type_Gran::Us)
						{
							#pragma omp critical (vxod) 
							{
								N_vxod[iH - 1] += gr->AMR[iH - 1][ii]->Size();
								N1[iH - 1]++;
							}
						}
					}

				}
			}


			// Считаем сразу поток атомов через грань (только для входящей функции)
			if (true)
			{
				if (gr->type == Type_Gran::Us)
				{
					short int sort_H = -1;
					for (auto& ai : gr->AMR)
					{
						sort_H++;
						ai[ni2]->Culk_SpotokV(gr->area[0]);

						#pragma omp critical (erfgwerweS) 
						{
							S += ai[ni2]->SpotokV;
							SS[sort_H] += ai[ni2]->SpotokV;
						}
						gr->MK_Potok += ai[ni2]->SpotokV;
						gr->MK_Potok_on_sort[sort_H] += ai[ni2]->SpotokV;
					}
				}
				else // Вручную посчитаем поток с границы расчётной области
				{
					for (size_t ijk = 0; ijk < this->phys_param->num_H; ijk++)
					{
						gr->AMR[ijk][ni2]->SpotokV = 0.0;
					}

					Eigen::Vector3d n;
					n << -gr->AMR[3][ni2]->Vn[0], -gr->AMR[3][ni2]->Vn[1], -gr->AMR[3][ni2]->Vn[2];
					// Так как нормаль должна быть внешняя к грани (а надо внутреннюю)

					//cout << "POTOK = " << Get_Spotok_inf(n) << " " << Get_Spotok_inf(-n) << endl;                                        // !DELETE
					double sjv = Get_Spotok_inf(n);
					gr->AMR[3][ni2]->SpotokV = sjv * gr->area[0];
					#pragma omp critical (erfgwerweS) 
					{
						S += gr->AMR[3][ni2]->SpotokV;
						SS[3] += gr->AMR[3][ni2]->SpotokV;
					}
					gr->MK_Potok += gr->AMR[3][ni2]->SpotokV;
					gr->MK_Potok_on_sort[3] += gr->AMR[3][ni2]->SpotokV;
				}
			}

			// Удалим сразу все входящие функции распределения (так как цель была посчитать поток)
			for (short int iH = 1; iH <= gr->AMR.size(); iH++)
			{
				gr->AMR[iH - 1][ni2]->Delete();
				delete gr->AMR[iH - 1][ni2];
				gr->AMR[iH - 1][ni2] = nullptr;
			}

		}
		//exit(-1);                                                          // !DELETE
		//cout << "Izmelcheno  " << NNall << "  yacheek" << endl;
		std::ofstream file1("info_AMR_size.txt", std::ios::app);
		for (short int iH = 0; iH < 9; iH++)
		{
			if (N1[iH] == 0) N1[iH] = 1;
			// Если надо писать в файл информацию про входные грани (так как они не меняются, я решил не писать)
			if (false)
			{
				file1 << "Zone:  " << zone_MK << "   H = " << iH + 1 <<
					"   vxod size: " << 1.0 * N_vxod[iH] / N1[iH] << std::endl;
			}
		}
		file1.close();
		cout << "End: Zagruzka AMR" << endl;
		this->MK_Potoks[zone_MK - 1] = S; // Входящий поток через всю границу зоны
		for (short int iH = 0; iH < this->phys_param->num_H; iH++)
		{
			this->MK_Potoks_on_sort[zone_MK - 1][iH] = SS[iH]; // Входящий поток через всю границу зоны отдельно для каждого сорта
		}
	}

	//pause_seconds(15);

	// Можно удалить файлы граничных граней, так как они только место занимают
	// Уже удалил, больше не должны появляться
	if (false)
	{
		cout << "Start: Ydalenie lishnix failov" << endl;
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			if (gr->type != Type_Gran::Us)
			{
				for (short int ii = 0; ii <= 1; ii++)
				{
					for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
					{
						if (gr->AMR[iH - 1][ii] == nullptr)
						{
							string name_f = "func_grans_AMR_" + to_string(ii) + "_H" +
								to_string(iH) + "_" + to_string(gr->number) + ".bin";
							if (file_exists(this->phys_param->AMR_folder + "/" + name_f))
							{
								std::filesystem::remove(name_f);
							}
						}
					}
				}
			}
		}
		cout << "END: Ydalenie lishnix failov" << endl;
	}

	// Заполняем граничные условия для AMR сетки (на границе расчётной области)
	// В том случае, если среди граней зоны есть граничные грани, на которых надо задавать максвелл
	// можно не задавать, так как я сделал аналитическое разыгрывание скорости
	if (false)
	{
		//for (auto& gr : this->All_boundary_Gran)
#pragma omp parallel for schedule(dynamic)
		//for (auto& gr : this->MK_Grans[zone_MK - 1])
		for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
		{
			auto& gr = this->MK_Grans[zone_MK - 1][idx];
			if (gr->type == Type_Gran::Outer_Hard || gr->type == Type_Gran::Outer_Soft)
			//if(true)  // тестирование, задаём максвел для всех граней
			//if (gr->type2 == Type_Gran_surf::BS)
			{
				if (gr->AMR.size() != 0)
				{
					short int ni = 0;  // определяем входящую функцию распределения
					if (gr->cells[0]->MK_zone == zone_MK)
					{
						ni = 1;
					}

					/*gr->AMR[3][ni]->Print_all_center_Tecplot(gr->AMR[3][ni], "1");
					gr->AMR[3][ni]->Culk_SpotokV(gr->area[0]);
					cout << "1 POTOK = " << gr->AMR[3][ni]->SpotokV/ gr->area[0] << endl;
					cout << gr->AMR[3][ni]->Vn[0] << " " << gr->AMR[3][ni]->Vn[1] << "  " <<
						gr->AMR[3][ni]->Vn[2] << endl;*/

					// Здесь надо задать граничные условия для четвёртого сорта 
					// а также помельчить сетку, если необходимо
					//cout << "_----------------------_" << endl;
					gr->AMR[3][ni]->Fill_maxwel_inf(this->phys_param->Velosity_inf);
					/*cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->de_Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->de_Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					cout << gr->AMR[3][ni]->de_Refine() << endl;
					cout << "Size = " << gr->AMR[3][ni]->Size() << endl;
					exit(-1);*/



					// <<1>> - внутрь, <<3>> - 4 сорт водорода
					//gr->AMR[3][ni]->Print_all_center_Tecplot(gr->AMR[3][ni], "2");
					//gr->AMR[3][ni]->Culk_SpotokV(gr->area[0]);
					//cout << "2 POTOK = " << gr->AMR[3][ni]->SpotokV/ gr->area[0] << endl;
					//exit(-1);
				}
			}
		}
	}

	// Теперь для каждой функции распределения вычисляем поток через неё
	if (false)
	{
		double S = 0.0;
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			//gr->Culc_measure(0); // Вычисляем площадь грани (на всякий случай ещё раз)
			// Нужно вычислять поток только у входящей части функции распределения
			short int ni = 0; // Определяем входящую функцию
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 1;
			}
			gr->MK_Potok = 0.0;

			if (gr->type == Type_Gran::Us)
			{
				for (auto& ai : gr->AMR)
				{
					ai[ni]->Culk_SpotokV(gr->area[0]);

					S += ai[ni]->SpotokV;
					gr->MK_Potok += ai[ni]->SpotokV;
				}
			}
			else // Вручную посчитаем поток с границы расчётной области
			{
				gr->AMR[0][ni]->SpotokV = 0.0;
				gr->AMR[1][ni]->SpotokV = 0.0;
				gr->AMR[2][ni]->SpotokV = 0.0;
				Eigen::Vector3d n;
				n << -gr->AMR[3][ni]->Vn[0], -gr->AMR[3][ni]->Vn[1], -gr->AMR[3][ni]->Vn[2];
				// Так как нормаль должна быть внешняя к грани
				double sjv = Get_Spotok_inf(n);
				gr->AMR[3][ni]->SpotokV = sjv * gr->area[0];
				S += gr->AMR[3][ni]->SpotokV;
				gr->MK_Potok += gr->AMR[3][ni]->SpotokV;

				//cout << " Potok = " << sjv << endl;
				//exit(-1);
			}


		}
		this->MK_Potoks[zone_MK - 1] = S; // Входящий поток через всю границу зоны
	}

	// Считаем необходимые геометрические параметры для МК
	// И добавляем переменные в ячейки (заполняем нулями)
	if (true)
	{
		for (auto& i : this->All_Cell)
		{
			i->Set_Cell_Geo_for_MK();
			if (i->MK_zone == zone_MK && this->phys_param->culc_cell_moments == true)
			//if(true)
			{
				for (const auto& nam : this->phys_param->MK_param)
				{
					i->parameters[0][nam] = 0.0;
				}
			}
		}

		for (auto& i : this->All_Gran)
		{
			i->Set_Gran_Geo_for_MK();
		}
	}

	// Обнулим функции распределения, в которые будем накапливать информацию
	if (false)
	{
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			// gr->Culc_measure(0); // Вычисляем площадь грани (на всякий случай ещё раз)
			short int ni = 1;
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 0;
			}

			for (short int iH = 1; iH <= this->phys_param->num_H; iH++)
			{
				auto funk = gr->AMR[iH - 1][ni];

				// Это для того, чтобы помельчить функцию как надо
				/*if (iH == 4)
				{
					funk->Fill_maxwel_inf(this->phys_param->Velosity_inf);
				}*/

				funk->Fill_null();
			}
			
		}
	}


	// Считаем и создадим S+  S-  в ячейках
	if (this->phys_param->MK_source_S == true)
	{
		for (auto& i : this->All_Cell)
		{
			if (i->MK_zone == zone_MK)
			{
				i->Init_S(this->phys_param->num_pui, this->phys_param->pui_nW);
			}
		}
	}

	// Считаем и создадим массивы для поглощения  в ячейках
	if (this->phys_param->culc_pogl == true)
	{
		for (auto& i : this->All_Cell)
		{
			if (i->MK_zone == zone_MK)
			{
				i->Init_mas_pogl(this->phys_param->pogl_n, this->phys_param->num_H);
			}
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

	// Мельчим AMR сетку, если надо
	if (false)
	{
		if (this->phys_param->de_refine_AMR == true)
		{
			cout << "start de_refine_AMR" << endl;
			unsigned int nmnm = 0;
			unsigned int num_ = 0;
			unsigned int sr_num = 0;
#pragma omp parallel for schedule(dynamic)
			for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
			{
				auto& gr = this->MK_Grans[zone_MK - 1][idx];
				if (gr->type != Type_Gran::Us) continue;

				short int ni = 1;
				if (gr->cells[0]->MK_zone == zone_MK)
				{
					ni = 0;
				}

				for (short int iH = 1; iH <= gr->AMR.size(); iH++)
				{
					int iki = gr->AMR[iH - 1][ni]->de_Refine(iH);

#pragma omp critical (third) 
					{
						nmnm += iki;
						num_++;
						sr_num += gr->AMR[iH - 1][ni]->Size();
					}

				}
			}
			cout << "Ydaleno  " << nmnm << "  yacheek" << endl;
			cout << "Srednee chislo yacheek =  " << (1.0 * sr_num) / num_ << endl;
		}
	}

	// Записываем AMR сетку для граней
	if (false)
	{
		// Нужно сохрянять только выходящие грани!
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			short int ni = 1;  // Выходящая функция распределения
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 0;
			}

			for (short int ii = 0; ii <= 1; ii++) 
			{
				for (short int iH = 1; iH <= gr->AMR.size(); iH++)
				{
					if (ii == ni)
					{
						string name_f = "func_grans_AMR_" + to_string(ii) + "_H" +
							to_string(iH) + "_" + to_string(gr->number) + ".bin";

						if (this->phys_param->save_AMR == true)
						{
							gr->AMR[iH - 1][ii]->Save(this->phys_param->AMR_folder + "/" + name_f);
						}

						gr->AMR[iH - 1][ii]->Delete();
						delete gr->AMR[iH - 1][ii];
					}
					else
					{
						gr->AMR[iH - 1][ii]->Delete();
						delete gr->AMR[iH - 1][ii];
					}
				}
			}
			gr->AMR.clear();
			gr->AMR.resize(0);
		}

		//this->MK_Grans.clear();  // Не надо удалять, чтобы сразу считать следующую зону
	}

	// Записываем моменты
	if (this->phys_param->culc_cell_moments == true)
	{
		// Так как в файле хранятся параметры во всех ячейках (даже в тех, которые были за пределом рассчитанной зоны)
		// Нужно сначала скачать все кроме текущей зоны (так как они только что посчитаны), а потом записать все
		if (file_exists(this->phys_param->MK_file))
		{
			this->Download_cell_MK_parameters(this->phys_param->MK_file, zone_MK);
		}
		this->Save_cell_MK_parameters(this->phys_param->MK_file);
	}

	// Считаем и создадим S+  S-  в ячейках
	if (this->phys_param->MK_source_S == true)
	{
		for (auto& i : this->All_Cell)
		{
			if (i->MK_zone == zone_MK)
			{
				i->write_S_ToFile();
				// Освобождаем память
				i->pui_Sm.resize(0);
				i->pui_Sp.resize(0, 0);
			}
		}
	}

	if (this->phys_param->culc_pogl == true)
	{
		for (auto& i : this->All_Cell)
		{
			if (i->MK_zone == zone_MK)
			{
				i->write_mas_pogl_ToFile(this->phys_param);
				i->Delete_mas_pogl();
			}
		}
	}

	cout << "END MK_delete" << endl;
}

void Setka::MK_go(short int zone_MK, int N_per_gran, Interpol* Interpol, Setka*& S_main)
{
	// Как работает алгоритм пошагово:
	// 1) Загружаются выходящие функции распределения только для сортов, которые рождаются в данной области
	// 2) Далее следует цикл по сортам водорода
	//       это не очень удобно для понимания, но позволяет экономить память в разы
	// {
	//    3) Загружаются выходящие функции распределения для данного сорта водорода для всех граней области
    //          в случае, если они ещё не загружены в шаге 1.
	//    4) Далее следует цикл по всем граням области
	//    {
	//       5) Для каждой грани загружается входящая функция распределения для текущего сорта
	//       6) С этой грани запускается нужное число частиц
	//		 7) Входящая функция распределения для данной грани и для данного сорта удаляется из памяти
	//    }
	//    8) Сохраняются в файл и удаляются ненужные выходящие функции распределения (предварительно нормируются):
	//          это функции, загруженные на шаге 3, то есть для тех сортов, которые не рождаются в области.
	// }
	// 9) Нормируются оставшиеся функции распределения (загруженные на шаге 1), сохраняются в файл и удаляются
	// 10) Нормируются моменты в ячейках


	auto start = std::chrono::high_resolution_clock::now();
	cout << "Start MK_go " << zone_MK << "   N_on_gran = " << N_per_gran << endl;
	int N_on_gran = N_per_gran;   // Сколько запускаем частиц на грань в среднем
	
	double mu_expect = 0.0;
	mu_expect = this->MK_Potoks[zone_MK - 1] / 
		(1.0 * N_on_gran * this->MK_Grans[zone_MK - 1].size());

	vector<double>mu_expect_per_sort(this->phys_param->num_H);
	for (short int i = 0; i < this->phys_param->num_H; i++)
	{
		mu_expect_per_sort[i] = max(this->MK_Potoks_on_sort[zone_MK - 1][i], this->MK_Potoks[zone_MK - 1]/1000.0) /
			(1.0 * N_on_gran * this->MK_Grans[zone_MK - 1].size());
	}

	cout << "All potok = " << this->MK_Potoks[zone_MK - 1] << endl;
	//exit(-1);

	unsigned int N2[9];
	unsigned int N_vixod[9];
	for (short int i = 0; i < 9; i++)
	{
		N2[i] = 0;
		N_vixod[i] = 0;
	}

	// Подготовка нужных массивов 
	if (true)
	{
		cout << "Start download_1" << endl;
		// 1. Надо загрузить выходящие функции распределения только для сортов, которые рождаются в данной области
		for (size_t j = 0; j < this->phys_param->num_H; j++)
		{
			cout << "For sort: " << j + 1;
			if (this->MK_zone_H(zone_MK - 1, j) == true)
			{
				cout << "  Nado" << endl;
				// Значит нам надо загрузить выходящую функцию для водорода сорта j
				for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
				{
					auto& gr = this->MK_Grans[zone_MK - 1][idx];
					if (gr->type != Type_Gran::Us) continue;

					short int ni = 1; // Номер "выходящей" функции распределения
					if (gr->cells[0]->MK_zone == zone_MK)
					{
						ni = 0;
					}
					gr->Read_AMR(ni, j + 1, this->phys_param, this->phys_param->refine_AMR);
					gr->AMR[j][ni]->Fill_null();
					N_vixod[j] += gr->AMR[j][ni]->Size();
					N2[j]++;

					gr->AMR[j][ni]->Partially_free_space();
				}
			}
			else
			{
				cout << "  Ne nado" << endl;
			}
		}

		cout << "End download_1" << endl;
	}

	unsigned int ALL_N = 0;  // Общее число запущенных в итоге частиц
	unsigned int k1 = 0;

	// 2. Разыгрываем каждый сорт отдельно, так как для него нужны свои массивы
	for (short int nh_ = 0; nh_ < this->phys_param->num_H; ++nh_)
	//for (short int nh_ = 3; nh_ <= 3; ++nh_)                                                // DELETE
	{
		// Для каждого запускаемого сорта надо загрузить выходяющии функции распределения на всех гранях
		// и входящую функуию только для текущей грани

		cout << "Start download_2  for sort " << nh_ + 1 << endl;
		// 3. Загружаем выходящие функции распределения данного сорта аодорода для всех граней 
		// если они ещё не загружены на предыдущем шаге
		if (this->MK_zone_H(zone_MK - 1, nh_) == false)
		{
			for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
			{
				auto& gr = this->MK_Grans[zone_MK - 1][idx];
				if (gr->type != Type_Gran::Us) continue;

				short int ni = 1; // Номер "выходящей" функции распределения
				if (gr->cells[0]->MK_zone == zone_MK)
				{
					ni = 0;
				}
				gr->Read_AMR(ni, nh_ + 1, this->phys_param, this->phys_param->refine_AMR);
				gr->AMR[nh_][ni]->Fill_null();
				N_vixod[nh_] += gr->AMR[nh_][ni]->Size();
				N2[nh_]++;

				gr->AMR[nh_][ni]->Partially_free_space();
			}
		}
		cout << "End download_2  for sort " << nh_ + 1 << endl;

		// 4. Теперь бежим по граням и делаем основной алгоритм
		k1 = 0;
		#pragma omp parallel for schedule(dynamic)                                                   // DELETE
		for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
		//for (size_t idx = 2700; idx < 2701; ++idx)
		{
			auto& gr = this->MK_Grans[zone_MK - 1][idx];
			
			Eigen::Vector3d n;
			Eigen::Vector3d t;
			Eigen::Vector3d m;

			#pragma omp critical (first) 
			{
				k1++;
				if (k1 % 500 == 0 || k1 == 100)
				{
					cout << "Gran = " << k1 << "    Iz: " << this->MK_Grans[zone_MK - 1].size() << "  sort " << nh_ + 1 << endl;
				}
			}

			//if (gr->type2 != Type_Gran_surf::BS) continue;                             // DELETE
			//if (gr->type != Type_Gran::Outer_Hard) continue;                             // DELETE


			// Выбираем конкретный номер датчика случайных чисел
			unsigned int sens_num1 = 2 * omp_get_thread_num();
			unsigned int sens_num2 = 2 * omp_get_thread_num() + 1;

			short int ni = 0; // Номер "входящей" функции распределения
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 1;
			}
			
			//double full_gran_potok = gr->MK_Potok;
			double full_gran_potok = gr->MK_Potok_on_sort[nh_];

			// Критерий запуска частиц с грани
			// Если полный поток через грань очень мал, то можно не разыгрывать
			//if (full_gran_potok < 0.0000001 * MK_Potoks[zone_MK - 1] / this->MK_Grans[zone_MK - 1].size())
			if (this->MK_Potoks_on_sort[zone_MK - 1][nh_] < 0.000001 || 
				(full_gran_potok < 0.0001 * this->MK_Potoks_on_sort[zone_MK - 1][nh_] / this->MK_Grans[zone_MK - 1].size()))
			{
				continue;
			}

			// Получаем нормаль для граничных граней (это будет внешняя нормаль)
			if (gr->type != Type_Gran::Us)
			{
				n << gr->normal[0][0], gr->normal[0][1], gr->normal[0][2];
				get_bazis(n, t, m);
			}

			auto& func = gr->AMR[nh_][ni];

			if (gr->type == Type_Gran::Us)
			{
				gr->Read_AMR(ni, nh_ + 1, this->phys_param, false);
				func->Culk_SpotokV(gr->area[0]);
				func->Culc_gradients();                   // Считаем градиенты для второго порядка (minmod)
			}
			else
			{
				gr->Read_AMR(ni, nh_ + 1, this->phys_param, false);
				func->SpotokV = 0.0;
				if (nh_ == 3) // Так как поток есть только у атомов 4-го сорта
				{
					Eigen::Vector3d nn;
					nn << -func->Vn[0], -func->Vn[1], -func->Vn[2];
					// Так как нормаль должна быть внешняя к грани
					double sjv = Get_Spotok_inf(nn);
					func->SpotokV = sjv * gr->area[0];
				}
			}



			if (func->SpotokV < 0.0000001 * full_gran_potok)
			{
				// Функция распределния нулевая, можно не разыгрывать
				func->Delete();
				delete func;
				func = nullptr;
				continue;
			}

			// Расчитываем число запускаемых частиц
			
			//unsigned int N_particle = max(static_cast<int>(func->SpotokV / mu_expect) + 1,
			//	min(N_on_gran, 1000));


			unsigned int N_particle = max(static_cast<int>(func->SpotokV / mu_expect_per_sort[nh_]) + 1,
				min(N_on_gran, 1000));

			// Мало ли какой там поток, всё-равно больше планируемого числа запускать не надо
			if (N_particle > N_on_gran) N_particle = N_on_gran;



			double mu = func->SpotokV / N_particle; // Вес каждой частицы

			#pragma omp critical (second) 
			{
				ALL_N += N_particle;
			}

			// Запускаем каждую частицу
			for (unsigned int num = 0; num < N_particle; ++num)                                                  // DELETE
			//for (unsigned int num = 0; num < 1; ++num)                                                  // DELETE
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

				if (P.cel->MK_zone != zone_MK)
				{
					cout << "Error 9767653421" << endl;
					exit(-1);
				}

				P.mu = mu;                           // Вес частицы
				P.sort = nh_ + 1;                    // Сорт частицы

				Eigen::Vector3d poz;

				// Находим положение точки на грани
				gr->Get_Random_pozition(poz, this->Sensors[sens_num1]);
				P.Addcoord(poz);

				// Находим скорость частицы
				if (gr->type != Type_Gran::Us)
				{
					// Можно разыгрывать аналитическую функцию распределения на границе, а не табличную
					this->Velosity_initial(this->Sensors[sens_num1], poz, n, t, m);
				}
				else
				{
					func->Get_random_velosity(func, gr->area[0], poz, this->Sensors[sens_num1]);
				}

				P.AddVel(poz);

				// Некоторые проверки разыгрынной скорости частицы
				if (this->regim_otladki)
				{
					if (std::isnan(poz[0])) 
					{
						std::cout << "ERROR erg34t34t34f " << std::endl;
						cout << poz[0] << " " << poz[1] << " " << poz[2] << endl;
						exit(-1);
					}


					if (ni == 0)
					{
						if (scalarProductFast(poz(0), poz(1), poz(2),
							gr->normal[0][0], gr->normal[0][1], gr->normal[0][2]) < 0.0)
						{
							cout << "Error  8765656431" << endl;
							whach((int)gr->type);
							whach(gr->number);
							whach(gr->normal[0][0]);
							whach(gr->normal[0][1]);
							whach(gr->normal[0][2]);
							whach(poz(0));
							whach(poz(1));
							whach(poz(2));
							whach(nh_);
							whach(num);
							whach(func->Vn[0]);
							whach(func->Vn[1]);
							whach(func->Vn[2]);
							exit(-1);
						}
					}
					else
					{
						if (scalarProductFast(poz(0), poz(1), poz(2),
							gr->normal[0][0], gr->normal[0][1], gr->normal[0][2]) > 0.0)
						{
							cout << "Error  7411100090" << endl;
							whach(gr->number);
							whach(scalarProductFast(poz(0), poz(1), poz(2),
								gr->normal[0][0], gr->normal[0][1], gr->normal[0][2]));
							whach(poz(0));
							whach(poz(1));
							whach(poz(2));
							whach(nh_);
							whach(num);
							whach(gr->normal[0][0]);
							whach(gr->normal[0][1]);
							whach(gr->normal[0][2]);
							whach(gr->center[0][0]);
							whach(gr->center[0][1]);
							whach(gr->center[0][2]);
							for (int i = 0; i < gr->MK_type.size(); i++)
							{
								whach(i);
								whach(gr->MK_type[i]);
							}
							whach(gr->cells[0]->MK_zone);
							whach(gr->cells[1]->MK_zone);
							whach(func->Vn[0]);
							whach(func->Vn[1]);
							whach(func->Vn[2]);
							exit(-1);
						}
					}
				}

				// Надо проверить, (1) что точка находится в нужной ячейке и
				// (2) что она будет находиться в ней через время dt
				if (true)
				{
					Cell* previos = P.cel;
					Eigen::Vector3d Center_cell, Move;
					Center_cell << P.cel->center[0][0], P.cel->center[0][1], P.cel->center[0][2];
					double dt = previos->geo_parameters["l_size"] / P.Vel_norm() / 1000.0;

					Cell* ppp = this->Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, previos);

					short unsigned int klk = 0;
					while (ppp == nullptr || ppp != P.cel)
					{
						klk++;
						if (klk > 100)
						{
							cout << "Error 6439011209" << endl;
							cout << P.coord[0] << " " << P.coord[1] <<
								" " << P.coord[2] << endl;
							exit(-1);
						}

						if (klk > 97)
						{
							cout << "INFO! 457658u56ytg456347" << endl;
							P.coord[0] = Center_cell[0];
							P.coord[1] = Center_cell[1];
							P.coord[2] = Center_cell[2];
						}

						Move[0] = (-P.coord[0] + Center_cell[0]) / 1000.0;
						Move[1] = (-P.coord[1] + Center_cell[1]) / 1000.0;
						Move[2] = (-P.coord[2] + Center_cell[2]) / 1000.0;
						P.Move(Move);
						ppp = this->Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, previos);
					}

					ppp = this->Find_cell_point(P.coord[0] + P.Vel[0] * dt,
						P.coord[1] + P.Vel[1] * dt,
						P.coord[2] + P.Vel[2] * dt,
						0, previos);

					klk = 0;
					while (ppp == nullptr || ppp != P.cel)
					{
						klk++;
						if (klk > 100)
						{
							cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
							cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
							cout << previos->geo_parameters["l_size"] << " " << dt << endl;
							P.cel->Tecplot_print_cell();
							cout << "Error 3296411221" << endl;
							exit(-1);
						}

						if (klk == 99)
						{
							P.coord[0] = Center_cell[0];
							P.coord[1] = Center_cell[1];
							P.coord[2] = Center_cell[2];
						}
						else
						{
							Move[0] = (-P.coord[0] + Center_cell[0]) / 1000.0;
							Move[1] = (-P.coord[1] + Center_cell[1]) / 1000.0;
							Move[2] = (-P.coord[2] + Center_cell[2]) / 1000.0;
							P.Move(Move);
						}


						ppp = this->Find_cell_point(P.coord[0] + P.Vel[0] * dt,
							P.coord[1] + P.Vel[1] * dt,
							P.coord[2] + P.Vel[2] * dt,
							0, previos);
					}
				}

				Cell* previos = P.cel;
				Cell* ppp = this->Find_cell_point(P.coord[0], P.coord[1], P.coord[2], 0, previos);
				if (ppp != P.cel)
				{
					cout << "Error  8756121199" << endl;
					P.coord[0] = P.cel->center[0][0];
					P.coord[1] = P.cel->center[0][1];
					P.coord[2] = P.cel->center[0][2];
					whach(gr->number);
				}


				P.KSI = -log(1.0 - this->Sensors[sens_num1]->MakeRandom());
				P.I_do = 0.0;

				//cout << "FLY" << endl;


				this->MK_fly_immit(P, zone_MK, this->Sensors[sens_num2], Interpol, S_main); // Запускаем частицу в полёт   // !! Не написана
				//exit(-1);
				//cout << "END" << endl;
			}

			func->Delete();
			delete func;
			func = nullptr;
			
			
		}
	
		//exit(-1);                                         // DELETE

		// 8. Теперь надо сохранить ненужные выходящие функции распределения, но предварительно нормировать их
		// Ненужные функции распределения, это функции сорта nh_, при условии, что этот сорт не рождается в области
		// Сохраняем, удаляем, но перед этим нормируем
		cout << "Start delete_1  for sort " << nh_ + 1 << endl;
		if (this->MK_zone_H(zone_MK - 1, nh_) == false)
		{
			for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
			{
				auto& gr = this->MK_Grans[zone_MK - 1][idx];
				if (gr->type != Type_Gran::Us) continue;

				short int ni = 1; // Номер "выходящей" функции распределения
				if (gr->cells[0]->MK_zone == zone_MK)
				{
					ni = 0;
				}

				// Здесь надо для функции загрузить обратно ненужные массивы переменных
				//cout << "Re_Partially_free_space" << endl;
				gr->AMR[nh_][ni]->Re_Partially_free_space();

				//cout << "Normir_velocity_volume" << endl;
				gr->AMR[nh_][ni]->Normir_velocity_volume(gr->area[0]);
				if (this->phys_param->de_refine_AMR == true)
				{
					gr->AMR[nh_][ni]->de_Refine(nh_ + 1);
				}

				string name_f = "func_grans_AMR_" + to_string(ni) + "_H" +
					to_string(nh_ + 1) + "_" + to_string(gr->number) + ".bin";
				if (this->phys_param->save_AMR == true)
				{
					//cout << "Save" << endl;
					gr->AMR[nh_][ni]->Save(this->phys_param->AMR_folder + "/" + name_f);
				}
				//cout << "Delete" << endl;
				gr->AMR[nh_][ni]->Delete();
				delete gr->AMR[nh_][ni];
				gr->AMR[nh_][ni] = nullptr;
			}
		}
		cout << "End delete_1  for sort " << nh_ + 1 << endl;
	}


	cout << "**********************************" << endl;
	cout << "Obshee chislo chastic = " << ALL_N << endl;

	std::ofstream file1("info_AMR_size.txt", std::ios::app);
	for (short int iH = 0; iH < 9; iH++)
	{
		if (N2[iH] == 0) N2[iH] = 1;
		file1 << "Zone:  " << zone_MK << "   H = " << iH + 1 <<
			"   vixod size: " << 1.0 * N_vixod[iH] / N2[iH] << std::endl;
	}
	file1.close();



	// 9. Теперь нормируем оставшиеся функции распределения, сохраняем и удаляем
	if (true)
	{
		cout << "Start save_1" << endl;
		for (size_t j = 0; j < this->phys_param->num_H; j++)
		{
			if (this->MK_zone_H(zone_MK - 1, j) == true)
			{
				// Значит нам надо сохранить выходящую функцию для водорода сорта j
				for (size_t idx = 0; idx < this->MK_Grans[zone_MK - 1].size(); ++idx)
				{
					auto& gr = this->MK_Grans[zone_MK - 1][idx];
					if (gr->type != Type_Gran::Us) continue;

					short int ni = 1; // Номер "выходящей" функции распределения
					if (gr->cells[0]->MK_zone == zone_MK)
					{
						ni = 0;
					}

					gr->AMR[j][ni]->Re_Partially_free_space();

					gr->AMR[j][ni]->Normir_velocity_volume(gr->area[0]);
					if (this->phys_param->de_refine_AMR == true)
					{
						gr->AMR[j][ni]->de_Refine(j + 1);
					}
					string name_f = "func_grans_AMR_" + to_string(ni) + "_H" +
						to_string(j + 1) + "_" + to_string(gr->number) + ".bin";
					if (this->phys_param->save_AMR == true)
					{
						gr->AMR[j][ni]->Save(this->phys_param->AMR_folder + "/" + name_f);
					}
					gr->AMR[j][ni]->Delete();
					delete gr->AMR[j][ni];
					gr->AMR[j][ni] = nullptr;
				}
			}
		}

		cout << "End save_1" << endl;
	}



	// Нормировка Моментов в ячейках
	k1 = 0;
	cout << "Start: Normir moment in cells" << endl;
#pragma omp parallel for schedule(dynamic)
	for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
	{
		auto cell = this->All_Cell[idx];
		if (cell->MK_zone != zone_MK) continue;
		#pragma omp critical (first) 
		{
			k1++;
			if (k1 % 10000 == 0)
			{
				cout << "Cells = " << k1 << endl;
			}
		}

		cell->MK_normir_Moments(this->phys_param);

		if (this->phys_param->MK_source_S == true)
		{
			cell->MK_calc_Sm(this->phys_param);  // Нужно параллелить, так как эта функция долго обрабатывается
		}
	}
	cout << "End: Normir moment in cells" << endl;

	// Выведем одну функцию посмотреть что получилось)
	if (false)
	{
		for (auto& gr : this->MK_Grans[zone_MK - 1])
		{
			if (gr->type2 != Type_Gran_surf::BS) continue;

			if( fabs(fabs(gr->normal[0][0]) - 0.936503) > 0.0001) continue;
			if( fabs(fabs(gr->normal[0][1]) - 0.350126) > 0.0001) continue;
			if( fabs(fabs(gr->normal[0][2]) - 0.0193568) > 0.0001) continue;


			short int ni = 1; // Номер "выходящей" функции распределения
			if (gr->cells[0]->MK_zone == zone_MK)
			{
				ni = 0;
			}

			cout << gr->center[0][0] << " " << gr->center[0][1] << " " <<
				gr->center[0][2] << endl;
			cout << "Normal = " << gr->normal[0][0] << " " << gr->normal[0][1] << " " <<
				gr->normal[0][2] << endl;
			cout << "N_particle = " << gr->N_particle << endl;
			auto& func = gr->AMR[3][ni];
			cout << "n_func = " << func->parameters["n"] << endl;
			cout << "nn_func = " << func->parameters["nn"] << endl;
			cout << "Potok cherez funkcion = " << func->parameters["Smu"] << endl;
			double dF = func->parameters["nn"] / func->parameters["Smu"] * gr->area[0]
				- kv(func->parameters["n"] * gr->area[0] / func->parameters["Smu"]);
			cout << "dF = " << dF << endl;
			cout << "delta = " << sqrt(dF / func->parameters["Smu"]) / func->parameters["n"] << endl;
			func->Print_all_center_Tecplot(func, "test");
			func->Print_1D_Tecplot(func, -2.54327);
			func->Culk_SpotokV(gr->area[0]);
			cout << "Spotok = " << func->SpotokV / gr->area[0] << endl;
			break;
		}
	}


	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << "MK all time: " << duration.count() / 1000.0 / 60.0 << " minutes" << std::endl;
}

void Setka::MK_fly_immit(MK_particle& P, short int zone_MK, Sensor* Sens, Interpol* Interpol, Setka*& S_main)
{
	// S_main - это основная сетка (с большим числом ячеек)
	//          в ней хранятся массивы пикапов, их частот и т.д.


	/*cout << "______Start_MK_fly_immit___________" << endl;
	whach(P.coord[0]);
	whach(P.coord[1]);
	whach(P.coord[2]);
	cout << "_______________________________" << endl;*/
	//cout << "Start " << endl;
	//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;


	Eigen::Vector3d coord_init;
	Eigen::Vector3d Vel_init;

	coord_init << P.coord[0], P.coord[1], P.coord[2];
	Vel_init << P.Vel[0], P.Vel[1], P.Vel[2];


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

		//cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;                // DELETE
		if (P.cel->MK_zone != zone_MK)
		{
			cout << "Error 8767567487" << endl;
			cout << P.cel->MK_zone << "   " << zone_MK << endl;
			cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
			cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
			cout << k_cikl << endl;
			return;
			exit(-1);
		}

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
			if (P.cel->MK_zone != zone_MK)
			{
				cout << "Error 9871216655" << endl;
				cout << P.cel->MK_zone << "   " << zone_MK << endl;
				cout << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
				cout << P.Vel[0] << " " << P.Vel[1] << " " << P.Vel[2] << endl;
				cout << k_cikl << endl;
				exit(-1);
			}

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
					if(k1 < 10)
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
					perturbVectorKeepingMagnitude(P.Vel[0], P.Vel[1], P.Vel[2], 0.1 * norm2(P.Vel[0], P.Vel[1], P.Vel[2]));

					if (P.cel->MK_zone != zone_MK)
					{
						cout << "Error 9856312132" << endl;
					}


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

		// Получаем параметры плазмы в ячейке ----------------------------

		Cell_main = (*S_main).Find_cell_point(P.cel->center[0][0], P.cel->center[0][1], P.cel->center[0][2], 0, Cell_main_prev);
		if (Cell_main == nullptr)
		{
			Cell_main = (*S_main).Find_cell_point(P.cel->center[0][0] * 0.99, P.cel->center[0][1] * 0.99, P.cel->center[0][2] * 0.99, 0, Cell_main_prev);
			if (Cell_main == nullptr)
			{
				cout << "Error ihergiegufyiowehfvhuewygfiw" << endl;
				exit(-1);
			}
		}
		Cell_main_prev = Cell_main;

		short int zone = this->determ_zone(P.cel, 0);
		short int zone_main = this->determ_zone(Cell_main, 0);
		double ro, p, rho_He, cp, vx, vy, vz, rho_Th, p_Th;
		unordered_map<string, double> param;
		unordered_map<string, double> param2;
		double u, u1, u2, u3, skalar, nu_ex, sig = 0.0, nu_ex_pui_1, nu_ex_pui_2;
		double cp_sr = 0.0, u_sr = 0.0, u1_sr = 0.0, u2_sr = 0.0, u3_sr = 0.0, skalar_sr = 0.0;    // Средние параметры в ячейке при пролёте
		double vx_sr = 0.0, vy_sr = 0.0, vz_sr = 0.0;

		double I = P.I_do;
		double l = sqrt(kvv(time * P.Vel[0], time * P.Vel[1], time * P.Vel[2]));    // Расстояние, которое атом потенциально пролетает внутри ячейки
		double Vel_norm = sqrt(kvv(P.Vel[0], P.Vel[1], P.Vel[2]));                  // Модуль скорости атома
		
		
		//ro = P.cel->parameters[0]["rho"];
		//p = P.cel->parameters[0]["p"];
		//rho_He = P.cel->parameters[0]["rho_He"];
		//	// cp;// = sqrt(P.cel->parameters[0]["p"] / ro);
		//vx = P.cel->parameters[0]["Vx"];			// Скорости плазмы в ячейке
		//vy = P.cel->parameters[0]["Vy"];
		//vz = P.cel->parameters[0]["Vz"];

		//ro = Cell_main->parameters[0]["rho"];
		//p = Cell_main->parameters[0]["p"];
		rho_He = Cell_main->parameters[0]["rho_He"];
		// cp;// = sqrt(P.cel->parameters[0]["p"] / ro);
		vx = Cell_main->parameters[0]["Vx"];			// Скорости плазмы в ячейке
		vy = Cell_main->parameters[0]["Vy"];
		vz = Cell_main->parameters[0]["Vz"];

		vx_sr = vx;
		vy_sr = vy;
		vz_sr = vz;

		//Sootnosheniya(ro, p, rho_He, 0.0, 0.0, (int)(P.cel->type),
		//	rho_Th, rho_E, p_Th, p_Pui, T_Th, T_E);

		//this->phys_param->Plasma_components_1(zone, P.cel->parameters[0], param); // Это без пикапов
		S_main->phys_param->Plasma_components(zone_main, Cell_main->parameters[0], param, false);

		rho_Th = param["rho_Th"];
		p_Th = param["p_Th"];

		/*cout << "-------  " << zone << " | " << zone_main << "  | " << rho_Th << " |  " << p_Th << " |  " 
			<< rho_He << "  | " <<
			P.cel->parameters[0]["rho"] << "  | " << P.cel->parameters[0]["p"] << " |  " << 
			Cell_main->parameters[0]["rho"]
			<< "  | " << Cell_main->parameters[0]["p"] << endl;*/

		if (rho_Th <= 1e-8) rho_Th = 1e-8;
		if (p_Th <= 1e-8 / 2.0) p_Th = 1e-8 / 2.0;

		ro = rho_Th;
		cp = sqrt(2.0 * p_Th / rho_Th);

		// Постоянные поля для тестирования
		if (false)  
		{
			ro = 1.0;
			cp = 1.0;
			vx = this->phys_param->Velosity_inf;
			vy = 0.0;
			vz = 0.0;

			vx_sr = vx;
			vy_sr = vy;
			vz_sr = vz;
		}


		// ------------------------------
		// ------------------------------
		// ------------------------------
		// Находим частоты по перезарядке и другим процессам

		u = sqrt(kvv(P.Vel[0] - vx, P.Vel[1] - vy, P.Vel[2] - vz));
		u1 = vx - P.Vel[0];
		u2 = vy - P.Vel[1];
		u3 = vz - P.Vel[2];
		skalar = u1 * P.Vel[0] + u2 * P.Vel[1] + u3 * P.Vel[2];

		cp_sr = cp; 
		u_sr = u; 
		u1_sr = u1; 
		u2_sr = u2; 
		u3_sr = u3; 
		skalar_sr = skalar;

		if (u / cp > 7.0)
		{
			double uz = Velosity_1(u, cp);
			nu_ex = ro * uz * this->phys_param->sigma(uz) / this->phys_param->par_Kn;
		}
		else
		{
			nu_ex = (ro * this->phys_param->MK_int_1(u, cp)) / this->phys_param->par_Kn;  // Пробуем вычислять интеграллы численно
		}

		nu_ex_pui_1 = 0.0;
		nu_ex_pui_2 = 0.0;

		// Посчитаем частоты перезарядки на пикапах
		if (this->phys_param->is_PUI == true)
		{
			if (this->phys_param->pui_in_zone(zone - 1, 0) == true)
			{
				nu_ex_pui_1 = Cell_main->pui_get_nu(u, 0, this->phys_param->pui_wR) / this->phys_param->par_Kn;
				if (nu_ex_pui_1 < 0.0)
				{
					cout << "Error eijg9uerhguoiehg89pger  " << Cell_main->number << " " << nu_ex_pui_1 << " " <<
					u << endl;
					exit(-4);
				}
			}

			if (this->phys_param->pui_in_zone(zone - 1, 1) == true)
			{
				nu_ex_pui_2 = Cell_main->pui_get_nu(u, 1, this->phys_param->pui_wR) / this->phys_param->par_Kn;
				if (nu_ex_pui_2 < 0.0)
				{
					cout << "Error jtyu5rtygergergeg  " << Cell_main->number << " " << nu_ex_pui_2 << " " <<
						u << endl;
					exit(-4);
				}
			}
		}


		if (zone_main == 3)
		{
			if (nu_ex_pui_1 > nu_ex * 0.1 || nu_ex_pui_2 > nu_ex * 0.1)
			{
				cout << "Error uiehrgiehfouherf343   " << nu_ex << " " << nu_ex_pui_1 << " " << 
					nu_ex_pui_2 << endl;
				//exit(-1);
			}
		}

		if (std::isnan(cp_sr) || std::isnan(u1_sr))
		{
			std::cout << "ERROR frewrtgewr4 e4tewfwerfwf " << std::endl;
			cout << zone_main << " | " << rho_Th << " | " << p_Th << " | v =  " <<
				vx << " | " << vy << " | " << vz << " || " <<
				Cell_main->parameters[0]["rho"] << " | " << Cell_main->parameters[0]["p"]
				<< " | " << Cell_main->parameters[0]["rho_He"]
				<< " | MK_rho_Pui_1 = " << Cell_main->parameters[0]["MK_rho_Pui_1"]
				<< " | " << Cell_main->parameters[0]["MK_T_Pui_1"]
				<< " | " << Cell_main->parameters[0]["MK_rho_Pui_2"]
				<< " | " << Cell_main->parameters[0]["MK_T_Pui_2"] << endl;
			cout << cp_sr << " " << u1_sr << endl;
			exit(-1);
		}

		if (cp_sr > 10000000000.0)
		{
			std::cout << "ERROR egewgegegf egrerfegrvefg " << std::endl;
			cout << cp_sr << " " << u1_sr << endl;
			exit(-1);
		}

		double summ_nu = nu_ex + nu_ex_pui_1 + nu_ex_pui_2;

		if (summ_nu >= 0.000000001)
		{
			// Иначе если частота процессов нулевая, то в этой ячейке не произошло никакое событие
			sig = Vel_norm / summ_nu;
			I += l / sig;
		}


		if (P.cel->MK_zone != zone_MK)
		{
			cout << P.cel->number << " " << Cell_do->number << endl;
			cout << "Error 1654875068" << endl;
			exit(-1);
		}


		if (vtoroy_shans == false)
		{
			if (I < P.KSI)
			{
				P.I_do = I;  // В этом случае перезарядки в ячейке не произошло

				// Здесь записываем необходимые моменты в ячейку ---------------------
				if (this->phys_param->culc_cell_moments == true)
				{
					P.cel->MK_Add_particle(P, time, this->phys_param);
				}

				if (this->phys_param->culc_cell_source == true)
				{
					double kappa = nu_ex * time;
					double mu_ex = P.mu * (1.0 - exp(-kappa));   // Вес перезаряженного атома (фиктивная часть)
					P.cel->MK_Add_moment(P, cp_sr, u_sr, mu_ex, u1_sr, u2_sr, u3_sr, skalar_sr, this->phys_param);
				}


				short int zone = this->determ_zone(P.cel, 0);
				if (this->phys_param->MK_source_S == true)
				{
					P.cel->MK_Add_pui_source(P, u, nu_ex + nu_ex_pui_1 + nu_ex_pui_2, P.mu, time, this->phys_param, zone, 0);
				}
				// -------------------------------------------------------------------
			}
			else
			{
				double ksi = (P.KSI - P.I_do) * sig;
				double t_ex = ksi / Vel_norm;
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

				

				// Здесь записываем необходимые моменты в ячейку ---------------------
				if (this->phys_param->culc_cell_moments == true)
				{
					P.cel->MK_Add_particle(P, t_ex, this->phys_param);
				}

				if (this->phys_param->culc_cell_source == true)
				{
					double kappa = nu_ex * t_ex;
					double mu_ex = P.mu * (1.0 - exp(-kappa));   // Вес перезаряженного атома (фиктивная часть)
					P.cel->MK_Add_moment(P, cp_sr, u_sr, mu_ex, u1_sr, u2_sr, u3_sr, skalar_sr, this->phys_param);
				}

				if (this->phys_param->MK_source_S == true)
				{
					P.cel->MK_Add_pui_source(P, u_sr, nu_ex + nu_ex_pui_1 + nu_ex_pui_2, P.mu, t_ex, this->phys_param, zone, 0);
				}
				// -------------------------------------------------------------------
				// теперь нужно определить процесс, который произошёл
				double ksi_ = Sens->MakeRandom();

				// Матрицы взаимодействия сортов по областям
				Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>* hydrogen_arise;

				if (this->phys_param->is_PUI == true)
				{
					if (zone == 1)
					{
						hydrogen_arise = &this->phys_param->hydrogen_arise_1;
					}
					else if (zone == 2)
					{
						hydrogen_arise = &this->phys_param->hydrogen_arise_2;
					}
					else if (zone == 3)
					{
						hydrogen_arise = &this->phys_param->hydrogen_arise_3;
					}
					else if (zone == 4)
					{
						hydrogen_arise = &this->phys_param->hydrogen_arise_4;
					}
					else
					{
						cout << "Error geuihrbgeyurfge8754ty3tg" << endl;
						exit(-1);
					}
				}

				// В этом случае произошла перезарядка на тепловых протонах
				if (ksi_ <= nu_ex / summ_nu)
				{
					// Разыгрываем новую скорость
					double Ur, Uphi, Uthe;
					double Vr, Vphi, Vthe;
					double Wr, Wthe, Wphi;
					spherical_skorost(P.coord[0], P.coord[1], P.coord[2],
						vx_sr, vy_sr, vz_sr, Ur, Uphi, Uthe);
					spherical_skorost(P.coord[0], P.coord[1], P.coord[2],
						P.Vel[0], P.Vel[1], P.Vel[2], Vr, Vphi, Vthe);
					this->M_K_Change_Velosity(Sens, Ur / cp_sr, Uthe / cp_sr, Uphi / cp_sr,
						Vr / cp_sr, Vthe / cp_sr, Vphi / cp_sr, Wr, Wthe, Wphi, cp_sr);
					Wr *= cp_sr;
					Wthe *= cp_sr;
					Wphi *= cp_sr;

					dekard_skorost(P.coord[0], P.coord[1], P.coord[2],
						Wr, Wphi, Wthe, P.Vel[0], P.Vel[1], P.Vel[2]);

					P.sort = zone;
				}
				else if(ksi_ <=  (nu_ex + nu_ex_pui_1) / summ_nu) // Пикапы 1
				{
					if (this->phys_param->is_PUI == false)
					{
						cout << "Error rthryertgegrsy5ry45" << endl;
						exit(-1);
					}

					double uu, vv, ww;
					(*Cell_main).MK_pui_charge_exchange_velocity(Sens, S_main, this->phys_param,
						vx, vy, vz, P.Vel[0], P.Vel[1], P.Vel[2], uu, vv, ww, 0);
					P.Vel[0] = uu;
					P.Vel[1] = vv;
					P.Vel[2] = ww;

					if (hydrogen_arise->rows() < P.sort)
					{
						cout << "Error gergretg45ty45gerge" << endl;
						exit(-1);
					}

					int sss = (*hydrogen_arise)(P.sort - 1, 1);
					P.sort = sss;
					if (zone != 2 && sss == 6)
					{
						cout << "Error jegiurhguyeorf7893tf8er" << endl;
						exit(-1);
					}
				}
				else if (ksi_ <= (nu_ex + nu_ex_pui_1 + nu_ex_pui_2) / summ_nu) // Пикапы 1
				{
					if (this->phys_param->is_PUI == false)
					{
						cout << "Error htrhry45y45tyegergeg" << endl;
						exit(-1);
					}
					double uu, vv, ww;
					// Здесь не та ячейка!
					(*Cell_main).MK_pui_charge_exchange_velocity(Sens, S_main, this->phys_param,
						vx, vy, vz, P.Vel[0], P.Vel[1], P.Vel[2], uu, vv, ww, 1);
					P.Vel[0] = uu;
					P.Vel[1] = vv;
					P.Vel[2] = ww;

					if (hydrogen_arise->cols() < 3)
					{
						cout << "Error rewjhfgueyrghf983yt4r83" << endl;
						exit(-1);
					}

					if (hydrogen_arise->rows() < P.sort)
					{
						cout << "Error 56u5rhtruhrsy4r5y" << endl;
						exit(-1);
					}

					int sss = (*hydrogen_arise)(P.sort - 1, 2);
					P.sort = sss;
					if (zone != 2 && sss == 6)
					{
						cout << "Error rtyhy45t45tg54ye4gredgerg" << endl;
						exit(-1);
					}
				}
				else
				{
					cout << "Error e80u4g79eogh9p3h40983t" << endl;
					exit(-1);
				}

				
				P.KSI = -log(1.0 - Sens->MakeRandom());
				vtoroy_shans = false;
				if (P.cel->MK_zone != zone_MK)
				{
					cout << P.cel->number << " " << Cell_do->number << endl;
					cout << "Error 9088565453" << endl;
				}
				continue;
				//return this->MK_fly_immit(P, zone_MK, Sens);
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
			if (findSphereIntersectionTime(X, V, 1.01 * this->geo->R0,
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

		if (gran->Have_zone_number(zone_MK))
		{
			a1:
			// В этом случае долетели до границы, записываем что надо и выключаем частицу
			short int nn = 1;
			if (gran->cells[0]->MK_zone == zone_MK) nn = 0;
			auto AMR = gran->AMR[P.sort - 1][nn];

			if (this->phys_param->culc_AMR == true && gran->type == Type_Gran::Us)
			{
				// В грани на границе нет смысла ничего записывать
				AMR->Add_particle(P.Vel[0], P.Vel[1], P.Vel[2], P.mu); // мьютексы внутри
			}

			gran->mut.lock(); // Мьютекс для записи в гранб
			gran->N_particle++;
			gran->mut.unlock();

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
				if (gr->Have_zone_number(zone_MK))
				{
					gran = gr;
					goto a1;
				}
			}
		}

		//cout << "F " << endl;
		// В этом случае попали в следующую зону, пропустив граничную грань
		if (P.cel != nullptr && P.cel->MK_zone != zone_MK)
		{
			for (auto& gr : P.cel->grans)
			{
				if (gr->Have_zone_number(zone_MK))
				{
					gran = gr;
					goto a1;
				}
			}

			// В этом случае надо либо отключать ячейку (но мы потеряем часть массы)
			// либо запускать ей заново в этой ячейке (тогда наоборот получим лишнюю массу, так
			// как она уже записалась в данную ячейку

			if (vtoroy_shans == false)
			{
				P.coord[0] = coord_do[0] + (Cell_do->center[0][0] - coord_do[0]) / 300.0;
				P.coord[1] = coord_do[1] + (Cell_do->center[0][1] - coord_do[1]) / 300.0;
				P.coord[2] = coord_do[2] + (Cell_do->center[0][2] - coord_do[2]) / 300.0;
				P.cel = Cell_do;
				perturbVectorKeepingMagnitude(P.Vel[0], P.Vel[1], P.Vel[2], 0.01);
				vtoroy_shans = true;

				if (P.cel->MK_zone != zone_MK)
				{
					cout << "Error 1213563589" << endl;
				}

				continue;
			}

			cout << "Poteryal A" << endl;
			return;
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
					P.coord[i] += (gran->center[0][i] - P.coord[i])/200.0 + normal[i] * l_ / 200.0;
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
				perturbVector(P.coord[0], P.coord[1], P.coord[2], Cell_do->geo_parameters["l_size"]/100);
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
			return;
		}

		//cout << "G " << endl;
		//cout << "2 Soburausi otpravit " << endl;
		//cout << "P.coord = " << P.coord[0] << " " << P.coord[1] << " " << P.coord[2] << endl;
		vtoroy_shans = false;

		if (P.cel->MK_zone != zone_MK)
		{
			cout << "Error 4562190567" << endl;
		}

		continue;
		//return this->MK_fly_immit(P, zone_MK, Sens);
	}
}

double Setka::Get_Spotok_inf(const Eigen::Vector3d& n)
{
	Eigen::Vector3d V_inf_;
	V_inf_ << this->phys_param->Velosity_inf, 0.0, 0.0;

	double Ux = V_inf_.dot(n);

	return -(-exp(-kv(Ux)) + sqrtpi_ * Ux * erfc(Ux)) / (2.0 * sqrtpi_);
}

void Setka::Culc_f_pui_in_cell(Cell* Cel, Setka& S_MK, Interpol& SI_main, Interpol& SI_MK, bool Interpol_S)
{
	short int pui_nW = this->phys_param->pui_nW;
	double pui_wR = this->phys_param->pui_wR;
	double nH = this->phys_param->par_n_H_LISM;

	bool main_interpol = false;

	const short int dstep = 2;  // На сколько дробим шаг по времени (по сравнению со временем пролёта 1 ае)

	Eigen::VectorXd mas_w0(pui_nW);
	Eigen::VectorXd mas_w(pui_nW);
	Eigen::VectorXd f0_pui(pui_nW);
	Eigen::VectorXd f1_pui(pui_nW);
	Eigen::VectorXd mas_Sm(pui_nW);

	vector<double> mas_Sm_(pui_nW);
	vector<double> mas_Sp1_(pui_nW);
	vector<double> mas_Sp2_(pui_nW);

	Eigen::Vector3d r;
	r[0] = Cel->center[0][0];
	r[1] = Cel->center[0][1];
	r[2] = Cel->center[0][2];

	double dt = 0.001;
	double rho0 = Cel->parameters[0]["rho"];
	double qInt = 0.0;               // Интеграл от источника массы при ионизации
	double q1, rho, rho_do;
	Cell* prev = nullptr;
	Cell* A = nullptr;
	Cell* A_do = nullptr;
	short int zone_now;

	for (short int i = 0; i < pui_nW; i++)
	{
		mas_w0[i] = ((i + 0.5) * pui_wR / pui_nW);
		f0_pui[i] = 0.0;
		f1_pui[i] = 0.0;
		mas_Sm[i] = 0.0;
	}

	mas_w = mas_w0;

	short int zone = determ_zone(Cel, 0);
	unsigned int step = 0;

	std::array<Cell_handle, 6> prev_cell;
	std::array<Cell_handle, 6> next_cell;
	for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();

	Cell_handle prev_cell_ = Cell_handle();
	Cell_handle next_cell_ = nullptr;

	std::unordered_map<string, double> parameters;
	double x_do, y_do, z_do;

	//cout << "Start zone = " << zone << endl;

	if (zone == 1)
	{
		while (true)
		{
			step++; if (step > 1000000) { cout << "Infiniti cycle ERROR t8y5y65" << endl; exit(-1);}

			if (step > 990000)
			{
				cout << "step > 990000 " << r[0] << " " << r[1] << " " << r[2] << endl;
			}


			A = Find_cell_point(r[0], r[1], r[2], 0, prev);
			x_do = r[0];
			y_do = r[1];
			z_do = r[2];

			if (A == nullptr)
			{
				cout << "Error 9089h45hgtine5gg" << endl;
				exit(-1);
			}

			zone_now = determ_zone(A, 0);
			if (zone_now >= 2)
			{
				cout << "Error 8i67u5u4h4rhuytyjkt6" << endl;
				exit(-1);
			}

			//if (step % 20 == 0)  // Иногда корректируем шаг по времени
			//{
			//	dt = this->geo->R0 / norm2(A->parameters[0]["Vx"], A->parameters[0]["Vy"], A->parameters[0]["Vz"]) / dstep;
			//}

			q1 = 0.0;  // Если есть ионизация, надо сюда дописывать

			//cout << "C1 " << endl;
			bool bb = false;
			if(main_interpol) bb = SI_main.Get_param(r[0], r[1], r[2], parameters, prev_cell, next_cell);
			//cout << "C2 " << endl;

			if (bb == true)
			{
				for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i];

				rho = parameters["rho"];
				r[0] -= parameters["Vx"] * dt;
				r[1] -= parameters["Vy"] * dt;
				r[2] -= parameters["Vz"] * dt;
			}
			else
			{
				rho = A->parameters[0]["rho"];
				r[0] -= A->parameters[0]["Vx"] * dt;
				r[1] -= A->parameters[0]["Vy"] * dt;
				r[2] -= A->parameters[0]["Vz"] * dt;
			}

			qInt = qInt + dt * q1 / rho;

			//cout << "C3 " << endl;

			bool b;
			if (Interpol_S == true)
			{
				b = this->Get_pui_SS(mas_Sm_, mas_Sp1_, mas_Sp2_, 1, x_do, y_do, z_do,
					S_MK, SI_MK, prev_cell_, next_cell_);
			}

			for (short int iw = 0; iw < pui_nW; iw++)
			{
				short int numw = min(int(mas_w[iw] / pui_wR * pui_nW), pui_nW - 1);
				if (mas_w(iw) < pui_wR && mas_w(iw) > 0)
				{
					double Sm, Sp;
					// Сначала интерполируем Sm и Sp
					if (Interpol_S == true)
					{
						Sm = mas_Sm_[numw];
						if (b == false)
						{
							//cout << "ERROR wergvwevrtgewte44" << endl;
							//Sm = A->pui_Sm[numw];
							Sm = 0.0;
						}

						Sp = mas_Sp1_[numw];
						if (b == false)
						{
							//cout << "ERROR werfg345tb4t345" << endl;
							//Sp = A->pui_Sp(0, numw);
							Sp = 0.0;
						}
					}
					else
					{
						Sm = A->pui_Sm[numw];
						Sp = A->pui_Sp(0, numw);
					}

					mas_Sm(iw) = mas_Sm(iw) + nH * Sm * dt;
					f0_pui(iw) = f0_pui(iw) + nH * Sp * dt * exp(-mas_Sm(iw));  // Это S + , просто сразу накапливаем в функцию распределения
				}
				mas_w(iw) = mas_w0(iw) / (pow((rho0 / rho), (1.0 / 3.0)) * exp(-1.0 / 3.0 * qInt));
			}
			//cout << "C4 " << endl;

			if (r.norm() < 1.01 * this->geo->R0) break;
		}

		for (short int iw = 0; iw < pui_nW; iw++)
		{
			Cel->f_pui_1[iw] = f0_pui[iw];
		}
	}
	else if (zone == 2)
	{
		unsigned int step_in_cell = 0;
		while (true)
		{
			if(A!= nullptr) A_do = A;
			step++; if (step > 1000000) { cout << "Infiniti cycle ERROR ertert34634rt34tgewrg" << endl; 
			cout << A->parameters[0]["Vx"] << " " << A->parameters[0]["Vy"] << " " << A->parameters[0]["Vz"] << endl; 
			exit(-1); }

			if (step > 990000)
			{
				cout << "AB step > 990000 " << r[0] << " " << r[1] << " " << r[2] << "   " << zone << " " << zone_now << endl;
			}


			A = Find_cell_point(r[0], r[1], r[2], 0, prev);
			if (A == nullptr)
			{
				r = r / 1.0001;
				continue;
				//cout << "Error hrtgth45t4twtiet5" << endl;
				//cout << r[0] << " " << r[1] << " " << r[2] << endl;
				//exit(-1);
			}
			x_do = r[0];
			y_do = r[1];
			z_do = r[2];

			if (A == A_do)
			{
				step_in_cell++;
			}
			else
			{
				step_in_cell = 0;
			}

			zone_now = determ_zone(A, 0);
			if (zone_now >= 3 || step_in_cell > 1000)
			{
				if (r[0] > 0 || sqrt(kv(r[1]) + kv(r[2])) < 20.0)
				{
					r = r / 1.03;
				}
				else
				{
					r[1] = r[1] / 1.03;
					r[2] = r[2] / 1.03;
				}
				continue;
			}

			//if (step % 20 == 0)  // Иногда корректируем шаг по времени
			//{
			//	dt = this->geo->R0 / norm2(A->parameters[0]["Vx"], A->parameters[0]["Vy"], A->parameters[0]["Vz"]) / dstep;
			//}

			//dt = this->geo->R0 / norm2(A->parameters[0]["Vx"], A->parameters[0]["Vy"], A->parameters[0]["Vz"]) / dstep;

			if (zone_now == 1) break;

			bool bb = false;
			if (main_interpol) bb = SI_main.Get_param(r[0], r[1], r[2], parameters, prev_cell, next_cell);


			if (bb == true)
			{
				for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i];

				rho = parameters["rho"];
				r[0] -= parameters["Vx"] * dt;
				r[1] -= parameters["Vy"] * dt;
				r[2] -= parameters["Vz"] * dt;
			}
			else
			{
				rho = A->parameters[0]["rho"];
				r[0] -= A->parameters[0]["Vx"] * dt;
				r[1] -= A->parameters[0]["Vy"] * dt;
				r[2] -= A->parameters[0]["Vz"] * dt;
			}

			q1 = 0.0;  // Если есть ионизация, надо сюда дописывать
			rho_do = rho;
			qInt = qInt + dt * q1 / rho;

			bool b;
			if (Interpol_S == true)
			{

				b = this->Get_pui_SS(mas_Sm_, mas_Sp1_, mas_Sp2_, 2, x_do, y_do, z_do,
					S_MK, SI_MK, prev_cell_, next_cell_);
			}

			for (short int iw = 0; iw < pui_nW; iw++)
			{
				short int numw = min(int(mas_w[iw] / pui_wR * pui_nW), pui_nW - 1);
				if (mas_w(iw) < pui_wR && mas_w(iw) > 0)
				{
					double Sm, Sp1, Sp2;
					// Сначала интерполируем Sm и Sp
					if (Interpol_S == true)
					{
						Sm = mas_Sm_[numw];
						if (b == false)
						{
							//cout << "ERROR fwerfe4tf34tf345t34" << endl;
							//Sm = A->pui_Sm[numw];
							Sm = 0.0;
						}

						Sp1 = mas_Sp1_[numw];
						if (b == false)
						{
							//cout << "ERROR wergvwevrtg34tr34tr34" << endl;
							Sp1 = 0.0;
							//Sp1 = A->pui_Sp(0, numw);
						}

						Sp2 = mas_Sp2_[numw];
						if (b == false)
						{
							//cout << "ERROR wergvwevertertgewte44" << endl;
							Sp2 = 0.0;
							//Sp2 = A->pui_Sp(1, numw);
						}
					}
					else
					{
						Sm = A->pui_Sm[numw];
						Sp1 = A->pui_Sp(0, numw);
						Sp2 = A->pui_Sp(1, numw);
					}

					mas_Sm(iw) = mas_Sm(iw) + nH * Sm * dt;
					f0_pui(iw) = f0_pui(iw) + nH * Sp1 * dt * exp(-mas_Sm(iw));  // Это S + , просто сразу накапливаем в функцию распределения
					f1_pui(iw) = f1_pui(iw) + nH * Sp2 * dt * exp(-mas_Sm(iw));  // Это S + , просто сразу накапливаем в функцию распределения
				}
				mas_w(iw) = mas_w0(iw) / (pow((rho0 / rho), (1.0 / 3.0)) * exp(-1.0 / 3.0 * qInt));
			}
		}

		// Ищем грань-TS в текущей ячейке
		Gran* G = nullptr;
		for (auto& gr : A->grans)
		{
			if (gr->type2 == Type_Gran_surf::TS)
			{
				G = gr;
				break;
			}
		}

		if (G == nullptr) 
		{ 
			for (auto& gr : A_do->grans)
			{
				if (gr->type2 == Type_Gran_surf::TS)
				{
					G = gr;
					break;
				}
			}

			if (G == nullptr)
			{
				cout << "ERRROR  3i4j34hfg34tf3" << endl;
				cout << r[0] << " " << r[1] << " " << r[2] << endl;
				cout << dt << " " << A->parameters[0]["Vx"] << " " << A->parameters[0]["Vy"] << 
					" " << A->parameters[0]["Vz"] << endl;
				exit(-2);
			}
		}

		double s = rho_do / rho;
		Eigen::Vector3d B;
		Eigen::Vector3d normal;
		B << A->parameters[0]["Bx"], A->parameters[0]["By"], A->parameters[0]["Bz"];
		normal << G->normal[0][0], G->normal[0][1], G->normal[0][2];
		double cospsi = B.dot(normal) / (B.norm());
		double AA = sqrt(kv(cospsi) + kv(s) * (1.0 - kv(cospsi)));
		double BB = kv(s) / kv(AA);
		double C = (2.0 * AA + BB) / 3.0;

		rho0 = rho;
		qInt = 0.0;
		mas_w0 = mas_w / sqrt(C);
		mas_w = mas_w0;

		// Теперь бежим до Солнца
		step = 0;
		while (true)
		{
			step++; if (step > 1000000) { cout << "Infiniti cycle ERROR gegeyteg45hytik979687" << endl; exit(-1); }
			A = Find_cell_point(r[0], r[1], r[2], 0, prev);
			if (A == nullptr)
			{
				cout << "Error 9089h45hgtine5gg" << endl;
				exit(-1);
			}

			zone_now = determ_zone(A, 0);
			if (zone_now >= 2)
			{
				cout << "Error 8i67u5u4h4rhuytyjkt6" << endl;
				exit(-1);
			}

			if (step % 20 == 0)  // Иногда корректируем шаг по времени
			{
				//dt = this->geo->R0 / norm2(A->parameters[0]["Vx"], A->parameters[0]["Vy"], A->parameters[0]["Vz"]) / dstep;
			}

			x_do = r[0];
			y_do = r[1];
			z_do = r[2];

			bool bb = false;
			if (main_interpol) bb = SI_main.Get_param(r[0], r[1], r[2], parameters, prev_cell, next_cell);

			if (bb == true)
			{
				for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i];

				rho = parameters["rho"];
				r[0] -= parameters["Vx"] * dt;
				r[1] -= parameters["Vy"] * dt;
				r[2] -= parameters["Vz"] * dt;
			}
			else
			{
				rho = A->parameters[0]["rho"];
				r[0] -= A->parameters[0]["Vx"] * dt;
				r[1] -= A->parameters[0]["Vy"] * dt;
				r[2] -= A->parameters[0]["Vz"] * dt;
			}


			q1 = 0.0;  // Если есть ионизация, надо сюда дописывать
			qInt = qInt + dt * q1 / rho;

			bool b;
			if (Interpol_S == true)
			{
				b = this->Get_pui_SS(mas_Sm_, mas_Sp1_, mas_Sp2_, 1, x_do, y_do, z_do,
					S_MK, SI_MK, prev_cell_, next_cell_);
			}

			for (short int iw = 0; iw < pui_nW; iw++)
			{
				short int numw = min(int(mas_w[iw] / pui_wR * pui_nW), pui_nW - 1);
				if (mas_w(iw) < pui_wR && mas_w(iw) > 0)
				{
					double Sm, Sp1;
					// Сначала интерполируем Sm и Sp
					if (Interpol_S == true)
					{
						double Sm = mas_Sm_[numw];
						if (b == false)
						{
							//cout << "ERROR wergvwevrtgewtertferfefewee44" << endl;
							Sm = 0.0;
						}

						double Sp1 = mas_Sp1_[numw];
						if (b == false)
						{
							//cout << "ERROR wergvwevr3453452345tgewte44" << endl;
							Sp1 = 0.0;
						}
					}
					else
					{
						Sm = A->pui_Sm[numw];
						Sp1 = A->pui_Sp(0, numw);
					}


					mas_Sm(iw) = mas_Sm(iw) + nH * Sm * dt;
					f0_pui(iw) = f0_pui(iw) + nH * Sp1 * dt * exp(-mas_Sm(iw)) * s / pow(C, 1.5);  // Это S + , просто сразу накапливаем в функцию распределения
				}
				mas_w(iw) = mas_w0(iw) / (pow((rho0 / rho), (1.0 / 3.0)) * exp(-1.0 / 3.0 * qInt));
			}

			if (r.norm() < 1.01 * this->geo->R0) break;
		}

		for (short int iw = 0; iw < pui_nW; iw++)
		{
			Cel->f_pui_1[iw] = f0_pui[iw];
			Cel->f_pui_2[iw] = f1_pui[iw];
		}
	}
	else if (zone >= 3)
	{
		dt = 0.001;
		while (true)
		{
			step++; if (step > 1000000) { cout << "Infiniti cycle ERROR  geget34" << endl; exit(-1); }
			if (step > 990000)
			{
				cout << "AC step > 990000 " << r[0] << " " << r[1] << " " << r[2] << "   " << zone << " " << zone_now << endl;
			}

			A = Find_cell_point(r[0], r[1], r[2], 0, prev);
			if (A == nullptr) break;

			zone_now = determ_zone(A, 0);
			if (zone_now == 1)
			{
				cout << "Error jytk7896yh4y4h4t345tg" << endl;
				exit(-1);
			}

			if (zone_now == 2)
			{
				if (r[0] > 0)
				{
					r = r * 1.0005;
				}
				else
				{
					r[1] = r[1] * 1.0005;
					r[2] = r[2] * 1.0005;
				}
				continue;
			}

			if (step % 20 == 0)  // Иногда корректируем шаг по времени
			{
				//dt = this->geo->R0 / norm2(A->parameters[0]["Vx"], A->parameters[0]["Vy"], A->parameters[0]["Vz"]) / dstep;
			}

			x_do = r[0];
			y_do = r[1];
			z_do = r[2];

			bool bb = false;
			if (main_interpol) bb = SI_main.Get_param(r[0], r[1], r[2], parameters, prev_cell, next_cell);


			if (bb == true)
			{
				for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i];

				rho = parameters["rho"];
				r[0] -= parameters["Vx"] * dt;
				r[1] -= parameters["Vy"] * dt;
				r[2] -= parameters["Vz"] * dt;
			}
			else
			{
				rho = A->parameters[0]["rho"];
				r[0] -= A->parameters[0]["Vx"] * dt;
				r[1] -= A->parameters[0]["Vy"] * dt;
				r[2] -= A->parameters[0]["Vz"] * dt;
			}

			q1 = 0.0;  // Если есть ионизация, надо сюда дописывать
			qInt = qInt + dt * q1 / rho;

			bool b;
			if (Interpol_S == true)
			{
				b = this->Get_pui_SS(mas_Sm_, mas_Sp1_, mas_Sp2_, 1, x_do, y_do, z_do,
					S_MK, SI_MK, prev_cell_, next_cell_);
			}

			for (short int iw = 0; iw < pui_nW; iw++)
			{
				short int numw = min(int(mas_w[iw] / pui_wR * pui_nW), pui_nW - 1);
				if (mas_w(iw) < pui_wR && mas_w(iw) > 0)
				{
					double Sm, Sp1;
					// Сначала интерполируем Sm и Sp
					if (Interpol_S == true)
					{
						double Sm = mas_Sm_[numw];
						if (b == false)
						{
							//cout << "ERROR wergvwevrtgewterwrwere44" << endl;
							//cout << x_do << " " <<  y_do << " " << z_do << endl;
							//Sm = A->pui_Sm[numw];
							Sm = 0.0;
						}

						double Sp1 = mas_Sp1_[numw];
						if (b == false)
						{
							//cout << "ERROR wergvwevrtgewtewerwqerwer232244" << endl;
							//Sp1 = A->pui_Sp(0, numw);
							Sp1 = 0.0;
						}
					}
					else
					{
						Sm = A->pui_Sm[numw];
						Sp1 = A->pui_Sp(0, numw);
					}

					mas_Sm(iw) = mas_Sm(iw) + nH * Sm * dt;
					f0_pui(iw) = f0_pui(iw) + nH * Sp1 * dt * exp(-mas_Sm(iw));  // Это S + , просто сразу накапливаем в функцию распределения
				}
				mas_w(iw) = mas_w0(iw) / (pow((rho0 / rho), (1.0 / 3.0)) * exp(-1.0 / 3.0 * qInt));
			}
		}

		for (short int iw = 0; iw < pui_nW; iw++)
		{
			Cel->f_pui_1[iw] = f0_pui[iw];
		}
	}

}

void Setka::mas_pogl_Culc(const double& ex, const double& ey, const double& ez, const string& name)
{
	ofstream fout;
	string name_f = "poglosh_" + name + ".txt";
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = u, ";

	for (int i = 0; i < this->phys_param->num_H; i++)
	{
		fout << "f_" + to_string(i + 1) << ", f_moment_" + to_string(i + 1) << ", ";
	}

	fout << "fAll, fALL_moment" << endl;

	Eigen::Vector3d e;
	Eigen::Vector3d r;
	Cell* A, *prev;
	prev = nullptr;
	// Для трёх способов расчёта поглощения создаём массивы
	Eigen::MatrixXd mas_pogl;      // (sort, n)    Массив поглощения
	Eigen::MatrixXd mas_pogl2;      // (sort, n)    Массив поглощения для флюидов (я его убрал - его надо считать на основной сетке)
	Eigen::MatrixXd mas_pogl3;      // (sort, n)    Массив поглощения для моментов водорода
	mas_pogl.resize(this->phys_param->num_H, this->phys_param->pogl_n);
	mas_pogl.setZero();
	mas_pogl2.resize(this->phys_param->num_H, this->phys_param->pogl_n);
	mas_pogl2.setZero();
	mas_pogl3.resize(this->phys_param->num_H, this->phys_param->pogl_n);
	mas_pogl3.setZero();

	e << ex, ey, ez;
	double ee = e.norm();
	e /= ee;

	r = e * phys_param->R_0 * 1.1;  // 1.1
	//r = e * 14.0;  // 1.1
	double dr = phys_param->R_0 / 5.0;
	double dv = (this->phys_param->pogl_R - this->phys_param->pogl_L) / this->phys_param->pogl_n;
	double u1, u2, u3, c, n, p;
	double u1_MK, u2_MK, u3_MK, c_MK, n_MK, p_MK;

	double My_S = 0.0;

	while (true)
	{
		r += e * dr;
		A = this->Find_cell_point(r[0], r[1], r[2], 0, prev);

		if (A == nullptr) break;


		for (int i = 0; i < this->phys_param->num_H; i++)
		{
			/*u1 = A->parameters[0]["Vx_H" + to_string(i + 1)];
			u2 = A->parameters[0]["Vy_H" + to_string(i + 1)];
			u3 = A->parameters[0]["Vz_H" + to_string(i + 1)];
			n = A->parameters[0]["rho_H" + to_string(i + 1)];
			p = A->parameters[0]["p_H" + to_string(i + 1)];
			c = sqrt(2.0 * p / n);*/

			u1_MK = A->parameters[0]["MK_Vx_H" + to_string(i + 1)];
			u2_MK = A->parameters[0]["MK_Vy_H" + to_string(i + 1)];
			u3_MK = A->parameters[0]["MK_Vz_H" + to_string(i + 1)];
			n_MK = A->parameters[0]["MK_n_H" + to_string(i + 1)];
			p_MK = A->parameters[0]["MK_T_H" + to_string(i + 1)];
			c_MK = sqrt(p_MK);

			for (int j = 0; j < this->phys_param->pogl_n; j++)
			{
				double v = this->phys_param->pogl_L + dv * (j + 0.5);

				mas_pogl(i, j) += A->mas_pogl(i, j);
				My_S += A->mas_pogl(i, j);
				
				//mas_pogl2(i, j) +=  n *
				//	exp(-(kv(v - u1 * e[0] - u2 * e[1] - u3 * e[2])) / kv(c)) / (sqrt_pi * c);

				mas_pogl3(i, j) += n_MK *
					exp(-(kv(v - u1_MK * e[0] - u2_MK * e[1] - u3_MK * e[2])) / kv(c_MK)) / (sqrt_pi * c_MK);
			}
		}

		if (r[0] > this->phys_param->R_MK_Max) break;
		//if (r[0] > 14.0) break;
	}

	cout << "My_S = " << My_S << endl;
	cout << this->phys_param->par_n_H_LISM << " " << this->phys_param->par_poglosh << " " << dr / dv << endl;

	for (int j = 0; j < this->phys_param->pogl_n; j++)
	{
		fout << this->phys_param->pogl_L + dv * (j + 0.5) << " ";
		double S = 0.0;
		double SS = 0.0;
		double SSS = 0.0;
		for (int i = 0; i < this->phys_param->num_H; i++)
		{
			mas_pogl(i, j) *= this->phys_param->par_n_H_LISM * this->phys_param->par_poglosh * dr / dv;
			//mas_pogl2(i, j) *= this->phys_param->par_n_H_LISM * this->phys_param->par_poglosh * dr;
			mas_pogl3(i, j) *= this->phys_param->par_n_H_LISM * this->phys_param->par_poglosh * dr;
			fout << exp(-mas_pogl(i, j)) << " " << exp(-mas_pogl3(i, j)) << " ";
			S += mas_pogl(i, j);
			//SS += mas_pogl2(i, j);
			SSS += mas_pogl3(i, j);
		}
		fout << exp(-S) << " " << exp(-SSS) << " " << endl;
	}
	
	fout.close();
}


void Setka::Velosity_initial(Sensor* s, Eigen::Vector3d& V,
	const Eigen::Vector3d& n, const Eigen::Vector3d& t, 
	const Eigen::Vector3d& m)
{
	Eigen::Vector3d V_inf_;
	Eigen::Vector3d V_inf;
	V_inf_ << this->phys_param->Velosity_inf, 0.0, 0.0;

	V_inf[0] = V_inf_.dot(n);
	V_inf[1] = V_inf_.dot(t);
	V_inf[2] = V_inf_.dot(m);

	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = fabs(V_inf[0]) * sqrtpi_ /
		(1.0 + fabs(V_inf[0]) * sqrtpi_);

	do
	{
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();

		if (p1 > ksi3)
		{
			z = cos(const_pi * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			if (ksi4 <= 0.5)
			{
				z = -sqrt(-log(2.0 * ksi4));
			}
			else
			{
				z = sqrt(-log(2.0 * (1.0 - ksi4)));
			}
		}
	} while (fabs(z + V_inf[0]) / (fabs(V_inf[0]) + fabs(z)) < ksi6 || z > -V_inf[0]);

	double V1 = z + V_inf[0];

	double ksi1 = s->MakeRandom();
	double ksi2 = s->MakeRandom();
	double a = sqrt(-log(ksi2));
	double V2 = V_inf[1] + a * cos(2.0 * const_pi * ksi1);
	double V3 = V_inf[2] + a * sin(2.0 * const_pi * ksi1);

	V = V1 * n + V2 * t + V3 * m;

	return;
}

bool Setka::Time_to_vilet(const MK_particle& P, double& time, Gran*& gran)
{
	Cell* C = P.cel;
	Eigen::Vector3d R, V;

	R[0] = P.coord[0];
	R[1] = P.coord[1];
	R[2] = P.coord[2];
	V[0] = P.Vel[0];
	V[1] = P.Vel[1];
	V[2] = P.Vel[2];

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


void Setka::M_K_Change_Velosity(Sensor* sens, const double& Ur, const double& Uthe,
	const double& Uphi, const double& Vr, const double& Vthe, 
	const double& Vphi, double& Wr, double& Wthe, double& Wphi, const double& cp)
{
	double X = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double p4 = 0.5 * sqrtpi_ * X / (1.0 + 0.5 * sqrtpi_ * X);
	double om1, om2, om3, lo;
	double y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h;

	double gg = 0.0;
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6;

	do
	{
		ksi1 = sens->MakeRandom();
		ksi2 = sens->MakeRandom();
		ksi3 = sens->MakeRandom();
		ksi4 = sens->MakeRandom();
		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * const_pi * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * const_pi * ksi5);
			// Более экономичный алгоритм   --  выйгрыша нет вроде от него
			/*do
			{
				om2 = 1.0 - 2.0 * sens->MakeRandom();
				om3 = 1.0 - 2.0 * sens->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1.0);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;*/

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(const_pi * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * const_pi * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * const_pi * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uuu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
		h = ((uuu * this->phys_param->sigma2(uuu, cp)) / (this->phys_param->sigma2(X, cp) * (X + yy)));
	} while (h < ksi6); 


	Wr = v1;
	Wthe = v2;
	Wphi = v3;


	return;
}