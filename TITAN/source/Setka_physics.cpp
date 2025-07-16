#include "Setka.h"

#include <omp.h>


void Setka::Set_Gran_par_for_interpolate(void)
{
	cout << "Start: Set_Gran_par_for_interpolate" << endl;

	auto names = this->phys_param->param_names;

	// Заполняем параметры на грани
	for (auto& gr : this->All_Gran)
	{
		if (gr->cells.size() == 1)
		{
			for (auto& nam : names)
			{
				gr->parameters[nam] = gr->cells[0]->parameters[0][nam];
			}
		}
		else
		{
			Eigen::Vector3d vec1, vec2;
			vec1 << gr->center[0][0] - gr->cells[0]->center[0][0],
				gr->center[0][1] - gr->cells[0]->center[0][1],
				gr->center[0][2] - gr->cells[0]->center[0][2];
			vec2 << gr->center[0][0] - gr->cells[1]->center[0][0],
				gr->center[0][1] - gr->cells[1]->center[0][1],
				gr->center[0][2] - gr->cells[1]->center[0][2];

			for (auto& nam : names)
			{
				gr->parameters[nam] = (gr->cells[0]->parameters[0][nam] * vec2.norm() +
					gr->cells[1]->parameters[0][nam] * vec1.norm()) / (vec1.norm() + vec2.norm());
			}
		}
	}

	cout << "A" << endl;
	// Заполняем параметры в узлах
	for (auto& yz : this->All_Yzel)
	{
		double xc, yc, zc;
		xc = yz->coord[0][0];
		yc = yz->coord[0][1];
		zc = yz->coord[0][2];


		if (yz->grans.size() == 0)
		{
			cout << "Error  9764397648" << endl;
			continue;
		}

		for (auto& nam : names)
		{

			double sum_weights = 0.0;
			double sum_weighted_values = 0.0;

			for (auto& gr : yz->grans)
			{
				if (gr->parameters.find(nam) == gr->parameters.end()) continue;

				double dist = norm2(gr->center[0][0] - xc,
					gr->center[0][1] - yc, gr->center[0][2] - zc);
				double weight = 1.0 / std::pow(dist, 1.0);
				sum_weights += weight;
				sum_weighted_values += weight * gr->parameters[nam];
			}


			yz->parameters[nam] = sum_weighted_values / sum_weights;
		}
	}

	cout << "B" << endl;
	// Находим интерполяцию в ячейке
	for (auto& cell : this->All_Cell)
	{
		std::vector<Eigen::Vector3d> points;
		std::vector<double> values;

		for (auto& p : cell->yzels)
		{
			Eigen::Vector3d ppp;
			ppp << p->coord[0][0], p->coord[0][1],
				p->coord[0][2];
			points.push_back(ppp);
			values.push_back(p->parameters["rho"]);
		}

		size_t n = points.size();
		Eigen::MatrixXd A(n, n);

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				A(i, j) = rbfKernel((points[i] - points[j]).norm());
			}
		}

		cell->interpolate_alpha["rho"] = A.colPivHouseholderQr().solve(
			Eigen::Map<const Eigen::VectorXd>(values.data(), n));

	}


	cout << "End: Set_Gran_par_for_interpolate" << endl;
}


void Setka::Init_boundary_grans(void)
{
	this->Calculating_measure(0);

	int k = 0;
	// Выделяем внутренние узлы
	for (auto& i : All_Luch)
	{
		if (i->type == "A_Luch" || i->type == "A2_Luch" || i->type == "B_Luch" ||
			i->type == "C_Luch" || i->type == "C2_Luch")
		{
			for (size_t j = 0; j < this->geo->M0 + 1 + 7; j++) // this->geo->M0 + 1
			{
				i->Yzels[j]->is_inner = true;
				k++;
			}
		}
	}
	whach(k);

	// Разделяем внутренние и внешние ячейки
	for (auto& i : this->All_Cell)
	{
		bool aa = true;  // внутренняя ли ячейка
		for (auto& j : i->yzels)
		{
			if (j->is_inner == false)
			{
				aa = false;
				break;
			}
		}

		if (aa == true)
		{
			i->is_inner = true;
			this->Cell_inner_area.push_back(i);
		}
		else
		{
			i->is_inner = false;
			this->Cell_outer_area.push_back(i);
		}
	}

	whach(this->Cell_inner_area.size());
	whach(this->Cell_outer_area.size());
	whach(this->All_Cell.size());

	if (this->Cell_inner_area.size() + this->Cell_outer_area.size() != this->All_Cell.size())
	{
		cout << "Error 8564532456" << endl;
		exit(-1);
	}

	// Заполняем внутренние и внешние грани
	for (auto& i : this->All_Gran)
	{
		if (i->cells.size() == 1)
		{
			if (i->cells[0]->is_inner == true)
			{
				this->Gran_inner_area.push_back(i);
			}
			else
			{
				this->Gran_outer_area.push_back(i);
			}
		}
		else
		{
			if (i->cells[0]->is_inner == true || i->cells[1]->is_inner == true)
			{
				this->Gran_inner_area.push_back(i);
			}

			if (i->cells[0]->is_inner == false || i->cells[1]->is_inner == false)
			{
				this->Gran_outer_area.push_back(i);
			}
		}
	}

	whach(this->Gran_inner_area.size());
	whach(this->Gran_outer_area.size());
	whach(this->All_Gran.size());

	if (this->Gran_inner_area.size() > this->All_Gran.size() ||
		this->Gran_outer_area.size() > this->All_Gran.size() ||
		this->Gran_inner_area.size() + this->Gran_outer_area.size() <= this->All_Gran.size())
	{
		cout << "Error 6567876120" << endl;
		exit(-1);
	}

	// Находим граничные грани
	size_t oh = 0;
	size_t ih = 0;
	size_t os = 0;
	for (auto& i : this->All_Gran)
	{
		i->type = Type_Gran::Us;

		if (i->cells.size() != 1) continue;

		if (i->func_R(0) <= this->geo->R0 + 0.00001)
		{
			i->type = Type_Gran::Inner_Hard;
			ih++;
			this->All_boundary_Gran.push_back(i);
		}
		//else if(i->center[0][0] > 0.0)
		else if(i->center[0][0] > this->geo->L7 + 0.00001)
		{
			i->type = Type_Gran::Outer_Hard;
			oh++;
			this->All_boundary_Gran.push_back(i);
		}
		else
		{
			i->type = Type_Gran::Outer_Soft;
			os++;
			this->All_boundary_Gran.push_back(i);
		}
	}

	cout << "Granichniy grans: " << endl;
	whach(oh);
	whach(ih);
	whach(os);
}

void Setka::Init_physics(void)
{
	this->Calculating_measure(0);
	this->Calculating_measure(1);
	double x, y, z, r, the;
	Eigen::Vector3d vec, cc, vv;
	double BR, BPHI, V1, V2, V3, mV;

	// Редактирование каких-то переменных
	if (false)
	{
		for (auto& i : this->All_Cell)
		{
			if (i->center[0][0] > 200)
			{
				i->parameters[0]["rho"] = 1.60063; // 1.0 * (1.0 + this->phys_param->mn_He_inf);
				i->parameters[0]["n_He"] = 0.6; // this->phys_param->mn_He_inf;
				i->parameters[0]["p"] = 1.15; // 1 + (i->parameters["n_He"]) /
					//(i->parameters["rho"] - i->parameters["n_He"]);
				i->parameters[0]["Vx"] = this->phys_param->Velosity_inf;
				i->parameters[0]["Vy"] = 0.0;
				i->parameters[0]["Vz"] = 0.0;
				i->parameters[0]["Bx"] = -this->phys_param->B_inf * cos(this->phys_param->alphaB_inf);
				i->parameters[0]["By"] = -this->phys_param->B_inf * sin(this->phys_param->alphaB_inf);
				i->parameters[0]["Bz"] = 0.0;
				i->parameters[0]["Q"] = 100.0 * i->parameters[0]["rho"];
			}

		}
	}

	// Если ввели какие-то новые переменные, их надо заполнить
	if (false)
	{
		for (auto& i : this->All_Cell)
		{
			i->parameters[0]["n_He"] = 0.000001;
			i->parameters[0]["n_He"] = 0.000001;

			for (short unsigned int j = 1; j < i->parameters.size(); j++)
			{
				i->parameters[j] = i->parameters[0];
			}
		}
	}

	bool tt1 = false;
	// Проверяем наличие всех необходимых переменных
	for (auto& i : this->All_Cell)
	{
		for (auto& num : this->phys_param->param_names)
		{
			if (i->parameters[0].find(num) == i->parameters[0].end())
			{
				//cout << "Parameters: " << num << "    ne opredelen v cells" << endl;
				tt1 = true;
			}
		}
		if (tt1 == true)
		{
			//cout << "Error 4532896514" << endl;
			//exit(-1);
		}
	}

	// Задаём начальные условия на сетке
	if (false)
	{
		for (auto& i : this->All_Cell)
		{
			x = i->center[0][0];
			y = i->center[0][1];
			z = i->center[0][2];
			r = i->func_R(0);

			if (r < this->geo->R0 * 10.0)
			{
				vec << x, y, z;

				cc = this->phys_param->Matr2 * vec;
				the = acos(cc(2) / r);

				BR = -this->phys_param->B_0 * kv(this->phys_param->R_0 / r);
				BPHI = -BR * sin(the) * (r / this->phys_param->R_0);

				dekard_skorost(cc(2), cc(0), cc(1), BR, BPHI, 0.0, V3, V1, V2);

				vv << V1, V2, V3;

				the = -the + const_pi / 2.0;   // Т.к.в данных по СВ на 1 а.е.угол от - 90 до 90 у Алексашова
				cc = this->phys_param->Matr * vv;

				mV = this->phys_param->Get_v_0(the);

				i->parameters[0]["rho"] = this->phys_param->Get_rho_0(the) * pow(this->phys_param->R_0 /r, 2);
				i->parameters[0]["p"] = this->phys_param->p_0 * pow(this->phys_param->R_0 / r, 2 * this->phys_param->gamma);
				i->parameters[0]["Vx"] = mV * vec(0);
				i->parameters[0]["Vy"] = mV * vec(1);
				i->parameters[0]["Vz"] = mV * vec(2);
				i->parameters[0]["Bx"] = cc(0);
				i->parameters[0]["By"] = cc(1);
				i->parameters[0]["Bz"] = cc(2);
				i->parameters[0]["Q"] = i->parameters[0]["rho"];
			}
			else
			{
				i->parameters[0]["rho"] = 1.0;
				i->parameters[0]["p"] = 1.0;
				i->parameters[0]["Vx"] = this->phys_param->Velosity_inf;
				i->parameters[0]["Vy"] = 0.0;
				i->parameters[0]["Vz"] = 0.0;
				i->parameters[0]["Bx"] = -this->phys_param->B_inf * cos(this->phys_param->alphaB_inf);
				i->parameters[0]["By"] = -this->phys_param->B_inf * sin(this->phys_param->alphaB_inf);
				i->parameters[0]["Bz"] = 0.0;
				i->parameters[0]["Q"] = 100.0 * i->parameters[0]["rho"];
			}

			for (short unsigned int j = 1; j < i->parameters.size(); j++)
			{
				i->parameters[j] = i->parameters[0];
			}
		}
	}

	// Задаём граничные условия (на граничных гранях)
	if (true)
	{
		for (auto& i : this->All_boundary_Gran)
		{
			if (i->type == Type_Gran::Outer_Hard)
			{
				i->parameters["rho"] = 1.60063; // 1.0 * (1.0 + this->phys_param->mn_He_inf);
				i->parameters["n_He"] = 0.6; // this->phys_param->mn_He_inf;
				i->parameters["p"] = 1.15; // 1 + (i->parameters["n_He"]) /
					//(i->parameters["rho"] - i->parameters["n_He"]);
				i->parameters["Vx"] = this->phys_param->Velosity_inf;
				i->parameters["Vy"] = 0.0;
				i->parameters["Vz"] = 0.0;
				i->parameters["Bx"] = -this->phys_param->B_inf * cos(this->phys_param->alphaB_inf);
				i->parameters["By"] = -this->phys_param->B_inf * sin(this->phys_param->alphaB_inf);
				i->parameters["Bz"] = 0.0;
				i->parameters["Q"] = 100.0 * i->parameters["rho"];

				i->parameters["rho_H4"] = 1.0;
				i->parameters["Vx_H4"] = this->phys_param->Velosity_inf;
				i->parameters["Vy_H4"] = 0.0;
				i->parameters["Vz_H4"] = 0.0;
				i->parameters["p_H4"] = 0.5;
			}
			else if(i->type == Type_Gran::Inner_Hard)
			{
				x = i->center[0][0];
				y = i->center[0][1];
				z = i->center[0][2];
				r = norm2(x, y, z);

				vec << x, y, z;

				cc = this->phys_param->Matr2 * vec;
				the = acos(cc(2) / r);

				BR = -this->phys_param->B_0 * kv(this->phys_param->R_0 / r);
				BPHI = -BR * sin(the) * (r / this->phys_param->R_0);

				dekard_skorost(cc(2), cc(0), cc(1), BR, BPHI, 0.0, V3, V1, V2);

				vv << V1, V2, V3;

				the = -the + const_pi / 2.0;   // Т.к.в данных по СВ на 1 а.е.угол от - 90 до 90 у Алексашова
				cc = this->phys_param->Matr * vv;

				mV = this->phys_param->Get_v_0(the);

				double Tp = this->phys_param->Get_T_0(the); // Температура

				double np = this->phys_param->Get_rho_0(the);

				double rho = (this->phys_param->mep + 2.0 * this->phys_param->mn_He_0 * 
					this->phys_param->mep + 
					1.0 + 4.0 * this->phys_param->mn_He_0) * np;



				i->parameters["rho"] = rho * pow(this->phys_param->R_0 / r, 2);
				i->parameters["n_He"] = 4.0 * this->phys_param->mn_He_0 * np * pow(this->phys_param->R_0 / r, 2);
				i->parameters["p"] = (1.0 + 3.0 * this->phys_param->mn_He_0 /2.0) * np * Tp *
					pow(this->phys_param->R_0 / r, 2 * this->phys_param->gamma);
				i->parameters["Vx"] = mV * vec(0);
				i->parameters["Vy"] = mV * vec(1);
				i->parameters["Vz"] = mV * vec(2);
				i->parameters["Bx"] = cc(0);
				i->parameters["By"] = cc(1);
				i->parameters["Bz"] = cc(2);
				i->parameters["Q"] = i->parameters["rho"];

				i->parameters["rho_H1"] = 0.00001;
				i->parameters["Vx_H1"] = mV * vec(0);
				i->parameters["Vy_H1"] = mV * vec(1);
				i->parameters["Vz_H1"] = mV * vec(2);
				i->parameters["p_H1"] = 0.00001;
			}
		}
	}

	// Заполняем центральную фиктивную ячейку значениями
	if (true)
	{
		for (auto& num : this->phys_param->param_names)
		{
			this->Cell_Center->parameters[0][num] = 0.0;
		}
		int nk = 0;

		for (auto& i : this->All_boundary_Gran)
		{
			if (i->type != Type_Gran::Inner_Hard) continue;

			nk++;
			for (auto& num : this->phys_param->param_names)
			{
				//this->Cell_Center->parameters[0][num] += i->parameters[num];
				this->Cell_Center->parameters[0][num] += i->cells[0]->parameters[0][num];
			}
		}

		for (auto& num : this->phys_param->param_names)
		{
			this->Cell_Center->parameters[0][num] /= nk;
		}

		this->Cell_Center->parameters[1] = this->Cell_Center->parameters[0];
	}
}

void Setka::Init_TVD(void)
{
	// Эта функция работает после определения нормалей, центра грани и центров ячеек
	// Пробегаемся по всем граням, добавляем ТВД соседей 
	for (auto& gr : this->All_Gran)
	{
		Eigen::Vector3d vec;
		Eigen::Vector3d normalizedCurrent;

		vec << gr->normal[0][0], gr->normal[0][1], gr->normal[0][2];

		auto A = gr->cells[0];
		double maxOppositeScore = 2.0; // Максимальное значение косинуса
		size_t bestIndex = -1;
		size_t Index = -1;
		for (auto& ggr : A->grans)
		{
			Index++;
			if (ggr == gr) continue;

			auto cel = Get_Sosed(A, ggr); // Сосед ячейки A через грань ggr
			if(cel == nullptr) continue;

			normalizedCurrent << cel->center[0][0] - gr->center[0][0],
				cel->center[0][1] - gr->center[0][1], cel->center[0][2] - gr->center[0][2];
			normalizedCurrent = normalizedCurrent.normalized();

			double cosine = vec.dot(normalizedCurrent);

			// Ищем вектор с минимальным косинусом (ближе всего к 180 градусам)
			if (cosine < maxOppositeScore) {
				maxOppositeScore = cosine;
				bestIndex = Index;
			}
		}

		if (this->regim_otladki && bestIndex == -1)
		{
			cout << "Error 7689098909" << endl;
			exit(-1);
		}

		gr->cells_TVD.push_back(Get_Sosed(A, A->grans[bestIndex]));


		if (gr->cells.size() == 1) continue;
		A = gr->cells[1];
		vec = vec * -1.0;
		maxOppositeScore = 2.0; // Максимальное значение косинуса
		bestIndex = -1;
		Index = -1;
		for (auto& ggr : A->grans)
		{
			Index++;
			if (ggr == gr) continue;

			auto cel = Get_Sosed(A, ggr); // Сосед ячейки A через грань ggr
			if (cel == nullptr) continue;

			normalizedCurrent << cel->center[0][0] - gr->center[0][0],
				cel->center[0][1] - gr->center[0][1], cel->center[0][2] - gr->center[0][2];
			normalizedCurrent = normalizedCurrent.normalized();

			double cosine = vec.dot(normalizedCurrent);

			// Ищем вектор с минимальным косинусом (ближе всего к 180 градусам)
			if (cosine < maxOppositeScore) {
				maxOppositeScore = cosine;
				bestIndex = Index;
			}
		}

		if (this->regim_otladki && bestIndex == -1)
		{
			cout << "Error 8764531245" << endl;
			exit(-1);
		}

		gr->cells_TVD.push_back(Get_Sosed(A, A->grans[bestIndex]));

		// Теперь проверим, если ячейка совсем в другой стороне, то удалим её
		vec << gr->normal[0][0], gr->normal[0][1], gr->normal[0][2];
		Eigen::Vector3d k1, k2, k3;
		k1 << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		k2 << gr->cells_TVD[0]->center[0][0],
			gr->cells_TVD[0]->center[0][1], gr->cells_TVD[0]->center[0][2];
		k3 = k2 - k1;
		k3.normalize();
		if (-vec.dot(k3) < 0.1) gr->cells_TVD[0] = nullptr;

		k2 << gr->cells_TVD[1]->center[0][0],
			gr->cells_TVD[1]->center[0][1], gr->cells_TVD[1]->center[0][2];
		k3 = k2 - k1;
		k3.normalize();
		if (vec.dot(k3) < 0.1) gr->cells_TVD[1] = nullptr;
	}
}

void Setka::Calc_sourse_MF(Cell* C, boost::multi_array<double, 2>& SOURSE, 
	short int now, short int zone)
{
	short int N1 = SOURSE.shape()[0];
	short int N2 = SOURSE.shape()[1];

	std::vector<double> U_M_H(N1 - 1);
	std::vector<double> U_H(N1 - 1);
	std::vector<double> sigma(N1 - 1);
	std::vector<double> nu(N1 - 1);
	std::vector<double> kk(N1 - 1);

	for (int i = 0; i < N1 - 1; i++)
	{
		kk[i] = 0.0;
	}

	// Следующая настройка только вручную    &INIT&
	// zone   1, 2, 3, 4
	kk[zone - 1] = 1.0;


	double rho_Th;
	double rho_E; 
	double p_Th;
	double p_Pui;
	double T_Th; 
	double T_E;

	Sootnosheniya(C->parameters[now]["rho"], C->parameters[now]["p"],
		C->parameters[now]["n_He"],
		0.0, 0.0, zone,
		rho_Th, rho_E, p_Th, p_Pui, T_Th, T_E);



	double S1 = 0.0;
	double S2 = 0.0;


	int i = 0;
	for (auto& nam : this->phys_param->H_name)
	{
		U_M_H[i] = sqrt(kv(C->parameters[now]["Vx"] - C->parameters[now]["Vx" + nam])
			+ kv(C->parameters[now]["Vy"] - C->parameters[now]["Vy" + nam])
			+ kv(C->parameters[now]["Vz"] - C->parameters[now]["Vz" + nam])
			+ (64.0 / (9.0 * const_pi)) *
			(p_Th / rho_Th
				+ 2.0 * C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam]));

		U_H[i] = sqrt(kv(C->parameters[now]["Vx"] - C->parameters[now]["Vx" + nam])
			+ kv(C->parameters[now]["Vy"] - C->parameters[now]["Vy" + nam])
			+ kv(C->parameters[now]["Vz"] - C->parameters[now]["Vz" + nam])
			+ (4.0 / const_pi) *
			(p_Th / rho_Th
				+ 2.0 * C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam]));


		sigma[i] = kv(1.0 - this->phys_param->par_a_2 * log(U_M_H[i]));
		nu[i] = rho_Th *
			C->parameters[now]["rho" + nam] * U_M_H[i] * sigma[i];

		i++;
	}

	SOURSE[0][0] = 0.0;

	i = 0;
	for (auto& nam : this->phys_param->H_name)
	{
		SOURSE[0][1] += nu[i] * (C->parameters[now]["Vx" + nam] - C->parameters[now]["Vx"]);
		SOURSE[0][2] += nu[i] * (C->parameters[now]["Vy" + nam] - C->parameters[now]["Vy"]);
		SOURSE[0][3] += nu[i] * (C->parameters[now]["Vz" + nam] - C->parameters[now]["Vz"]);
		SOURSE[0][4] += nu[i] * ((kvv(C->parameters[now]["Vx" + nam], C->parameters[now]["Vy" + nam],
			C->parameters[now]["Vz" + nam]) - kvv(C->parameters[now]["Vx"], C->parameters[now]["Vy"]
				, C->parameters[now]["Vz"])) / 2.0 + (U_H[i] / U_M_H[i]) *
			(2.0 * C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam] - p_Th / rho_Th));
		i++;
	}

	double ddp = (this->phys_param->par_n_H_LISM / this->phys_param->par_Kn);
	SOURSE[0][1] *= ddp;
	SOURSE[0][2] *= ddp;
	SOURSE[0][3] *= ddp;
	SOURSE[0][4] *= ddp;


	i = 0;
	for (auto& nam : this->phys_param->H_name)
	{
		S1 = S1 + nu[i];
		S2 = S2 + nu[i] *
			(kvv(C->parameters[now]["Vx"], C->parameters[now]["Vy"],
				C->parameters[now]["Vz"]) / 2.0
				+ (U_H[i] / U_M_H[i]) * (p_Th / rho_Th));
		i++;
	}

	i = 0;
	for (auto& nam : this->phys_param->H_name)
	{
		double VHx = C->parameters[now]["Vx" + nam];
		double VHy = C->parameters[now]["Vy" + nam];
		double VHz = C->parameters[now]["Vz" + nam];

		SOURSE[i + 1][0] = ddp * (kk[i] * S1 - nu[i]);
		SOURSE[i + 1][1] = ddp * (kk[i] * S1 * C->parameters[now]["Vx"] - nu[i] * C->parameters[now]["Vx" + nam]);
		SOURSE[i + 1][2] = ddp * (kk[i] * S1 * C->parameters[now]["Vy"] - nu[i] * C->parameters[now]["Vy" + nam]);
		SOURSE[i + 1][3] = ddp * (kk[i] * S1 * C->parameters[now]["Vz"] - nu[i] * C->parameters[now]["Vz" + nam]);
		SOURSE[i + 1][4] = ddp * (kk[i] * S2 - nu[i] * (kvv(VHx, VHy, VHz) / 2.0 +
			(U_H[i] / U_M_H[i]) * 2.0 * (C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam])));
		i++;
	}
}

int Setka::determ_zone(Cell* C, short int now)
{
	double rho = C->parameters[now]["rho"];
	double p = C->parameters[now]["p"];
	double u = C->parameters[now]["Vx"];
	double v = C->parameters[now]["Vy"];
	double w = C->parameters[now]["Vz"];
	double M = norm2(u, v, w) / sqrt(this->phys_param->gamma * p / rho);
	

	if (C->parameters[now]["Q"] < 50.0)
	{
		double x = C->center[now][0];
		double y = C->center[now][1];
		double z = C->center[now][2];
		double r = norm2(x, y, z);

		if (M > 1 && r < 45)
		{
			return 1;
		}
		else
		{
			return 2;
		}
	}
	else
	{
		if (M > 1)
		{
			return 4;
		}
		else
		{
			return 3;
		}
	}

	if (this->regim_otladki)
	{
		cout << "Error 7642343620" << endl;
		whach(rho);
		whach(p);
		whach(u);
		whach(v);
		whach(w);
		whach(M);
		whach(C->parameters[now]["Q"]);
		exit(-1);
	}
	return 0;
}

void Setka::Go(bool is_inner_area, size_t steps__, short int metod)
{
	// заполняем коррдинаты узлов на другом временном слое
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}

	whach(this->Gran_TS.size());
	whach(this->Gran_HP.size());
	whach(this->Gran_BS.size());



	this->Test_geometr();
	// Все настройки расчёта считываются из файла Setter.txt
	// Сначала реализовываем расчёт без движения сетки
	unsigned int steps = steps__;

	unsigned short int now1 = 1;
	unsigned short int now2 = 0;

	vector<Gran*>* gran_list;
	vector<Cell*>* cell_list;

	// Если хотим отдельно считать внутреннюю и наружнюю области
	if (false)
	{
		if (is_inner_area == true)
		{
			gran_list = &this->Gran_inner_area;
			cell_list = &this->Cell_inner_area;
		}
		else
		{
			gran_list = &this->Gran_outer_area;
			cell_list = &this->Cell_outer_area;
		}
	}
	else
	{
		gran_list = &this->All_Gran;
		cell_list = &this->All_Cell;
	}

	Cell* A, B;
	double dsr, dsc, dsl;


	double time = 0.000001;  // Текущий шаг по времени
	double loc_time = 0.000001;  // Текущий шаг по времени

	double xc_min = 0.0, yc_min = 0.0, zc_min = 0.0;
	string name_min_time = "___";


	for (unsigned int step = 1; step <= steps; step++)
	{
		if (step % 25 == 0)
		{
			cout << "Global step = " << step << endl;
			whach(time);
			whach(xc_min);
			whach(yc_min);
			whach(zc_min);
			whach(name_min_time);
			cout << "__________________" << endl;
		}

		now2 = now1;
		now1 = (now1 + 1)%2;

		time = loc_time;
		loc_time = 100000000.0;

		

		//omp_set_num_threads(1); // 32


		// Здесь задаётся выходящее условие для водорода 1-3 на Type_Gran::Outer_Hard
		if (false)
		{
			for (auto& gran : this->All_boundary_Gran)
			{
				if (gran->type == Type_Gran::Outer_Hard)
				{
					double u1, v1, w1;
					for (auto& nam : this->phys_param->H_name)
					{
						if (nam == "_H4") continue;

						u1 = gran->cells[0]->parameters[now1]["Vx" + nam];
						v1 = gran->cells[0]->parameters[now1]["Vy" + nam];
						w1 = gran->cells[0]->parameters[now1]["Vz" + nam];
						if (u1 * gran->normal[now1][0] + v1 * gran->normal[now1][1] +
							w1 * gran->normal[now1][2] > 0.0)
						{
							gran->parameters["rho" + nam] = gran->cells[0]->parameters[now1]["rho" + nam];
							gran->parameters["p" + nam] = gran->cells[0]->parameters[now1]["p" + nam];
							gran->parameters["Vx" + nam] = u1;
							gran->parameters["Vy" + nam] = v1;
							gran->parameters["Vz" + nam] = w1;
						}
						else
						{
							gran->parameters["rho" + nam] = gran->cells[0]->parameters[now1]["rho" + nam];
							gran->parameters["p" + nam] = gran->cells[0]->parameters[now1]["p" + nam];
							gran->parameters["Vx" + nam] = 0.1 * gran->normal[now1][0];
							gran->parameters["Vy" + nam] = 0.1 * gran->normal[now1][1];
							gran->parameters["Vz" + nam] = 0.1 * gran->normal[now1][2];
						}
					}
				}
			}
		}
		

		// Считаем скорости граней и сразу передвигаем опорные узлы
		if (is_inner_area == false)
		{
			this->Culc_Velocity_surface(now1, time, 1);

			// Перестраиваем сетку
			for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
			{
				auto lu = this->All_Luch[i_step];
				lu->dvigenie(now2);
			}
			this->Calculating_measure(now2);
		}

		// Расчитываем потоки через грани
		// в private не добавляются нормально vectora, надо либо обычные массивы делать, либо 
		// создавать их внутри в каждом потоке
#pragma omp parallel for reduction(min:loc_time)
		for(int i_step = 0; i_step < gran_list->size(); i_step++)
		{
			//whach(GG->parameters["rho_H4"]);
			// omp_get_thread_num()  - номер нити
			// omp_get_num_threads() - всего потоков
			//cout << "i_step = " << i_step << "  potok = " << omp_get_thread_num() << "  Vsego = " << 
			//	omp_get_num_threads() << endl;

			auto& gran = (*gran_list)[i_step];

			string nmnm;
			double ntnt = this->Culc_Gran_Potok(gran, now1, metod, nmnm, time);  // Считает потоки через данную грань (записывает результат в параметры грани)

			if (ntnt < loc_time)
			{
				#pragma omp critical 
				{
					if (ntnt < loc_time)
					{
						loc_time = min(loc_time, ntnt);
						xc_min = gran->center[now1][0];
						yc_min = gran->center[now1][1];
						zc_min = gran->center[now1][2];
						name_min_time = nmnm;
					}
				}
			}
		}

#pragma omp barrier
		//cout << "barrier" << endl;

		bool print_p_less_0 = false;

		// Расчитываем законы сохранения в ячейках
#pragma omp parallel for
		for (size_t i_step = 0; i_step < cell_list->size(); i_step++)
		{
			auto& cell = (*cell_list)[i_step];
			double Volume = cell->volume[now1];
			double Volume2 = cell->volume[now2];
			int zone;

			std::vector<double> POTOK;
			POTOK.resize(9);

			std::vector<double> POTOK2; // Для конвективных параметров
			POTOK2.resize(2);

			boost::multi_array<double, 2> POTOK_F(boost::extents[this->phys_param->H_name.size()][5]);
			for (short int i = 0; i < this->phys_param->H_name.size(); ++i) {
				for (short int j = 0; j < 5; ++j) {
					POTOK_F[i][j] = 0.0;
				}
			}

			boost::multi_array<double, 2> SOURSE(boost::extents[this->phys_param->H_name.size() + 1][5]);
			
			for (short int i = 0; i < this->phys_param->H_name.size() + 1; ++i) 
			{
				for (short int j = 0; j < 5; ++j) {
					SOURSE[i][j] = 0.0;
				}
			}

			POTOK[0] = POTOK[1] = POTOK[2] = POTOK[3] = POTOK[4] = POTOK[5] = 
				POTOK[6] = POTOK[7] = POTOK[8] = 0.0;

			POTOK2[0] = 0.0;
			POTOK2[1] = 0.0;

			for (auto& gran : cell->grans)
			{
				short int sign_potok = 1;
				if (gran->cells[0] != cell) sign_potok = -1;

				if (this->phys_param->culc_plasma == true)
				{
					POTOK[0] += sign_potok * gran->parameters["Pm"];
					POTOK[1] += sign_potok * gran->parameters["PVx"];
					POTOK[2] += sign_potok * gran->parameters["PVy"];
					POTOK[3] += sign_potok * gran->parameters["PVz"];
					POTOK[4] += sign_potok * gran->parameters["Pe"];
					POTOK[5] += sign_potok * gran->parameters["PBx"];
					POTOK[6] += sign_potok * gran->parameters["PBy"];
					POTOK[7] += sign_potok * gran->parameters["PBz"];
					POTOK[8] += sign_potok * gran->parameters["PdivB"];

					POTOK2[0] += sign_potok * gran->parameters["PQ"];
					POTOK2[1] += sign_potok * gran->parameters["Pn_He"];
				}

				for (short int i = 0; i < this->phys_param->H_name.size(); i++)
				{
					POTOK_F[i][0] += sign_potok * gran->parameters["Pm" + this->phys_param->H_name[i]];
					POTOK_F[i][1] += sign_potok * gran->parameters["PVx" + this->phys_param->H_name[i]];
					POTOK_F[i][2] += sign_potok * gran->parameters["PVy" + this->phys_param->H_name[i]];
					POTOK_F[i][3] += sign_potok * gran->parameters["PVz" + this->phys_param->H_name[i]];
					POTOK_F[i][4] += sign_potok * gran->parameters["Pe" + this->phys_param->H_name[i]];
				}
			}

			zone = this->determ_zone(cell, now1);
			this->Calc_sourse_MF(cell, SOURSE, now1, zone);

			double rho3, u3, v3, w3, bx3, by3, bz3, p3, Q3, n_He3, n_He;
			double rho, vx, vy, vz, p, bx, by, bz, dsk, Q;

			if (this->phys_param->culc_plasma == true)
			{
				rho = cell->parameters[now1]["rho"];
				if (cell->parameters[now1].find("Q") != cell->parameters[now1].end())
				{
					Q = cell->parameters[now1]["Q"];
				}
				else
				{
					Q = 0.0;
				}

				if (cell->parameters[now1].find("n_He") != cell->parameters[now1].end())
				{
					n_He = cell->parameters[now1]["n_He"];
				}
				else
				{
					n_He = 0.0;
				}



				vx = cell->parameters[now1]["Vx"];
				vy = cell->parameters[now1]["Vy"];
				vz = cell->parameters[now1]["Vz"];
				p = cell->parameters[now1]["p"];
				bx = cell->parameters[now1]["Bx"];
				by = cell->parameters[now1]["By"];
				bz = cell->parameters[now1]["Bz"];
				dsk = scalarProductFast(vx, vy, vz, bx, by, bz);

				rho3 = rho * Volume / Volume2 - time * POTOK[0] / Volume2;

				if (cell->parameters[now1].find("Q") != cell->parameters[now1].end())
				{
					Q3 = Q * Volume / Volume2 - time * POTOK2[0] / Volume2;
				}
				else
				{
					Q3 = 0.0;
				}

				if (cell->parameters[now1].find("n_He") != cell->parameters[now1].end())
				{
					n_He3 = n_He * Volume / Volume2 - time * POTOK2[1] / Volume2;
				}
				else
				{
					n_He3 = 0.0;
				}

				if (rho3 < 0.00000001)
				{
					rho3 = 0.0001;
					Q3 = Q / rho * rho3;
					cout << "Plasma  rho < 0" << endl;
				}

				u3 = (rho * vx * Volume / Volume2 - time * (POTOK[1] + (bx / cpi4) * POTOK[8]) / Volume2
					+ time * SOURSE[0][1]) / rho3;
				v3 = (rho * vy * Volume / Volume2 - time * (POTOK[2] + (by / cpi4) * POTOK[8]) / Volume2
					+ time * SOURSE[0][2]) / rho3;
				w3 = (rho * vz * Volume / Volume2 - time * (POTOK[3] + (bz / cpi4) * POTOK[8]) / Volume2
					+ time * SOURSE[0][3]) / rho3;

				bx3 = bx * Volume / Volume2 - time * (POTOK[5] + vx * POTOK[8]) / Volume2;
				by3 = by * Volume / Volume2 - time * (POTOK[6] + vy * POTOK[8]) / Volume2;
				bz3 = bz * Volume / Volume2 - time * (POTOK[7] + vz * POTOK[8]) / Volume2;

				p3 = (((p / this->phys_param->g1 + 0.5 * rho * kvv(vx, vy, vz) + kvv(bx, by, bz) / 25.13274122871834590768) * Volume / Volume2
					- time * (POTOK[4] + (dsk / cpi4) * POTOK[8]) / Volume2 + time * SOURSE[0][4]) -
					0.5 * rho3 * kvv(u3, v3, w3) - kvv(bx3, by3, bz3) / 25.13274122871834590768) * this->phys_param->g1;


				if (std::isnan(rho3) || std::fpclassify(rho3) == FP_SUBNORMAL || 
					std::isnan(Q3) || std::fpclassify(Q3) == FP_SUBNORMAL)
				{
					cout << "Error  9851234578" << endl;
					whach(Volume);
					whach(Volume2);
					whach(rho);
					whach(time);
					whach(Q);
					whach(n_He);

					for (short unsigned int ik = 0; ik < 9; ik++)
					{
						cout << "ik: " << POTOK[ik] << endl;
					}
					exit(-1);
				}


				if (p3 < 0.00000001)
				{
					p3 = 0.000001;
					if (step % 25 == 0 && print_p_less_0 == false)
					{
						cout << "Plasma  p < 0" << endl;
						cout << cell->center[now2][0] << " " <<
							cell->center[now2][1] << " " <<
							cell->center[now2][2] << endl;
					}
					print_p_less_0 = true;
				}
				

				cell->parameters[now2]["rho"] = rho3;
				if (cell->parameters[now1].find("Q") != cell->parameters[now1].end())
				{
					cell->parameters[now2]["Q"] = Q3;
				}
				if (cell->parameters[now1].find("n_He") != cell->parameters[now1].end())
				{
					cell->parameters[now2]["n_He"] = n_He3;
				}
				cell->parameters[now2]["Vx"] = u3;
				cell->parameters[now2]["Vy"] = v3;
				cell->parameters[now2]["Vz"] = w3;
				cell->parameters[now2]["Bx"] = bx3;
				cell->parameters[now2]["By"] = by3;
				cell->parameters[now2]["Bz"] = bz3;
				cell->parameters[now2]["p"] = p3;
			}

			/*if (norm2(cell->center[now1][0], cell->center[now1][1], cell->center[now1][2]) < 3.0)
			{
				SOURSE[2][0] = 0.0;
				SOURSE[2][1] = 0.0;
				SOURSE[2][2] = 0.0;
				SOURSE[2][3] = 0.0;
				SOURSE[2][4] = 0.0;
			}*/

			// Теперь считаем для остальных жидкостей
			int i = 0;
			for (auto& nam : this->phys_param->H_name)
			{
				rho = cell->parameters[now1]["rho" + nam];
				vx = cell->parameters[now1]["Vx" + nam];
				vy = cell->parameters[now1]["Vy" + nam];
				vz = cell->parameters[now1]["Vz" + nam];
				p = cell->parameters[now1]["p" + nam];

				rho3 = rho * Volume / Volume2 - time * POTOK_F[i][0] / Volume2 + time * SOURSE[i + 1][0];

				if (rho3 < 0.00000001)
				{
					rho3 = 0.00000001;
					//cout << "Hidrogen  rho < 0  " << nam << endl;
				}

				u3 = (rho * vx * Volume / Volume2 - time * (POTOK_F[i][1]) / Volume2
					+ time * SOURSE[i + 1][1]) / rho3;
				v3 = (rho * vy * Volume / Volume2 - time * (POTOK_F[i][2]) / Volume2
					+ time * SOURSE[i + 1][2]) / rho3;
				w3 = (rho * vz * Volume / Volume2 - time * (POTOK_F[i][3]) / Volume2
					+ time * SOURSE[i + 1][3]) / rho3;


				p3 = (((p / this->phys_param->g1 + 0.5 * rho * kvv(vx, vy, vz)) * Volume / Volume2
					- time * (POTOK_F[i][4]) / Volume2 + time * SOURSE[i + 1][4]) -
					0.5 * rho3 * kvv(u3, v3, w3)) * this->phys_param->g1;

				if (p3 < 0.00000001)
				{
					p3 = 0.00000001;
					//cout << "Hidrogen  p < 0  " << nam << endl;
					//cout << "Center = " << cell->center[now2][0] << " " <<
					//	cell->center[now2][1] << " " <<
					//	cell->center[now2][2] << endl;
				}

				cell->parameters[now2]["rho" + nam] = rho3;
				cell->parameters[now2]["Vx" + nam] = u3;
				cell->parameters[now2]["Vy" + nam] = v3;
				cell->parameters[now2]["Vz" + nam] = w3;
				cell->parameters[now2]["p" + nam] = p3;

				i++;
			}
		}

#pragma omp barrier

		// Считаем фиктивную центральную ячейку
		if(true) 
			{
				auto& cell = this->Cell_Center;
				double Volume = 4.0 * const_pi * kyb(this->geo->R0)/3.0;

				boost::multi_array<double, 2> POTOK_F(boost::extents[this->phys_param->H_name.size()][5]);
				for (short int i = 0; i < this->phys_param->H_name.size(); ++i) {
					for (short int j = 0; j < 5; ++j) {
						POTOK_F[i][j] = 0.0;
					}
				}

				for (auto& gran : this->All_boundary_Gran)
				{
					if (gran->type != Type_Gran::Inner_Hard) continue;

					short int sign_potok = -1;

					for (short int i = 1; i < this->phys_param->H_name.size(); i++)
					{
						POTOK_F[i][0] += sign_potok * gran->parameters["Pm" + this->phys_param->H_name[i]];
						POTOK_F[i][1] += sign_potok * gran->parameters["PVx" + this->phys_param->H_name[i]];
						POTOK_F[i][2] += sign_potok * gran->parameters["PVy" + this->phys_param->H_name[i]];
						POTOK_F[i][3] += sign_potok * gran->parameters["PVz" + this->phys_param->H_name[i]];
						POTOK_F[i][4] += sign_potok * gran->parameters["Pe" + this->phys_param->H_name[i]];
					}
				}

				double rho3, u3, v3, w3, bx3, by3, bz3, p3, Q3;
				double rho, vx, vy, vz, p, bx, by, bz, dsk, Q;

				// Теперь считаем для остальных жидкостей
				int i = 1;
				for (auto& nam : this->phys_param->H_name)
				{
					if (nam == "_H1") continue;

					rho = cell->parameters[now1]["rho" + nam];
					vx = cell->parameters[now1]["Vx" + nam];
					vy = cell->parameters[now1]["Vy" + nam];
					vz = cell->parameters[now1]["Vz" + nam];
					p = cell->parameters[now1]["p" + nam];

					rho3 = rho - time * POTOK_F[i][0] / Volume;

					if (rho3 < 0.0001)
					{
						rho3 = 0.0001;
					}

					u3 = (rho * vx - time * (POTOK_F[i][1]) / Volume) / rho3;
					v3 = (rho * vy - time * (POTOK_F[i][2]) / Volume) / rho3;
					w3 = (rho * vz - time * (POTOK_F[i][3]) / Volume) / rho3;


					p3 = (((p / this->phys_param->g1 + 0.5 * rho * kvv(vx, vy, vz))
						- time * (POTOK_F[i][4]) / Volume) -
						0.5 * rho3 * kvv(u3, v3, w3)) * this->phys_param->g1;

					if (p3 < 0.0001)
					{
						p3 = 0.0001;
					}

					cell->parameters[now2]["rho" + nam] = rho3;
					cell->parameters[now2]["Vx" + nam] = u3;
					cell->parameters[now2]["Vy" + nam] = v3;
					cell->parameters[now2]["Vz" + nam] = w3;
					cell->parameters[now2]["p" + nam] = p3;

					i++;
				}
			}

	}
}

double Setka::Culc_Gran_Potok(Gran* gr, unsigned short int now, short int metod, string& name, const double& time)
{
	double dsr, dsc, dsl;
	std::vector<double> qqq, qqq1, qqq2;
	qqq.resize(8);
	qqq1.resize(8);
	qqq2.resize(8);
	std::vector<double> konvect_left, konvect_right, konvect;
	PrintOptions Option = PrintOptions{};
	double area = gr->area[now];
	name = "______";
	short int metod_ = metod;

	unordered_map<string, double> par_left, par_right;
	double w = gr->culc_velosity(now, time);
	double loc_time = 10000000.0;

	double dist;
	if (gr->type == Type_Gran::Us)
	{
		// Это для расчёта шага по времени
		auto A = gr->cells[0];
		auto B = gr->cells[1];
		dist = norm2(A->center[now][0] - B->center[now][0],
			A->center[now][1] - B->center[now][1],
			A->center[now][2] - B->center[now][2]) / 2.0;
	}

	this->Snos_on_Gran(gr, par_left, par_right, now);

	Option.x = gr->center[now][0];
	Option.y = gr->center[now][1];
	Option.z = gr->center[now][2];

	if (this->phys_param->culc_plasma == true)
	{
		Option.x = gr->center[now][0];
		Option.y = gr->center[now][1];
		Option.z = gr->center[now][2];
		Option.fluid = "plasma";

		qqq1[0] = par_left["rho"];
		qqq1[1] = par_left["Vx"];
		qqq1[2] = par_left["Vy"];
		qqq1[3] = par_left["Vz"];
		qqq1[4] = par_left["p"];
		qqq1[5] = par_left["Bx"];
		qqq1[6] = par_left["By"];
		qqq1[7] = par_left["Bz"];

		qqq2[0] = par_right["rho"];
		qqq2[1] = par_right["Vx"];
		qqq2[2] = par_right["Vy"];
		qqq2[3] = par_right["Vz"];
		qqq2[4] = par_right["p"];
		qqq2[5] = par_right["Bx"];
		qqq2[6] = par_right["By"];
		qqq2[7] = par_right["Bz"];


		if (par_left.find("Q") != par_left.end())
		{
			konvect_left.push_back(par_left["Q"]);
			konvect_right.push_back(par_right["Q"]);
			konvect.push_back(0.0);
		}

		if (par_left.find("n_He") != par_left.end())
		{
			konvect_left.push_back(par_left["n_He"]);
			konvect_right.push_back(par_right["n_He"]);
			konvect.push_back(0.0);
		}

		// Если это контакт, записываем магнитное давление в обычное
		// И удаляем магнитные поля
		if (gr->type2 == Type_Gran_surf::HP && this->phys_param->bn_in_p_on_HP == true)
		{
			if (metod_ == 3) metod_ = 2;

			qqq1[4] += kvv(qqq1[5], qqq1[6], qqq1[7]) / (8.0 * const_pi);
			qqq1[5] = qqq1[6] = qqq1[7] = 0.0;

			qqq2[4] += kvv(qqq2[5], qqq2[6], qqq2[7]) / (8.0 * const_pi);
			qqq2[5] = qqq2[6] = qqq2[7] = 0.0;
		}

		metod_ = gr->Get_method();


		//metod_ = 2;
		
		if (gr->type2 == Type_Gran_surf::HP && this->phys_param->bn_in_p_on_HP == true)
		{
			std::vector<double> n(3);
			n[0] = gr->normal[now][0];
			n[1] = gr->normal[now][1];
			n[2] = gr->normal[now][2];

			this->phys_param->Godunov_Solver_Alexashov(qqq1, qqq2,//
				n, qqq, dsl, dsr, dsc, w, true);

			qqq[5] = qqq[6] = qqq[7] = 0.0;
			konvect[0] = 0.0;
			konvect[1] = 0.0;
		}
		else
		{
			bool left_ydar = false;
			if (gr->type2 == Type_Gran_surf::TS) left_ydar = true;

			this->phys_param->chlld(metod_, gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 1,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option, left_ydar);
		}


		if (gr->type == Type_Gran::Us)
		{
			double dnt = this->phys_param->KFL * dist
				/ (max(fabs(dsl), fabs(dsr)) + fabs(w));

			if (dnt < loc_time)
			{
				loc_time = min(loc_time, dnt);
				name = "plasma";
			}
		}

		gr->parameters["Pm"] = qqq[0] * area;
		gr->parameters["PVx"] = qqq[1] * area;
		gr->parameters["PVy"] = qqq[2] * area;
		gr->parameters["PVz"] = qqq[3] * area;
		gr->parameters["Pe"] = qqq[4] * area;
		gr->parameters["PBx"] = qqq[5] * area;
		gr->parameters["PBy"] = qqq[6] * area;
		gr->parameters["PBz"] = qqq[7] * area;
		gr->parameters["PdivB"] = 0.5 * scalarProductFast(gr->normal[now][0],
			gr->normal[now][1], gr->normal[now][2],
			qqq1[5] + qqq2[5], qqq1[6] + qqq2[6], qqq1[7] + qqq2[7]) * area;

		// Для контакта поток Bn равен нулю
		if (gr->type2 == Type_Gran_surf::HP && this->phys_param->bn_in_p_on_HP == true)
		{
			gr->parameters["PdivB"] = 0.0;
		}

		if (par_left.find("Q") != par_left.end())
		{
			gr->parameters["PQ"] = konvect[0] * area; // Может быть неправльный порядок при других конвективных переменных
		}

		if (par_left.find("n_He") != par_left.end())
		{
			gr->parameters["Pn_He"] = konvect[1] * area; // Может быть неправльный порядок при других конвективных переменных
		}
	}

	konvect_left.clear();
	konvect_right.clear();
	konvect.clear();

	for (auto& nam : this->phys_param->H_name)
	{
		Option.fluid = nam;
		qqq1[0] = par_left["rho" + nam];
		qqq1[1] = par_left["Vx" + nam];
		qqq1[2] = par_left["Vy" + nam];
		qqq1[3] = par_left["Vz" + nam];
		qqq1[4] = par_left["p" + nam];
		qqq1[5] = 0.0;
		qqq1[6] = 0.0;
		qqq1[7] = 0.0;

		qqq2[0] = par_right["rho" + nam];
		qqq2[1] = par_right["Vx" + nam];
		qqq2[2] = par_right["Vy" + nam];
		qqq2[3] = par_right["Vz" + nam];
		qqq2[4] = par_right["p" + nam];
		qqq2[5] = 0.0;
		qqq2[6] = 0.0;
		qqq2[7] = 0.0;

		// metod
		this->phys_param->chlld(0, gr->normal[now][0], gr->normal[now][1],
			gr->normal[now][2],
			w, qqq1, qqq2, qqq, false, 1,
			konvect_left, konvect_right, konvect, dsr, dsc, dsl,
			Option);

		
		if (gr->type == Type_Gran::Us)
		{
			double dnt = this->phys_param->KFL * dist
				/ (max(fabs(dsl), fabs(dsr)) + fabs(w));

			if (dnt < loc_time)
			{
				loc_time = min(loc_time, dnt);
				name = nam;
			}
		}

		gr->parameters["Pm" + nam] = qqq[0] * area;
		gr->parameters["PVx" + nam] = qqq[1] * area;
		gr->parameters["PVy" + nam] = qqq[2] * area;
		gr->parameters["PVz" + nam] = qqq[3] * area;
		gr->parameters["Pe" + nam] = qqq[4] * area;

	}


	return loc_time;
}

void Setka::Save_for_interpolate(string filename)
{
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		cout << "Error 097564537  Can not open file to writing: " + filename << endl;
		exit(-1);
	}

	// Записываем количество строк
	size_t size = this->phys_param->param_names.size();
	out.write(reinterpret_cast<const char*>(&size), sizeof(size));

	// Записываем каждую строку
	for (const auto& str : this->phys_param->param_names) {
		// Сначала записываем длину строки
		size_t str_size = str.size();
		out.write(reinterpret_cast<const char*>(&str_size), sizeof(str_size));
		// Затем саму строку
		out.write(str.data(), str_size);
	}

	// Считаем сколько дополнительных ячеек будет на внешней границе
	unsigned int gr_b = 0;
	for (const auto& gr : this->All_boundary_Gran)
	{
		if (gr->type != Type_Gran::Inner_Hard)
		{
			gr_b++;
		}
	}

	// Записываем количество ячеек
	size = this->All_Cell.size() + gr_b;
	out.write(reinterpret_cast<const char*>(&size), sizeof(size));

	for (const auto& Cel : this->All_Cell)
	{
		double aa = Cel->center[0][0];
		double bb = Cel->center[0][1];
		double cc = Cel->center[0][2];
		out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
		out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
		out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

		for (const auto& i : this->phys_param->param_names)
		{
			aa = 0.0;

			if (Cel->parameters[0].find(i) != Cel->parameters[0].end())
			{
				aa = Cel->parameters[0][i];
			}
			
			out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
		}
	}

	// Записываем дополнительные точки (на небольшом удалении от внешней границы, чтоб 
	// убрать артефакты в интерполяции
	Eigen::Vector3d C1, C2, C3, C4;
	for (const auto& gr : this->All_boundary_Gran)
	{
		if (gr->type != Type_Gran::Inner_Hard)
		{
			auto A = gr->cells[0];
			C1 << A->center[0][0], A->center[0][1], A->center[0][2];
			C2 << gr->center[0][0], gr->center[0][1], gr->center[0][2];
			C3 = C2 - C1;
			C4 = C2 + C3;
			double aa = C4[0];
			double bb = C4[1];
			double cc = C4[2];
			out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
			out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

			for (const auto& i : this->phys_param->param_names)
			{
				aa = 0.0;
				if (A->parameters[0].find(i) != A->parameters[0].end())
				{
					aa = A->parameters[0][i];
				}

				out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
			}
		}
	}

	// Записываем центральную точку
	if (true)
	{
		double aa = 0.0;
		double bb = 0.0;
		double cc = 0.0;
		out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
		out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
		out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

		for (const auto& i : this->phys_param->param_names)
		{
			aa = 0.0;
			if (this->Cell_Center->parameters[0].find(i) != this->Cell_Center->parameters[0].end())
			{
				aa = this->Cell_Center->parameters[0][i];
			}
			out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
		}
	}

	for (size_t i = 0; i < 1000; i++)
	{
		bool aa = false;
		out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
	}

	out.close();
}

void Setka::Save_cell_parameters(string filename)
{
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		cout << "Error 097564537  Can not open file to writing: " + filename << endl;
		exit(-1);
	}

	for (auto& i : this->All_Cell)
	{
		// Записываем количество элементов
		size_t size = i->parameters[0].size();
		out.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

		// Записываем каждую пару ключ-значение
		for (const auto& pair : i->parameters[0]) {
			// Сначала записываем длину ключа и сам ключ
			size_t key_size = pair.first.size();
			out.write(reinterpret_cast<const char*>(&key_size), sizeof(size_t));
			out.write(pair.first.c_str(), key_size);

			// Затем записываем значение
			out.write(reinterpret_cast<const char*>(&pair.second), sizeof(double));
		}
	}

	bool bb;

	
	// Записываем координаты узлов
	bb = true;
	out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
	size_t size = this->All_Yzel.size();
	out.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
	for (auto& yz : this->All_Yzel)
	{
		double x = yz->coord[0][0];
		double y = yz->coord[0][1];
		double z = yz->coord[0][2];
		out.write(reinterpret_cast<const char*>(&x), sizeof(double));
		out.write(reinterpret_cast<const char*>(&y), sizeof(double));
		out.write(reinterpret_cast<const char*>(&z), sizeof(double));
	}



	// Записываем для будующих считываний
	bb = false;
	for (int i = 0; i < 1000; i++)
	{

		out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
	}
}

void Setka::Save_cell_MK_parameters(string filename)
{
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		cout << "Error 097564537  Can not open file to writing: " + filename << endl;
		exit(-1);
	}

	size_t size = this->phys_param->MK_param.size();
	out.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

	for (const auto& pair : this->phys_param->MK_param)
	{
		size_t key_size = pair.size();
		out.write(reinterpret_cast<const char*>(&key_size), sizeof(size_t));
		out.write(pair.c_str(), key_size);
	}

	for (auto& i : this->All_Cell)
	{
		// Записываем каждую пару ключ-значение
		for (const auto& pair : this->phys_param->MK_param)
		{
			if (i->parameters[0].find(pair) == i->parameters[0].end())
			{
				cout << "Error 8765509090" << endl;
				cout << pair << endl;
				cout << "-------------------" << endl;
				i->parameters[0][pair] = 0.0;
			}

			// записываем значение
			out.write(reinterpret_cast<const char*>(&i->parameters[0][pair]), sizeof(double));
		}
	}

	bool bb;

	// Записываем для будующих считываний
	bb = false;
	for (int i = 0; i < 1000; i++)
	{
		out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
	}
}

void Setka::Download_cell_parameters(string filename)
{
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		cout << "Error 6545478564  Can not open file to reading: " + filename << endl;
		exit(-1);
	}


	for (auto& ii : this->All_Cell)
	{
		// Читаем количество элементов
		size_t size;
		in.read(reinterpret_cast<char*>(&size), sizeof(size_t));

		for (size_t i = 0; i < size; ++i) {
			// Читаем ключ
			size_t key_size;
			in.read(reinterpret_cast<char*>(&key_size), sizeof(size_t));

			std::vector<char> key_buffer(key_size);
			in.read(key_buffer.data(), key_size);
			std::string key(key_buffer.begin(), key_buffer.end());

			// Читаем значение
			double value;
			in.read(reinterpret_cast<char*>(&value), sizeof(double));

			if (std::find(this->phys_param->MK_param.begin(),
				this->phys_param->MK_param.end(), key) !=
				this->phys_param->MK_param.end())
			{
				ii->parameters[0][key] = 0.0;
			}
			else
			{
				ii->parameters[0][key] = value;
				ii->parameters[1][key] = value;
			}
		}
	}


	bool bb;

	// Читаем координаты узлов
	in.read(reinterpret_cast<char*>(&bb), sizeof(bb));
	if (bb == true)
	{
		size_t size;
		in.read(reinterpret_cast<char*>(&size), sizeof(size_t));
		double x, y, z;

		// Чтение координат каждого узла
		for (size_t i = 0; i < size; ++i) 
		{
			in.read(reinterpret_cast<char*>(&x), sizeof(double));
			in.read(reinterpret_cast<char*>(&y), sizeof(double));
			in.read(reinterpret_cast<char*>(&z), sizeof(double));

			// Заполняем координаты узла
			auto yz = this->All_Yzel[i];
			yz->coord[0][0] = x;
			yz->coord[0][1] = y;
			yz->coord[0][2] = z;

			yz->coord[1][0] = x;
			yz->coord[1][1] = y;
			yz->coord[1][2] = z;
		}
	}


	in.close();
}

void Setka::Download_cell_MK_parameters(string filename, short int zone_except)
{
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		cout << "Error 6545478564  Can not open file to reading: " + filename << endl;
		exit(-1);
	}

	vector<string> param_file;
	cout << "MK_parameters: ";
	size_t size;
	in.read(reinterpret_cast<char*>(&size), sizeof(size_t));
	for (size_t i = 0; i < size; ++i)
	{
		// Читаем ключ
		size_t key_size;
		in.read(reinterpret_cast<char*>(&key_size), sizeof(size_t));

		std::vector<char> key_buffer(key_size);
		in.read(key_buffer.data(), key_size);
		std::string key(key_buffer.begin(), key_buffer.end());

		cout << key << " ";
		param_file.push_back(key);
	}
	cout << endl;


	for (auto& ii : this->All_Cell)
	{

		for (const auto& pair : param_file)
		{
			// Читаем значение
			double value;
			in.read(reinterpret_cast<char*>(&value), sizeof(double));

			if (ii->MK_zone != zone_except)
			{
				ii->parameters[0][pair] = value;
			}
		}
	}


	bool bb;


	in.close();
}

void Setka::Culc_rotors_in_cell(void)
{
	cout << "Start: Culc_rotor_in_cell" << endl;

	this->phys_param->param_names.push_back("rotB/b2_x");
	this->phys_param->param_names.push_back("rotB/b2_y");
	this->phys_param->param_names.push_back("rotB/b2_z");
	// Добавили переменную для интерполяции
	unsigned int k1 = 0;

#pragma omp parallel for schedule(dynamic)
	for (auto& cell : this->All_Cell)
	{
		#pragma omp critical (first) 
		{
			k1++;
			if (k1 % 10000 == 0)
			{
				cout << "Gran = " << k1 << "    Iz: " << this->All_Cell.size() << endl;
			}
		}

		int n = 0;  // Число элементов матрицы или граней в ячейке
		n = cell->grans.size();
		Eigen::MatrixXd M(n, 3);
		Eigen::VectorXd F(n);

		short int ig = -1;
		for (auto gr : cell->grans)
		{
			ig++;
			Eigen::Vector3d center_gr;
			Eigen::Vector3d normal;
			double B_on_gran = 0.0;
			center_gr << gr->center[0][0], gr->center[0][1], gr->center[0][2];
			normal << gr->normal[0][0], gr->normal[0][1], gr->normal[0][2];

			if (gr->cells[0] != cell) normal = -normal;


			for (auto ed : gr->edges)
			{
				Eigen::Vector3d ll;
				Eigen::Vector3d ed_center;
				Eigen::Vector3d Vec;
				Eigen::Vector3d Vec_all;
				int v_all = 0;

				Vec_all << 0.0, 0.0, 0.0;

				ll << (ed->A->coord[0][0] - ed->B->coord[0][0]),
					(ed->A->coord[0][1] - ed->B->coord[0][1]),
					(ed->A->coord[0][2] - ed->B->coord[0][2]);

				ed_center << (ed->A->coord[0][0] + ed->B->coord[0][0]) / 2.0,
					(ed->A->coord[0][1] + ed->B->coord[0][1]) / 2.0,
					(ed->A->coord[0][2] + ed->B->coord[0][2]) / 2.0;

				if (true)
				{
					for (auto gr_ : ed->grans)
					{
						if (gr_->cells.size() == 1)
						{
							if (gr_->cells[0]->type != cell->type)
							{
								continue;
							}
						}
						else
						{
							if (gr_->cells[0]->type != cell->type &&
								gr_->cells[1]->type != cell->type)
							{
								continue;
							}
						}

						unordered_map<string, double> par_left, par_right;
						if (gr_->type == Type_Gran::Us)
						{
							this->Snos_on_Gran(gr_, par_left, par_right, 0);

							if (gr_->type2 == Type_Gran_surf::Us)
							{
								Vec << (par_left["Bx"] + par_right["Bx"]) / 2.0,
									(par_left["By"] + par_right["By"]) / 2.0,
									(par_left["Bz"] + par_right["Bz"]) / 2.0;
							}
							else
							{
								if (gr_->cells[0]->type == cell->type)
								{
									Vec << par_right["Bx"],
										par_right["By"],
										par_right["Bz"];
								}
								else
								{
									Vec << par_left["Bx"],
										par_left["By"],
										par_left["Bz"];
								}
							}
						}
						else
						{
							Vec << cell->parameters[0]["Bx"], cell->parameters[0]["By"], cell->parameters[0]["Bz"];
						}

						v_all++;
						Vec_all += Vec;
					}
				}
				else
				{
					for (auto ce : ed->cells)
					{
						v_all++;
						Vec_all[0] += ce->parameters[0]["Bx"];
						Vec_all[1] += ce->parameters[0]["By"];
						Vec_all[2] += ce->parameters[0]["Bz"];
					}
				}


				Vec_all /= v_all;

				Eigen::Vector3d nn = (ed_center - center_gr).cross(ll);

				double norm = Vec_all.norm();
				Vec_all /= (kv(norm));

				if (nn.dot(normal) > 0)
				{
					B_on_gran += Vec_all.dot(ll);
				}
				else
				{
					B_on_gran -= Vec_all.dot(ll);
				}
			}

			B_on_gran /= gr->area[0];


			M(ig, 0) = normal[0];
			M(ig, 1) = normal[1];
			M(ig, 2) = normal[2];
			F[ig] = B_on_gran;
		}

		Eigen::Vector3d vvv = M.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(F);

		cell->parameters[0]["rotB/b2_x"] = vvv[0];
		cell->parameters[0]["rotB/b2_y"] = vvv[1];
		cell->parameters[0]["rotB/b2_z"] = vvv[2];

		/*cout << F[0] << " " <<
			F[1] << " " <<
			F[2] << " " <<
			F[3] << " " <<
			F[4] << " " <<
			F[5] << endl;

		cout << M(0, 0) << " " << M(0, 1) << " " << M(0, 2) << endl;
		cout << M(1, 0) << " " << M(1, 1) << " " << M(1, 2) << endl;
		cout << M(2, 0) << " " << M(2, 1) << " " << M(2, 2) << endl;
		cout << M(3, 0) << " " << M(3, 1) << " " << M(3, 2) << endl;
		cout << M(4, 0) << " " << M(4, 1) << " " << M(4, 2) << endl;
		cout << M(5, 0) << " " << M(5, 1) << " " << M(5, 2) << endl;

		cout << vvv[0] << " " << vvv[1] << " " << vvv[2] << endl;

		exit(-1);*/
	}

	cout << "End: Culc_rotor_in_cell" << endl;
}


void Setka::Culc_divergence_in_cell(void)
{
	cout << "Start: Culc_divergence_in_cell" << endl;

	this->phys_param->param_names.push_back("divV");
	this->phys_param->param_names.push_back("div_eB");
	this->phys_param->param_names.push_back("div_et");
	this->phys_param->param_names.push_back("gradK_x");
	this->phys_param->param_names.push_back("gradK_y");
	this->phys_param->param_names.push_back("gradK_z");
	this->phys_param->param_names.push_back("et_x");
	this->phys_param->param_names.push_back("et_y");
	this->phys_param->param_names.push_back("et_z");
	// Добавили переменную для интерполяции

#pragma omp parallel for
	for (size_t i_step = 0; i_step < this->All_Cell.size(); i_step++)
	{
		unordered_map<string, double> par_left, par_right;
		Eigen::Vector3d normal, Vel;
		Eigen::Vector3d eB, gradK;

		auto& cell = this->All_Cell[i_step];
		double divV = 0.0;
		double div_eB = 0.0;
		double div_et = 0.0;
		gradK << 0.0, 0.0, 0.0;

		// пробегаемся по всем граням 
		for (const auto& gr : cell->grans)
		{
			eB << 0.0, 0.0, 0.0;
			Vel << 0.0, 0.0, 0.0;


			if (gr->cells[0] == cell)
			{
				normal << gr->normal[0][0], gr->normal[0][1], gr->normal[0][2];
			}
			else
			{
				normal << -gr->normal[0][0], -gr->normal[0][1], -gr->normal[0][2];
			}

			if (gr->type == Type_Gran::Us)
			{
				this->Snos_on_Gran(gr, par_left, par_right, 0);

				if (gr->type2 == Type_Gran_surf::Us)
				{
					Vel << (par_left["Vx"] + par_right["Vx"])/2.0, 
						(par_left["Vy"] + par_right["Vy"]) / 2.0,
						(par_left["Vz"] + par_right["Vz"]) / 2.0;
					eB << (par_left["Bx"] + par_right["Bx"]) / 2.0,
						(par_left["By"] + par_right["By"]) / 2.0,
						(par_left["Bz"] + par_right["Bz"]) / 2.0;
				}
				else
				{
					if (gr->cells[0] == cell)
					{
						Vel << par_left["Vx"],
							par_left["Vy"],
							par_left["Vz"];
						eB << par_left["Bx"],
							par_left["By"],
							par_left["Bz"];
					}
					else
					{
						Vel << par_right["Vx"],
							par_right["Vy"],
							par_right["Vz"];
						eB << par_right["Bx"],
							par_right["By"],
							par_right["Bz"];
					}
				}
			}
			else
			{
				Vel << cell->parameters[0]["Vx"], cell->parameters[0]["Vy"], cell->parameters[0]["Vz"];
				eB << cell->parameters[0]["Bx"], cell->parameters[0]["By"], cell->parameters[0]["Bz"];
			}

			double norm = eB.norm();
			eB /= norm;

			Eigen::Vector3d t, m;
			get_bazis(eB, t, m);

			gradK += normal * gr->area[0]/ norm;
			
			divV += Vel.dot(normal) * gr->area[0];
			div_eB += eB.dot(normal) * gr->area[0];
			div_et += t.dot(normal) * gr->area[0];
		}

		cell->parameters[0]["divV"] = divV/cell->volume[0];
		cell->parameters[0]["div_eB"] = div_eB /cell->volume[0];
		gradK = gradK /cell->volume[0];

		cell->parameters[0]["div_et"] = div_et / cell->volume[0];

		cell->parameters[0]["gradK_x"] = gradK[0];
		cell->parameters[0]["gradK_y"] = gradK[1];
		cell->parameters[0]["gradK_z"] = gradK[2];

		Eigen::Vector3d t, m;
		eB << cell->parameters[0]["Bx"], cell->parameters[0]["By"],
			cell->parameters[0]["Bz"];
		eB.normalize();
		get_bazis(eB, t, m);

		cell->parameters[0]["et_x"] = t[0];
		cell->parameters[0]["et_y"] = t[1];
		cell->parameters[0]["et_z"] = t[2];
	}

	cout << "End: Culc_divergence_in_cell" << endl;
}
