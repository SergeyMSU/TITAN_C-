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
		else if(i->center[0][0] > 0.0)
		//else if(i->center[0][0] > this->geo->L7 + 0.0001)
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
			if ( sqrt(kv(i->center[0][1]) + kv(i->center[0][2])) > 300)
			{
				i->parameters[0]["rho"] = 1.0;
				i->parameters[0]["p"] = 1.0;
				i->parameters[0]["Vx"] = this->phys_param->Velosity_inf;
				i->parameters[0]["Vy"] = 0.0;
				i->parameters[0]["Vz"] = 0.0;
				i->parameters[0]["Bx"] = -this->phys_param->B_inf * cos(this->phys_param->alphaB_inf);
				i->parameters[0]["By"] = -this->phys_param->B_inf * sin(this->phys_param->alphaB_inf);
				i->parameters[0]["Bz"] = 0.0;
				i->parameters[0]["Q"] = 100.0;

				/*i->parameters[0]["rho_H4"] = 1.0;
				i->parameters[0]["Vx_H4"] = this->phys_param->Velosity_inf;
				i->parameters[0]["Vy_H4"] = 0.0;
				i->parameters[0]["Vz_H4"] = 0.0;
				i->parameters[0]["p_H4"] = 0.5;*/

				i->parameters[1] = i->parameters[0];
			}
		}
	}

	// Если ввели какие-то новые переменные, их надо заполнить
	if (false)
	{
		for (auto& i : this->All_Cell)
		{
			i->parameters[0]["rho_H1"] = 0.0001;
			i->parameters[0]["Vx_H1"] = 0.0;
			i->parameters[0]["Vy_H1"] = 0.0;
			i->parameters[0]["Vz_H1"] = 0.0;
			i->parameters[0]["p_H1"] = 0.0001;

			i->parameters[0]["rho_H2"] = 0.0001;
			i->parameters[0]["Vx_H2"] = 0.0;
			i->parameters[0]["Vy_H2"] = 0.0;
			i->parameters[0]["Vz_H2"] = 0.0;
			i->parameters[0]["p_H2"] = 0.0001;

			i->parameters[0]["rho_H3"] = 0.0001;
			i->parameters[0]["Vx_H3"] = 0.0;
			i->parameters[0]["Vy_H3"] = 0.0;
			i->parameters[0]["Vz_H3"] = 0.0;
			i->parameters[0]["p_H3"] = 0.0001;

			i->parameters[0]["rho_H4"] = 1.0;
			i->parameters[0]["Vx_H4"] = this->phys_param->Velosity_inf;
			i->parameters[0]["Vy_H4"] = 0.0;
			i->parameters[0]["Vz_H4"] = 0.0;
			i->parameters[0]["p_H4"] = 0.5;

			for (short unsigned int j = 1; j < i->parameters.size(); j++)
			{
				i->parameters[j] = i->parameters[0];
			}
		}
	}

	bool tt1 = false;
	// Проверяем наличие всех необходимых переменных (инициализируем их, если их нет)
	for (auto& i : this->All_Cell)
	{
		for (auto& num : this->phys_param->param_names)
		{
			if (i->parameters[0].find(num) == i->parameters[0].end())
			{
				cout << "Parameters: " << num << "    ne opredelen v cells" << endl;
				tt1 = true;
			}
		}
		if (tt1 == true)
		{
			cout << "Error 4532896514" << endl;
			exit(-1);
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
				i->parameters[0]["Q"] = 100.0;
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
				i->parameters["rho"] = 1.0;
				i->parameters["p"] = 1.0;
				i->parameters["Vx"] = this->phys_param->Velosity_inf;
				i->parameters["Vy"] = 0.0;
				i->parameters["Vz"] = 0.0;
				i->parameters["Bx"] = -this->phys_param->B_inf * cos(this->phys_param->alphaB_inf);
				i->parameters["By"] = -this->phys_param->B_inf * sin(this->phys_param->alphaB_inf);
				i->parameters["Bz"] = 0.0;
				i->parameters["Q"] = 100.0;

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

				i->parameters["rho"] = this->phys_param->Get_rho_0(the) * pow(this->phys_param->R_0 / r, 2);
				i->parameters["p"] = this->phys_param->p_0 * pow(this->phys_param->R_0 / r, 2 * this->phys_param->gamma);
				i->parameters["Vx"] = mV * vec(0);
				i->parameters["Vy"] = mV * vec(1);
				i->parameters["Vz"] = mV * vec(2);
				i->parameters["Bx"] = cc(0);
				i->parameters["By"] = cc(1);
				i->parameters["Bz"] = cc(2);
				i->parameters["Q"] = i->parameters["rho"];

				i->parameters["rho_H1"] = 0.0001;
				i->parameters["Vx_H1"] = mV * vec(0);
				i->parameters["Vy_H1"] = mV * vec(1);
				i->parameters["Vz_H1"] = mV * vec(2);
				i->parameters["p_H1"] = 0.0001;

				i->parameters["rho_H2"] = 0.0001;
				i->parameters["Vx_H2"] = mV * vec(0);
				i->parameters["Vy_H2"] = mV * vec(1);
				i->parameters["Vz_H2"] = mV * vec(2);
				i->parameters["p_H2"] = 0.0001;

				i->parameters["rho_H3"] = 0.0001;
				i->parameters["Vx_H3"] = mV * vec(0);
				i->parameters["Vy_H3"] = mV * vec(1);
				i->parameters["Vz_H3"] = mV * vec(2);
				i->parameters["p_H3"] = 0.0001;
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
				this->Cell_Center->parameters[0][num] += i->parameters[num];
			}
		}

		for (auto& num : this->phys_param->param_names)
		{
			this->Cell_Center->parameters[0][num] /= nk;
		}
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
	kk[zone - 1] = 1.0;



	double S1 = 0.0;
	double S2 = 0.0;


	int i = 0;
	for (auto& nam : this->phys_param->H_name)
	{
		U_M_H[i] = sqrt(kv(C->parameters[now]["Vx"] - C->parameters[now]["Vx" + nam])
			+ kv(C->parameters[now]["Vy"] - C->parameters[now]["Vy" + nam])
			+ kv(C->parameters[now]["Vz"] - C->parameters[now]["Vz" + nam])
			+ (64.0 / (9.0 * const_pi)) *
			(C->parameters[now]["p"] / C->parameters[now]["rho"]
				+ 2.0 * C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam]));

		U_H[i] = sqrt(kv(C->parameters[now]["Vx"] - C->parameters[now]["Vx" + nam])
			+ kv(C->parameters[now]["Vy"] - C->parameters[now]["Vy" + nam])
			+ kv(C->parameters[now]["Vz"] - C->parameters[now]["Vz" + nam])
			+ (4.0 / const_pi) *
			(C->parameters[now]["p"] / C->parameters[now]["rho"]
				+ 2.0 * C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam]));


		sigma[i] = kv(1.0 - this->phys_param->par_a_2 * log(U_M_H[i]));
		nu[i] = C->parameters[now]["rho"] *
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
			(2.0 * C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam] - C->parameters[now]["p"] / C->parameters[now]["rho"]));
		i++;
	}

	double ddp = (this->phys_param->par_n_p_LISM / this->phys_param->par_Kn);
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
				+ (U_H[i] / U_M_H[i]) * (C->parameters[now]["p"] / C->parameters[now]["rho"]));
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
		if (M > 1)
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



	this->Test_geometr();
	// Все настройки расчёта считываются из файла Setter.txt
	// Сначала реализовываем расчёт без движения сетки
	unsigned int steps = steps__;

	unsigned short int now1 = 1;
	unsigned short int now2 = 0;

	vector<Gran*>* gran_list;
	vector<Cell*>* cell_list;

	// Если хотим отдельно считать внутреннюю и наружную области
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


		// Зададим условие для 4-го сорта водорода вблизи звезды
		// И для остальных сортов мягкие на границе
		for (auto& gran : this->All_boundary_Gran)
		{
			// Пока заменил на фиктивную центральную ячейку, если она работает это можно убрать
			if (false)//(gran->type == Type_Gran::Inner_Hard)
			{
				double u1, v1, w1;
				u1 = gran->cells[0]->parameters[now1]["Vx_H4"];
				v1 = gran->cells[0]->parameters[now1]["Vy_H4"];
				w1 = gran->cells[0]->parameters[now1]["Vz_H4"];
				if (u1 * gran->normal[now1][0] + u1 * gran->normal[now1][1] +
					u1 * gran->normal[now1][2] > 0.0)
				{
					gran->parameters["rho_H4"] = gran->cells[0]->parameters[now1]["rho_H4"];
					gran->parameters["p_H4"] = gran->cells[0]->parameters[now1]["p_H4"];
					gran->parameters["Vx_H4"] = gran->cells[0]->parameters[now1]["Vx_H4"];
					gran->parameters["Vy_H4"] = gran->cells[0]->parameters[now1]["Vy_H4"];
					gran->parameters["Vz_H4"] = gran->cells[0]->parameters[now1]["Vz_H4"];
				}
				else
				{
					gran->parameters["rho_H4"] = 0.0001;
					gran->parameters["p_H4"] = 0.0001;
					gran->parameters["Vx_H4"] = gran->cells[0]->parameters[now1]["Vx"];
					gran->parameters["Vy_H4"] = gran->cells[0]->parameters[now1]["Vy"];
					gran->parameters["Vz_H4"] = gran->cells[0]->parameters[now1]["Vz"];
				}

				/*if (gran->parameters["rho_H4"] < 0.0001)
				{
					gran->parameters["rho_H4"] = 0.0001;
					gran->parameters["p_H4"] = 0.00005;
				}*/
			}
			else if (gran->type == Type_Gran::Outer_Hard)
			{
				double u1, v1, w1;
				for (auto& nam : this->phys_param->H_name)
				{
					if (nam == "_H4") continue;

					u1 = gran->cells[0]->parameters[now1]["Vx" + nam];
					v1 = gran->cells[0]->parameters[now1]["Vy" + nam];
					w1 = gran->cells[0]->parameters[now1]["Vz" + nam];
					if (u1 * gran->normal[now1][0] + u1 * gran->normal[now1][1] +
						u1 * gran->normal[now1][2] > 0.0)
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
						gran->parameters["Vx" + nam] = 0.0001 * gran->normal[now1][0];
						gran->parameters["Vy" + nam] = 0.0001 * gran->normal[now1][1];
						gran->parameters["Vz" + nam] = 0.0001 * gran->normal[now1][2];
					}
				}
			}
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
			double ntnt = this->Culc_Gran_Potok(gran, now1, metod, nmnm);

			if (ntnt < loc_time)
			{
				loc_time = min(loc_time, ntnt);  // Считает потоки через данную грань (записывает результат в параметры грани)
				xc_min = gran->center[now1][0];
				yc_min = gran->center[now1][1];
				zc_min = gran->center[now1][2];
				name_min_time = nmnm;
			}
		}

#pragma omp barrier
		//cout << "barrier" << endl;

		// Расчитываем законы сохранения в ячейках
#pragma omp parallel for
		for (size_t i_step = 0; i_step < cell_list->size(); i_step++)
		{
			auto& cell = (*cell_list)[i_step];
			double Volume = cell->volume[now1];
			double Volume2 = cell->volume[now2];
			int zone;


			if (Volume != Volume2)
			{
				cout << "Error 8767654534" << endl;
				exit(-1);
			}

			std::vector<double> POTOK;
			POTOK.resize(9);

			std::vector<double> POTOK2; // Для конвективных параметров
			POTOK2.resize(1);

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

			double rho3, u3, v3, w3, bx3, by3, bz3, p3, Q3;
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

				if (rho3 < 0.0000000001)
				{
					rho3 = 0.001;
					Q3 = Q / rho * rho3;
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
					whach(Volume);
					whach(Volume2);
					whach(rho);
					whach(time);
					whach(Q);

					for (short unsigned int ik = 0; ik < 9; ik++)
					{
						cout << "ik: " << POTOK[ik] << endl;
					}
					exit(-1);
				}


				if (p3 < 0.0000000001)
				{
					p3 = 0.0000001;
				}

				cell->parameters[now2]["rho"] = rho3;
				if (cell->parameters[now1].find("Q") != cell->parameters[now1].end())
				{
					cell->parameters[now2]["Q"] = Q3;
				}
				cell->parameters[now2]["Vx"] = u3;
				cell->parameters[now2]["Vy"] = v3;
				cell->parameters[now2]["Vz"] = w3;
				cell->parameters[now2]["Bx"] = bx3;
				cell->parameters[now2]["By"] = by3;
				cell->parameters[now2]["Bz"] = bz3;
				cell->parameters[now2]["p"] = p3;
			}

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

				if (rho3 < 0.0000000001)
				{
					rho3 = 0.01;
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

				if (p3 < 0.0000000001)
				{
					p3 = 0.005;
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

					if (rho3 < 0.0000000001)
					{
						rho3 = 0.00001;
					}

					u3 = (rho * vx - time * (POTOK_F[i][1]) / Volume) / rho3;
					v3 = (rho * vy - time * (POTOK_F[i][2]) / Volume) / rho3;
					w3 = (rho * vz - time * (POTOK_F[i][3]) / Volume) / rho3;


					p3 = (((p / this->phys_param->g1 + 0.5 * rho * kvv(vx, vy, vz))
						- time * (POTOK_F[i][4]) / Volume) -
						0.5 * rho3 * kvv(u3, v3, w3)) * this->phys_param->g1;

					if (p3 < 0.0000000001)
					{
						p3 = 0.00001;
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

double Setka::Culc_Gran_Potok(Gran* gr, unsigned short int now, short int metod, string& name)
{
	double dsr, dsc, dsl;
	std::vector<double> qqq, qqq1, qqq2;
	qqq.resize(8);
	qqq1.resize(8);
	qqq2.resize(8);
	std::vector<double> konvect_left, konvect_right, konvect;
	PrintOptions Option = PrintOptions{};
	double area = gr->area[now];


	if (gr->type == Type_Gran::Us) // Обычная грань
	{
		auto A = gr->cells[0];
		auto B = gr->cells[1];

		double loc_time = 10000000.0;
		double dist = norm2(A->center[now][0] - B->center[now][0],
			A->center[now][1] - B->center[now][1],
			A->center[now][2] - B->center[now][2]) / 2.0;

		if (this->phys_param->culc_plasma == true)
		{
			// Без TVD
			if (this->phys_param->TVD == false)
			{
				if (regim_otladki == true)
				{
					if (A->parameters[now].find("rho") == A->parameters[now].end() ||
						B->parameters[now].find("rho") == B->parameters[now].end())
					{
						cout << "Error  0956453978" << endl;
						exit(-1);
					}


					if (std::isnan(A->parameters[now]["rho"]) || std::fpclassify(A->parameters[now]["rho"]) == FP_SUBNORMAL)
					{
						cout << "Error 0907675453" << endl;
						for (const auto& [key, value] : A->parameters[now]) {
							std::cout << key << ": " << value << std::endl;
						}
						exit(-1);
					}
				}

				qqq1[0] = A->parameters[now]["rho"];
				qqq1[1] = A->parameters[now]["Vx"];
				qqq1[2] = A->parameters[now]["Vy"];
				qqq1[3] = A->parameters[now]["Vz"];
				qqq1[4] = A->parameters[now]["p"];
				qqq1[5] = A->parameters[now]["Bx"];
				qqq1[6] = A->parameters[now]["By"];
				qqq1[7] = A->parameters[now]["Bz"];

				qqq2[0] = B->parameters[now]["rho"];
				qqq2[1] = B->parameters[now]["Vx"];
				qqq2[2] = B->parameters[now]["Vy"];
				qqq2[3] = B->parameters[now]["Vz"];
				qqq2[4] = B->parameters[now]["p"];
				qqq2[5] = B->parameters[now]["Bx"];
				qqq2[6] = B->parameters[now]["By"];
				qqq2[7] = B->parameters[now]["Bz"];
			}
			else // TVD
			{

			}

			double w = 0.0;

			Option.x = gr->center[now][0];
			Option.y = gr->center[now][1];
			Option.z = gr->center[now][2];

			if(A->parameters[now].find("Q") != A->parameters[now].end())
			{
				konvect_left.push_back(A->parameters[now]["Q"]);
				konvect_right.push_back(B->parameters[now]["Q"]);
				konvect.push_back(0.0);
			}

			this->phys_param->chlld(metod, gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 3,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option);

			double dnt = this->phys_param->KFL * dist
				/ (max(fabs(dsl), fabs(dsr)) + fabs(w));

			if (dnt < loc_time)
			{
				loc_time = min(loc_time, dnt);
				name = "plasma";
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

			if (A->parameters[now].find("Q") != A->parameters[now].end())
			{
				gr->parameters["PQ"] = konvect[0] * area; // Может быть неправльный порядок при других конвективных переменных
			}
		}

		// Теперь считаем для водорода

		konvect_left.clear();
		konvect_right.clear(); 
		konvect.clear();

		for (auto& nam : this->phys_param->H_name)
		{
			if (this->phys_param->TVD == false)
			{
				qqq1[0] = A->parameters[now]["rho" + nam];
				qqq1[1] = A->parameters[now]["Vx" + nam];
				qqq1[2] = A->parameters[now]["Vy" + nam];
				qqq1[3] = A->parameters[now]["Vz" + nam];
				qqq1[4] = A->parameters[now]["p" + nam];
				qqq1[5] = 0.0;
				qqq1[6] = 0.0;
				qqq1[7] = 0.0;

				qqq2[0] = B->parameters[now]["rho" + nam];
				qqq2[1] = B->parameters[now]["Vx" + nam];
				qqq2[2] = B->parameters[now]["Vy" + nam];
				qqq2[3] = B->parameters[now]["Vz" + nam];
				qqq2[4] = B->parameters[now]["p" + nam];
				qqq2[5] = 0.0;
				qqq2[6] = 0.0;
				qqq2[7] = 0.0;
			}
			else // TVD
			{

			}

			double w = 0.0;

			this->phys_param->chlld(metod, gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 3,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option);

			double dnt = this->phys_param->KFL * dist
				/ (max(fabs(dsl), fabs(dsr)) + fabs(w));

			if (dnt < loc_time)
			{
				loc_time = min(loc_time, dnt);
				name = nam;
			}

			gr->parameters["Pm" + nam] = qqq[0] * area;
			gr->parameters["PVx" + nam] = qqq[1] * area;
			gr->parameters["PVy" + nam] = qqq[2] * area;
			gr->parameters["PVz" + nam] = qqq[3] * area;
			gr->parameters["Pe" + nam] = qqq[4] * area;

		}


		return loc_time;
		
	}
	else if (gr->type == Type_Gran::Inner_Hard || gr->type == Type_Gran::Outer_Hard)
	{
		if (this->phys_param->culc_plasma == true)
		{
			this->phys_param->Get_Potok(gr->parameters["rho"], gr->parameters["p"],
				gr->parameters["Vx"], gr->parameters["Vy"], gr->parameters["Vz"],
				gr->parameters["Bx"], gr->parameters["By"], gr->parameters["Bz"],
				gr->normal[now][0], gr->normal[now][1], gr->normal[now][2],
				qqq);

			gr->parameters["Pm"] = qqq[0] * area;
			gr->parameters["PVx"] = qqq[1] * area;
			gr->parameters["PVy"] = qqq[2] * area;
			gr->parameters["PVz"] = qqq[3] * area;
			gr->parameters["Pe"] = qqq[7] * area;
			gr->parameters["PBx"] = qqq[4] * area;
			gr->parameters["PBy"] = qqq[5] * area;
			gr->parameters["PBz"] = qqq[6] * area;
			if (gr->parameters.find("Q") != gr->parameters.end())
			{
				gr->parameters["PQ"] = gr->parameters["Q"] * qqq[0]/ gr->parameters["rho"] * area;
			}

			gr->parameters["PdivB"] = scalarProductFast(gr->normal[now][0],
				gr->normal[now][1], gr->normal[now][2],
				gr->parameters["Bx"], gr->parameters["By"],
				gr->parameters["Bz"]) * area;
		}

		for (auto& nam : this->phys_param->H_name)
		{
			if (gr->type == Type_Gran::Outer_Hard || (nam == "_H1") )
			{
				this->phys_param->Get_Potok(gr->parameters["rho" + nam], gr->parameters["p" + nam],
					gr->parameters["Vx" + nam], gr->parameters["Vy" + nam], gr->parameters["Vz" + nam],
					0.0, 0.0, 0.0,
					gr->normal[now][0], gr->normal[now][1], gr->normal[now][2],
					qqq);

				gr->parameters["Pm" + nam] = qqq[0] * area;
				gr->parameters["PVx" + nam] = qqq[1] * area;
				gr->parameters["PVy" + nam] = qqq[2] * area;
				gr->parameters["PVz" + nam] = qqq[3] * area;
				gr->parameters["Pe" + nam] = qqq[7] * area;
			}
			else // Для  H2  H3  H4    и  Inner_Hard  -  делаем через центральную фиктувную ячейку
			{
				auto A = gr->cells[0];
				auto B = this->Cell_Center;

				if (this->phys_param->TVD == false)
				{
					qqq1[0] = A->parameters[now]["rho" + nam];
					qqq1[1] = A->parameters[now]["Vx" + nam];
					qqq1[2] = A->parameters[now]["Vy" + nam];
					qqq1[3] = A->parameters[now]["Vz" + nam];
					qqq1[4] = A->parameters[now]["p" + nam];
					qqq1[5] = 0.0;
					qqq1[6] = 0.0;
					qqq1[7] = 0.0;

					qqq2[0] = B->parameters[now]["rho" + nam];
					qqq2[1] = B->parameters[now]["Vx" + nam];
					qqq2[2] = B->parameters[now]["Vy" + nam];
					qqq2[3] = B->parameters[now]["Vz" + nam];
					qqq2[4] = B->parameters[now]["p" + nam];
					qqq2[5] = 0.0;
					qqq2[6] = 0.0;
					qqq2[7] = 0.0;
				}
				else // TVD
				{

				}

				double w = 0.0;

				this->phys_param->chlld(metod, gr->normal[now][0], gr->normal[now][1],
					gr->normal[now][2],
					w, qqq1, qqq2, qqq, false, 3,
					konvect_left, konvect_right, konvect, dsr, dsc, dsl,
					Option);

				gr->parameters["Pm" + nam] = qqq[0] * area;
				gr->parameters["PVx" + nam] = qqq[1] * area;
				gr->parameters["PVy" + nam] = qqq[2] * area;
				gr->parameters["PVz" + nam] = qqq[3] * area;
				gr->parameters["Pe" + nam] = qqq[4] * area;
			}
		}

		return 1000000.0;   // возвращаем большой шаг по времени
	}
	else if (gr->type == Type_Gran::Outer_Soft)
	{
		auto C = gr->cells[0];
		double loc_time = 10000000.0;
		double dist = norm2(C->center[now][0] - gr->center[now][0],
			C->center[now][1] - gr->center[now][1],
			C->center[now][2] - gr->center[now][2]);

		if (this->phys_param->culc_plasma == true)
		{
			qqq1[0] = C->parameters[now]["rho"];
			qqq1[1] = C->parameters[now]["Vx"];
			qqq1[2] = C->parameters[now]["Vy"];
			qqq1[3] = C->parameters[now]["Vz"];
			qqq1[4] = C->parameters[now]["p"];
			qqq1[5] = C->parameters[now]["Bx"];
			qqq1[6] = C->parameters[now]["By"];
			qqq1[7] = C->parameters[now]["Bz"];

			qqq2[0] = qqq1[0];
			qqq2[1] = qqq1[1];
			qqq2[2] = qqq1[2];
			qqq2[3] = qqq1[3];
			qqq2[4] = qqq1[4];
			qqq2[5] = qqq1[5];
			qqq2[6] = qqq1[6];
			qqq2[7] = qqq1[7];

			if (C->parameters[now].find("Q") != C->parameters[now].end())
			{
				konvect_left.push_back(C->parameters[now]["Q"]);
				konvect_right.push_back(C->parameters[now]["Q"]);
				konvect.push_back(0.0);
			}

			// Запрещаем затекание жидкости через мягкие граничные условия
			if (qqq2[1] * gr->normal[now][0] + qqq2[2] * gr->normal[now][1] +
				qqq2[3] * gr->normal[now][2] < 0.0)
			{
				qqq2[1] = 0.0;
				qqq2[2] = 0.0;
				qqq2[3] = 0.0;
			}

			if (gr->normal[now][0] < -0.9) // Для задней границы
			{
				// Отсос
				if (qqq2[1] > this->phys_param->Velosity_inf / 7.0)
				{
					qqq2[1] = this->phys_param->Velosity_inf / 5.0;
				}
			}

			double w = 0.0;
			Option.x = gr->center[now][0];
			Option.y = gr->center[now][1];
			Option.z = gr->center[now][2];

			this->phys_param->chlld(metod, gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 3,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option);

			loc_time = min(loc_time, this->phys_param->KFL * dist / (max(fabs(dsl), fabs(dsr)) + fabs(w)));

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

			if (C->parameters[now].find("Q") != C->parameters[now].end())
			{
				gr->parameters["PQ"] = konvect[0] * area; // Может быть неправльный порядок при других конвективных переменных
			}

		}
			

		if (false)//(this->phys_param->culc_plasma == true)
		{
			gr->parameters["rho"] = C->parameters[now]["rho"];
			gr->parameters["p"] = C->parameters[now]["p"];
			gr->parameters["Vx"] = C->parameters[now]["Vx"];
			gr->parameters["Vy"] = C->parameters[now]["Vy"];
			gr->parameters["Vz"] = C->parameters[now]["Vz"];
			gr->parameters["Bx"] = C->parameters[now]["Bx"];
			gr->parameters["By"] = C->parameters[now]["By"];
			gr->parameters["Bz"] = C->parameters[now]["Bz"];

			if (C->parameters[now].find("Q") != C->parameters[now].end())
			{
				gr->parameters["Q"] = C->parameters[now]["Q"];
			}

			// Запрещаем затекание жидкости через мягкие граничные условия
			if (gr->parameters["Vx"] * gr->normal[now][0] + gr->parameters["Vy"] * gr->normal[now][1] +
				gr->parameters["Vz"] * gr->normal[now][2] < 0.0)
			{
				gr->parameters["Vx"] = 0.0;
				gr->parameters["Vy"] = 0.0;
				gr->parameters["Vz"] = 0.0;
			}

			if (gr->normal[now][0] < -0.9) // Для задней границы
			{
				// Отсос
				if (gr->parameters["Vx"] > this->phys_param->Velosity_inf / 7.0)
				{
					gr->parameters["Vx"] = this->phys_param->Velosity_inf / 5.0;
				}

				// Сильно большой отсос тоже плохо
				if (gr->parameters["Vx"] < 5.0 * this->phys_param->Velosity_inf)
				{
					gr->parameters["Vx"] = this->phys_param->Velosity_inf;
				}
			}

			this->phys_param->Get_Potok(gr->parameters["rho"], gr->parameters["p"],
				gr->parameters["Vx"], gr->parameters["Vy"], gr->parameters["Vz"],
				gr->parameters["Bx"], gr->parameters["By"], gr->parameters["Bz"],
				gr->normal[now][0], gr->normal[now][1], gr->normal[now][2],
				qqq);

			gr->parameters["Pm"] = qqq[0] * area;
			gr->parameters["PVx"] = qqq[1] * area;
			gr->parameters["PVy"] = qqq[2] * area;
			gr->parameters["PVz"] = qqq[3] * area;
			gr->parameters["Pe"] = qqq[7] * area;
			gr->parameters["PBx"] = qqq[4] * area;
			gr->parameters["PBy"] = qqq[5] * area;
			gr->parameters["PBz"] = qqq[6] * area;
			gr->parameters["PdivB"] = scalarProductFast(gr->normal[now][0],
				gr->normal[now][1], gr->normal[now][2],
				gr->parameters["Bx"], gr->parameters["By"],
				gr->parameters["Bz"]) * area;
		}

		for (auto& nam : this->phys_param->H_name)
		{
			gr->parameters["rho" + nam] = C->parameters[now]["rho" + nam];
			gr->parameters["p" + nam] = C->parameters[now]["p" + nam];
			gr->parameters["Vx" + nam] = C->parameters[now]["Vx" + nam];
			gr->parameters["Vy" + nam] = C->parameters[now]["Vy" + nam];
			gr->parameters["Vz" + nam] = C->parameters[now]["Vz" + nam];


			this->phys_param->Get_Potok(gr->parameters["rho" + nam], gr->parameters["p" + nam],
				gr->parameters["Vx" + nam], gr->parameters["Vy" + nam], gr->parameters["Vz" + nam],
				0.0, 0.0, 0.0,
				gr->normal[now][0], gr->normal[now][1], gr->normal[now][2],
				qqq);

			gr->parameters["Pm" + nam] = qqq[0] * area;
			gr->parameters["PVx" + nam] = qqq[1] * area;
			gr->parameters["PVy" + nam] = qqq[2] * area;
			gr->parameters["PVz" + nam] = qqq[3] * area;
			gr->parameters["Pe" + nam] = qqq[7] * area;

		}


		return loc_time;   // возвращаем большой шаг по времени
	}
	else
	{
		cout << "Error 6323145345" << endl;
		exit(-1);
	}
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

			ii->parameters[0][key] = value;
			ii->parameters[1][key] = value;
		}
	}
}
