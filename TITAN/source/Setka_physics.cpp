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

	this->All_boundary_Gran.clear();
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


	// Далее задаём узлам параметр - расстояние до HP (для особых схем вдоль HP)
	for (auto& i : this->All_Luch)
	{
		if (i->type == "A_Luch" || i->type == "A2_Luch" || i->type == "B_Luch")
		{
			auto a1 = i->get_yzel_near_opor(2, -1);
			a1->dist_from_HP = 1;

			a1 = i->get_yzel_near_opor(2, 1);
			a1->dist_from_HP = 1;

			a1 = i->get_yzel_near_opor(2, -2);
			a1->dist_from_HP = 2;

			a1 = i->get_yzel_near_opor(2, 2);
			a1->dist_from_HP = 2;

			a1 = i->get_yzel_near_opor(2, -3);
			a1->dist_from_HP = 3;

			a1 = i->get_yzel_near_opor(2, 3);
			a1->dist_from_HP = 3;

			a1 = i->get_yzel_near_opor(2, -4);
			a1->dist_from_HP = 4;

			a1 = i->get_yzel_near_opor(2, 4);
			a1->dist_from_HP = 4;
		}
		else if (i->type == "E_Luch")
		{
			auto a1 = i->Yzels[1];
			a1->dist_from_HP = 1;

			a1 = i->Yzels[2];
			a1->dist_from_HP = 2;

			a1 = i->Yzels[3];
			a1->dist_from_HP = 3;

			a1 = i->Yzels[4];
			a1->dist_from_HP = 4;
		}
		else if (i->type == "G_Luch")
		{
			auto a1 = i->get_yzel_near_opor(2, -1);
			a1->dist_from_HP = 1;

			a1 = i->get_yzel_near_opor(2, -2);
			a1->dist_from_HP = 2;

			a1 = i->get_yzel_near_opor(2, -3);
			a1->dist_from_HP = 3;

			a1 = i->get_yzel_near_opor(2, -4);
			a1->dist_from_HP = 4;
		}
		else if (i->type == "D_Luch")
		{
			auto a1 = i->get_yzel_near_opor(1, -1);
			if (a1->coord[0][0] < this->geo->L6) continue;
			a1->dist_from_HP = 1;

			a1 = i->get_yzel_near_opor(1, 1);
			a1->dist_from_HP = 1;

			a1 = i->get_yzel_near_opor(1, -2);
			a1->dist_from_HP = 2;

			a1 = i->get_yzel_near_opor(1, 2);
			a1->dist_from_HP = 2;

			a1 = i->get_yzel_near_opor(1, -3);
			a1->dist_from_HP = 3;

			a1 = i->get_yzel_near_opor(1, 3);
			a1->dist_from_HP = 3;

			a1 = i->get_yzel_near_opor(1, -4);
			a1->dist_from_HP = 4;

			a1 = i->get_yzel_near_opor(1, 4);
			a1->dist_from_HP = 4;
		}
	}
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
			
			double r = norm2(0.0, i->center[0][1], i->center[0][2]);

			int zone = determ_zone(i, 0);
			if (r > 330.0)
			{
				i->parameters[0]["rho"] = this->phys_param->rho_LISM;
				i->parameters[0]["rho_He"] = this->phys_param->rho_HE_LISM; 
				i->parameters[0]["p"] = this->phys_param->rho_p_LISM; 
				i->parameters[0]["Vx"] = this->phys_param->Velosity_inf;
				i->parameters[0]["Vy"] = 0.0;
				i->parameters[0]["Vz"] = 0.0;
				i->parameters[0]["Bx"] = this->phys_param->B_inf * cos(this->phys_param->alphaB_inf);
				i->parameters[0]["By"] = this->phys_param->B_inf * sin(this->phys_param->alphaB_inf);
				i->parameters[0]["Bz"] = 0.0;
				i->parameters[0]["Q"] = 100.0 * i->parameters[0]["rho"];
			}
			if (false)
			{
				i->parameters[0]["rho"] *= 0.81;
				i->parameters[0]["p"] *= 0.934;
				i->parameters[0]["rho_He"] *= 0.5;
				i->parameters[0]["Bx"] *= 1.56751;
				i->parameters[0]["By"] *= 1.56751;
				i->parameters[0]["Bz"] *= 1.56751;


				double angle = 20.0 * const_pi / 180.0;
				double cos_angle = cos(angle);
				double sin_angle = sin(angle);

				double Bx_old = i->parameters[0]["Bx"];
				double By_old = i->parameters[0]["By"];
				double Bz_old = i->parameters[0]["Bz"]; 

				
				// |cos(?)  -sin(?)   0| |Bx|
				// |sin(?)   cos(?)   0| |By|
				// |  0        0      1| |Bz|

				i->parameters[0]["Bx"] = cos_angle * Bx_old - sin_angle * By_old;
				i->parameters[0]["By"] = sin_angle * Bx_old + cos_angle * By_old;
				i->parameters[0]["Bz"] = Bz_old; 
			}

			i->parameters[1] = i->parameters[0];

		}
	}

	// Интерполяция
	if (false)
	{
		Interpol SS = Interpol("For_intertpolate_1.bin");
		std::unordered_map<string, double> parameters;
		bool fine_int;

		for (auto& i : this->All_Cell)
		{
			if (norm2(i->center[0][0], i->center[0][1], i->center[0][2]) < 10.0)
			{
				fine_int = SS.Get_param(i->center[0][0], i->center[0][1], i->center[0][2], 
					parameters);

				i->parameters[0]["rho_H1"] = parameters["rho_H1"];
				i->parameters[0]["Vx_H1"] = parameters["Vx_H1"];
				i->parameters[0]["Vy_H1"] = parameters["Vy_H1"];
				i->parameters[0]["Vz_H1"] = parameters["Vz_H1"];
				i->parameters[0]["p_H1"] = parameters["p_H1"];


				i->parameters[1] = i->parameters[0];
			}
		}
	}

	// Если ввели какие-то новые переменные, их надо заполнить
	if (false)
	{
		for (auto& i : this->All_Cell)
		{
			i->parameters[0]["rho_Pui_1"] = 1e-8;
			i->parameters[0]["p_Pui_1"] = 1e-8 / 2.0;

			i->parameters[0]["rho_Pui_2"] = 1e-8;
			i->parameters[0]["p_Pui_2"] = 1e-8 / 2.0;

			for (short unsigned int j = 1; j < i->parameters.size(); j++)
			{
				i->parameters[j] = i->parameters[0];
			}
		}
	}


	std::unordered_set<std::string> no_names;
	// Проверяем наличие всех необходимых переменных
	for (auto& i : this->All_Cell)
	{
		for (auto& num : this->phys_param->param_names)
		{
			if (i->parameters[0].find(num) == i->parameters[0].end())
			{
				i->parameters[0][num] = 1e-7;

				if (no_names.find(num) == no_names.end()) no_names.insert(num);
			}
		}

		i->parameters[1] = i->parameters[0];
	}

	if (!no_names.empty()) 
	{
		std::cout << "Ne bilo nekotorix peremennix v yacheikax, oni bili opredeleni nulem. Elements:" << std::endl;
		// Вариант 1: Через range-based for (C++11)
		for (const auto& str : no_names) {
			std::cout << " - " << str << std::endl;
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

			if (r < 2.0)
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

				mV = this->phys_param->Get_v_0(the / const_pi * 180.0);

				double Tp = this->phys_param->Get_T_0(the / const_pi * 180.0); // Температура

				double np = this->phys_param->Get_rho_0(the / const_pi * 180.0);

				double rho = (this->phys_param->mep + 2.0 * this->phys_param->mrho_He_0 *
					this->phys_param->mep +
					1.0 + 4.0 * this->phys_param->mrho_He_0) * np;



				/*i->parameters[0]["rho"] = rho * pow(this->phys_param->R_0 / r, 2);
				i->parameters[0]["rho_He"] = 4.0 * this->phys_param->mrho_He_0 * np * pow(this->phys_param->R_0 / r, 2);
				i->parameters[0]["p"] = (1.0 + 3.0 * this->phys_param->mrho_He_0 / 2.0) * np * Tp *
					pow(this->phys_param->R_0 / r, 2 * this->phys_param->gamma);
				i->parameters[0]["Vx"] = mV * vec(0)/r;
				i->parameters[0]["Vy"] = mV * vec(1)/r;
				i->parameters[0]["Vz"] = mV * vec(2)/r;
				i->parameters[0]["Bx"] = cc(0);
				i->parameters[0]["By"] = cc(1);
				i->parameters[0]["Bz"] = cc(2);
				i->parameters[0]["Q"] = i->parameters[0]["rho"];*/

				i->parameters[0]["rho_H1"] = 0.00001;
				i->parameters[0]["Vx_H1"] = mV * vec(0) / r;
				i->parameters[0]["Vy_H1"] = mV * vec(1) / r;
				i->parameters[0]["Vz_H1"] = mV * vec(2) / r;
				i->parameters[0]["p_H1"] = 0.00001;
			}
			else
			{
				/*i->parameters[0]["rho"] = 1.0;
				i->parameters[0]["p"] = 1.0;
				i->parameters[0]["Vx"] = this->phys_param->Velosity_inf;
				i->parameters[0]["Vy"] = 0.0;
				i->parameters[0]["Vz"] = 0.0;
				i->parameters[0]["Bx"] = -this->phys_param->B_inf * cos(this->phys_param->alphaB_inf);
				i->parameters[0]["By"] = -this->phys_param->B_inf * sin(this->phys_param->alphaB_inf);
				i->parameters[0]["Bz"] = 0.0;
				i->parameters[0]["Q"] = 100.0 * i->parameters[0]["rho"];*/
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
				i->parameters["rho"] = this->phys_param->rho_LISM; // 1.60063; // 1.0 * (1.0 + this->phys_param->mrho_He_inf);
				i->parameters["rho_He"] = this->phys_param->rho_HE_LISM; // 0.6; // this->phys_param->mrho_He_inf;
				i->parameters["p"] = this->phys_param->rho_p_LISM; // 1 + (i->parameters["rho_He"]) /
					//(i->parameters["rho"] - i->parameters["rho_He"]);
				i->parameters["Vx"] = this->phys_param->Velosity_inf;
				i->parameters["Vy"] = 0.0;
				i->parameters["Vz"] = 0.0;
				i->parameters["Bx"] = this->phys_param->B_inf * cos(this->phys_param->alphaB_inf);
				i->parameters["By"] = this->phys_param->B_inf * sin(this->phys_param->alphaB_inf);
				i->parameters["Bz"] = 0.0;
				i->parameters["Q"] = 100.0 * i->parameters["rho"];

				i->parameters["rho_H4"] = 1.0;
				i->parameters["Vx_H4"] = this->phys_param->Velosity_inf;
				i->parameters["Vy_H4"] = 0.0;
				i->parameters["Vz_H4"] = 0.0;
				i->parameters["p_H4"] = 0.5;

				if (this->phys_param->num_H >= 9)
				{
					i->parameters["rho_H9"] = 0.000001;
					i->parameters["Vx_H9"] = this->phys_param->Velosity_inf;
					i->parameters["Vy_H9"] = 0.0;
					i->parameters["Vz_H9"] = 0.0;
					i->parameters["p_H9"] = 0.000001;
				}

				if (this->phys_param->is_PUI == true)
				{
					for (const auto& nam : this->phys_param->pui_name)
					{
						i->parameters["rho" + nam] = 0.0;
						i->parameters["p" + nam] = 0.0;
					}
				}
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

				mV = this->phys_param->Get_v_0(the / const_pi * 180.0);

				double Tp = this->phys_param->Get_T_0(the / const_pi * 180.0); // Температура

				double np = this->phys_param->Get_rho_0(the / const_pi * 180.0);

				double rho = (this->phys_param->mep + 2.0 * this->phys_param->mrho_He_0 * 
					this->phys_param->mep + 
					1.0 + 4.0 * this->phys_param->mrho_He_0) * np;



				i->parameters["rho"] = rho * pow(this->phys_param->R_0 / r, 2);
				i->parameters["rho_He"] = 4.0 * this->phys_param->mrho_He_0 * np * pow(this->phys_param->R_0 / r, 2);
				i->parameters["p"] = (1.0 + 3.0 * this->phys_param->mrho_He_0 /2.0) * np * Tp *
					pow(this->phys_param->R_0 / r, 2 * this->phys_param->gamma);
				i->parameters["Vx"] = mV * vec(0)/r;
				i->parameters["Vy"] = mV * vec(1)/r;
				i->parameters["Vz"] = mV * vec(2)/r;
				i->parameters["Bx"] = cc(0);
				i->parameters["By"] = cc(1);
				i->parameters["Bz"] = cc(2);
				i->parameters["Q"] = i->parameters["rho"];

				i->parameters["rho_H1"] = 0.0001;
				i->parameters["Vx_H1"] = mV * vec(0) / r;
				i->parameters["Vy_H1"] = mV * vec(1) / r;
				i->parameters["Vz_H1"] = mV * vec(2) / r;
				i->parameters["p_H1"] = 0.0001;
				

				if (this->phys_param->num_H >= 5)
				{
					i->parameters["rho_H5"] = 0.0001;
					i->parameters["Vx_H5"] = mV * vec(0) / r;
					i->parameters["Vy_H5"] = mV * vec(1) / r;
					i->parameters["Vz_H5"] = mV * vec(2) / r;
					i->parameters["p_H5"] = 0.0001;
				}

				if (this->phys_param->is_PUI == true)
				{
					for (const auto& nam : this->phys_param->pui_name)
					{
						i->parameters["rho" + nam] = 0.0;
						i->parameters["p" + nam] = 0.0;
					}
				}
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

			i->cells[0]->is_TVD = false;

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

	// Для первых ячеек задаём магнитное поле
	if (true)
	{
		for (auto& i : this->All_Cell)
		{
			if (i->is_TVD == true) continue;

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

			cc = this->phys_param->Matr * vv;

			i->parameters[0]["Bx"] = cc(0);
			i->parameters[0]["By"] = cc(1);
			i->parameters[0]["Bz"] = cc(2);
			

			for (short unsigned int j = 1; j < i->parameters.size(); j++)
			{
				i->parameters[j] = i->parameters[0];
			}
		}
	}


	for (auto& i : this->All_Cell)
	{
		i->parameters[1] = i->parameters[0];
	}
}

void Setka::Init_physics_with_time(const double& time)
{
	this->Calculating_measure(0);
	this->Calculating_measure(1);
	double x, y, z, r, the;
	Eigen::Vector3d vec, cc, vv;
	double BR, BPHI, V1, V2, V3, mV;


	// Задаём граничные условия (на граничных гранях)
	if (true)
	{
		for (auto& i : this->All_boundary_Gran)
		{
			if (i->type == Type_Gran::Inner_Hard)
			{
				x = i->center[0][0];
				y = i->center[0][1];
				z = i->center[0][2];
				r = norm2(x, y, z);

				vec << x, y, z;

				cc = this->phys_param->Matr2 * vec;
				the = acos(cc(2) / r);

				the = -the + const_pi / 2.0;   // Т.к.в данных по СВ на 1 а.е.угол от - 90 до 90 у Алексашова

				mV = this->phys_param->Get_v_0(the / const_pi * 180.0);

				if (time < 0.0427075)  // 1 месяц
				{
					mV = mV * sqrt(1.3);
				}

				i->parameters["Vx"] = mV * vec(0) / r;
				i->parameters["Vy"] = mV * vec(1) / r;
				i->parameters["Vz"] = mV * vec(2) / r;
			}
		}
	}

	// Для первых ячеек задаём магнитное поле
	if (true)
	{
		for (auto& i : this->All_Cell)
		{
			if (i->is_TVD == true) continue;

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

			cc = this->phys_param->Matr * vv;

			i->parameters[0]["Bx"] = cc(0);
			i->parameters[0]["By"] = cc(1);
			i->parameters[0]["Bz"] = cc(2);


			for (short unsigned int j = 1; j < i->parameters.size(); j++)
			{
				i->parameters[j] = i->parameters[0];
			}
		}
	}

	for (auto& i : this->All_Cell)
	{
		i->parameters[1] = i->parameters[0];
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

double U_bera(const double& T1, const double& T2, const Eigen::Vector3d& V1, const Eigen::Vector3d& V2)
{
	double uu = (V1 - V2).norm();
	return sqrt((4.0 / const_pi) * (T1 + T2) + kv(uu));
}

double Um_bera(const double& T1, const double& T2, const Eigen::Vector3d& V1, const Eigen::Vector3d& V2)
{
	double uu = (V1 - V2).norm();
	return sqrt((16.0 / const_pi) * T1 + (9.0 * const_pi / 4.0) * T2 + 4.0 * kv(uu));
}

double UE_bera(const double& T1, const double& T2, const Eigen::Vector3d& V1, const Eigen::Vector3d& V2)
{
	double uu = (V1 - V2).norm();
	return sqrt((4.0 / const_pi) * T1 + (64.0 / 9.0 / const_pi) * T2 + kv(uu));
}



void Setka::Calc_sourse_MF_Bera(Cell* C, unordered_map<string, double>& SOURSE,
	short int now, short int zone)
{
	short int zone_ = zone - 1;
	// Названия жидкостей
	// p, Pui1, Pui2, ....
	// H1, H2, H3, H4, .....
	// 
	// Для каждой жидкости нужна rho, T
	// Также нужна средняя скорость плазмы и Скорости всех сортов водорода
	unordered_map<string, double> rho;  // "_p", "_Pui_1", .... , "_H1", ...
	unordered_map<string, double> T;    // "_p", "_Pui_1", ...., "_H1", ...
	unordered_map<string, Eigen::Vector3d> V;    // "_p", "_H1", ...

	// Получаем переменные -----------------------------------------------------------------
	unordered_map<string, double> param;

	this->phys_param->Plasma_components(zone_, C->parameters[now], param, true);

	if (this->regim_otladki == true)
	{
		if (param.find("rho_Th") == param.end())
		{
			cout << "Error 9875641209" << endl;
			exit(-1);
		}

		if (param.find("T_Th") == param.end())
		{
			cout << "Error 9875641209" << endl;
			exit(-1);
		}
	}

	rho["_p"] = param["rho_Th"];
	T["_p"] = param["T_Th"];
	V["_p"] = Eigen::Vector3d(C->parameters[now]["Vx"], 
		C->parameters[now]["Vy"], C->parameters[now]["Vz"]);

	if (rho["_p"] <= 1e-8) rho["_p"] = 1e-8;
	if (T["_p"] <= 1e-8) T["_p"] = 1e-8;

	if(T["_p"] > 1e5) T["_p"] = 1e5;

	//cout << "B1 " << endl;

	short int pui_n = -1;
	for (const auto& nam1 : this->phys_param->pui_name)
	{
		pui_n++;
		if (this->phys_param->pui_in_zone(zone_, pui_n) == false) continue;

		rho[nam1] = C->parameters[now]["rho" + nam1];
		T[nam1] = 2.0 * C->parameters[now]["p" + nam1] / rho[nam1];
		if (rho[nam1] <= 1e-8) rho[nam1] = 1e-8;
		if (T[nam1] <= 1e-8) T[nam1] = 1e-8;
	}
	//cout << "B2 " << endl;
	for (const auto& nam2 : this->phys_param->H_name)
	{
		rho[nam2] = C->parameters[now]["rho" + nam2];
		T[nam2] = 2.0 * C->parameters[now]["p" + nam2] / max(rho[nam2], 1e-9);
		V[nam2] = Eigen::Vector3d(C->parameters[now]["Vx" + nam2],
			C->parameters[now]["Vy" + nam2], C->parameters[now]["Vz" + nam2]);
	}

	//cout << "B3 " << endl;
	// Все переменные получены -------------------------------------------------------------


	// Скорости - это симметричные функции
	unordered_map<string, double> U;    // _p_H1, _p_H2, ...., _Pui_1_H1, ...
	unordered_map<string, double> Um;   // _p_H1, _p_H2, ...., _Pui_1_H1, ...
	unordered_map<string, double> UE;   // _p_H1, _p_H2, ...., _Pui_1_H1, ...
	unordered_map<string, double> sig;  // _p_H1, _p_H2, ...., _Pui_1_H1, ...

	if (true)
	{
		pui_n = -2;
		for (const auto& nam1 : this->phys_param->p_pui_name)
		{
			pui_n++;
			if (pui_n >= 0 && this->phys_param->pui_in_zone(zone_, pui_n) == false) continue;
			for (const auto& nam2 : this->phys_param->H_name)
			{
				U[nam1 + nam2] = U_bera(T[nam1], T[nam2], V["_p"], V[nam2]);
				U[nam2 + nam1] = U_bera(T[nam2], T[nam1], V[nam2], V["_p"]);

				Um[nam1 + nam2] = Um_bera(T[nam1], T[nam2], V["_p"], V[nam2]);
				Um[nam2 + nam1] = Um_bera(T[nam2], T[nam1], V[nam2], V["_p"]);

				UE[nam1 + nam2] = UE_bera(T[nam1], T[nam2], V["_p"], V[nam2]);
				UE[nam2 + nam1] = UE_bera(T[nam2], T[nam1], V[nam2], V["_p"]);

				sig[nam1 + nam2] = kv(1.0 - this->phys_param->par_a_2 * log(U[nam1 + nam2]));
				sig[nam2 + nam1] = kv(1.0 - this->phys_param->par_a_2 * log(U[nam2 + nam1]));

				if (this->regim_otladki == true)
				{
					if (std::isnan(U[nam1 + nam2]) || std::fpclassify(U[nam1 + nam2]) == FP_SUBNORMAL ||
						std::isnan(Um[nam2 + nam1]) || std::fpclassify(Um[nam2 + nam1]) == FP_SUBNORMAL ||
						std::isnan(UE[nam1 + nam2]) || std::fpclassify(UE[nam1 + nam2]) == FP_SUBNORMAL ||
						std::isnan(sig[nam1 + nam2]) || std::fpclassify(sig[nam1 + nam2]) == FP_SUBNORMAL)
					{
						cout << "Error 4365877546g" << endl;
						whach(pui_n);
						whach(zone);
						whach(nam1);
						whach(nam2);
						whach(U[nam1 + nam2]);
						whach(Um[nam2 + nam1]);
						whach(UE[nam1 + nam2]);
						whach(sig[nam1 + nam2]);
						whach(sig[nam2 + nam1]);

						whach(rho[nam1]);
						whach(rho[nam2]);
						whach(T[nam1]);
						whach(T[nam2]);
						whach(V[nam2]);
						whach(V["_p"]);

						whach(C->parameters[now]["rho"]);
						whach(C->parameters[now]["rho_He"]);

						exit(-1);
					}
				}
			}
		}

	}
	//cout << "B4 " << endl;
	unordered_map<string, double> Hrho;
	unordered_map<string, double> HE;
	unordered_map<string, Eigen::Vector3d> Hm;
	unordered_map<string, double> Hp;

	if (true)
	{
		pui_n = -2;
		for (const auto& nam1 : this->phys_param->p_pui_name)
		{
			pui_n++;
			if (pui_n >= 0 && this->phys_param->pui_in_zone(zone_, pui_n) == false) continue;
			for (const auto& nam2 : this->phys_param->H_name)
			{
				Hrho[nam1 + nam2] = sig[nam1 + nam2] * rho[nam1] * rho[nam2] * U[nam1 + nam2];
				Hrho[nam2 + nam1] = sig[nam2 + nam1] * rho[nam1] * rho[nam2] * U[nam2 + nam1];

				Hm[nam2 + nam1] = sig[nam2 + nam1] * rho[nam1] * rho[nam2] * (U[nam2 + nam1] * V["_p"] +
					T[nam1] / Um[nam2 + nam1] * (V["_p"] - V[nam2]));
				Hm[nam1 + nam2] = sig[nam1 + nam2] * rho[nam1] * rho[nam2] * (U[nam1 + nam2] * V[nam2] -
					T[nam2] / Um[nam1 + nam2] * (V["_p"] - V[nam2]));

				HE[nam2 + nam1] = sig[nam2 + nam1] * rho[nam1] * rho[nam2] * (0.5 * U[nam2 + nam1] * V["_p"].squaredNorm() +
					T[nam1] / Um[nam2 + nam1] * V["_p"].dot(V["_p"] - V[nam2]) + 
					3.0/4.0 * T[nam1] * UE[nam2 + nam1]);

				HE[nam1 + nam2] = sig[nam1 + nam2] * rho[nam1] * rho[nam2] * (0.5 * U[nam1 + nam2] * V[nam2].squaredNorm() -
					T[nam2] / Um[nam1 + nam2] * V[nam2].dot(V["_p"] - V[nam2]) +
					3.0 / 4.0 * T[nam2] * UE[nam1 + nam2]);

				Hp[nam2 + nam1] = sig[nam2 + nam1] * rho[nam1] * rho[nam2] * this->phys_param->g1 * 3.0 / 4.0 *
					T[nam1] * UE[nam2 + nam1];

				Hp[nam1 + nam2] = sig[nam1 + nam2] * rho[nam1] * rho[nam2] * this->phys_param->g1 *
					((0.5 * U[nam1 + nam2] + T[nam2] / Um[nam1 + nam2]) * (V["_p"] - V[nam2]).squaredNorm() +
						3.0 / 4.0 * T[nam2] * UE[nam1 + nam2]);

				if (this->regim_otladki == true)
				{
					if (std::isnan(Hrho[nam1 + nam2]) || std::fpclassify(Hrho[nam1 + nam2]) == FP_SUBNORMAL ||
						std::isnan(Hm[nam1 + nam2][0]) || std::fpclassify(Hm[nam1 + nam2][0]) == FP_SUBNORMAL ||
						std::isnan(HE[nam1 + nam2]) || std::fpclassify(HE[nam1 + nam2]) == FP_SUBNORMAL ||
						std::isnan(Hp[nam2 + nam1]) || std::fpclassify(Hp[nam2 + nam1]) == FP_SUBNORMAL)
					{
						cout << "Error 65474353ywerhrt" << endl;
						whach(pui_n);
						whach(zone);
						whach(nam1);
						whach(nam2);
						whach(Hrho[nam1 + nam2]);
						whach(Hm[nam1 + nam2][0]);
						whach(HE[nam1 + nam2]);
						whach(Hp[nam2 + nam1]);
						whach(rho[nam1]);
						whach(rho[nam2]);
						whach(T[nam1]);
						whach(T[nam2]);
						whach(V[nam2]);
						whach(V["_p"]);
						whach(U[nam1 + nam2]);
						whach(U[nam2 + nam1]);
						exit(-1);
					}
				}
			}
		}

	}
	//cout << "B5 " << endl;
	if (true)
	{
		SOURSE["rho"] = 0.0;
		SOURSE["m_x"] = 0.0;
		SOURSE["m_y"] = 0.0;
		SOURSE["m_z"] = 0.0;
		SOURSE["E"] = 0.0;

		// Заполняем общие суммарные источники для плазмы  --------------------------------------
		if (this->phys_param->culc_plasma == true)
		{
			pui_n = -2;
			for (const auto& nam1 : this->phys_param->p_pui_name)
			{
				pui_n++;
				if (pui_n >= 0 && this->phys_param->pui_in_zone(zone_, pui_n) == false) continue;

				for (const auto& nam2 : this->phys_param->H_name)
				{
					auto SS = (-Hm[nam2 + nam1] + Hm[nam1 + nam2]);
					SOURSE["m_x"] += SS[0];
					SOURSE["m_y"] += SS[1];
					SOURSE["m_z"] += SS[2];

					SOURSE["E"] += (-HE[nam2 + nam1] + HE[nam1 + nam2]);

					if (this->regim_otladki == true)
					{
						if (std::isnan(SOURSE["m_x"]) || std::fpclassify(SOURSE["m_x"]) == FP_SUBNORMAL)
						{
							cout << "Error 9865342345" << endl;
							whach(SOURSE["m_x"]);
							whach(SS[0]);
							whach(pui_n);
							whach(zone);
							whach(nam1);
							whach(nam2);
							whach(Hm[nam2 + nam1]);
							whach(Hm[nam1 + nam2]);
							whach(rho[nam1]);
							whach(rho[nam2]);
							whach(T[nam1]);
							whach(T[nam2]);
							whach(U[nam1 + nam2]);
							whach(U[nam2 + nam1]);
							whach(sig[nam2 + nam1]);
							whach(sig[nam1 + nam2]);
							whach(Um[nam1 + nam2]);
							whach(Um[nam2 + nam1]);
							whach(V[nam2]);
							whach(V["_p"]);
							whach(U_bera(T[nam2], T[nam1], V[nam2], V[nam1]));
							whach(U_bera(T[nam1], T[nam2], V[nam1], V[nam2]));
							double uu = (V[nam1] - V[nam2]).norm();
							whach(uu);
							whach(((4.0 / const_pi) * (T[nam1] + T[nam2]) + kv(uu)));
							whach(sqrt((4.0 / const_pi) * (T[nam1] + T[nam2]) + kv(uu)));
							exit(-1);
						}
					}

				}
			}
		}
		//cout << "B6 " << endl;
		// Заполняем источники водорода  --------------------------------------------------------------------
		for (const auto& nam2 : this->phys_param->H_name)
		{
			SOURSE["rho" + nam2] = 0.0;
			SOURSE["m_x" + nam2] = 0.0;
			SOURSE["m_y" + nam2] = 0.0;
			SOURSE["m_z" + nam2] = 0.0;
			SOURSE["E" + nam2] = 0.0;
		}
		//cout << "B7 " << endl;
		// Сначала потери
		for (const auto& nam2 : this->phys_param->H_name)
		{
			pui_n = -2;
			for (const auto& nam1 : this->phys_param->p_pui_name)
			{
				pui_n++;
				if (pui_n >= 0 && this->phys_param->pui_in_zone(zone_, pui_n) == false) continue;

				SOURSE["rho" + nam2] += (-Hrho[nam1 + nam2]);

				SOURSE["m_x" + nam2] += (-Hm[nam1 + nam2][0]);
				SOURSE["m_y" + nam2] += (-Hm[nam1 + nam2][1]);
				SOURSE["m_z" + nam2] += (-Hm[nam1 + nam2][2]);

				SOURSE["E" + nam2] += (-HE[nam1 + nam2]);

				if (this->regim_otladki == true)
				{
					if (std::isnan(SOURSE["rho" + nam2]) || std::fpclassify(SOURSE["rho" + nam2]) == FP_SUBNORMAL ||
						std::isnan(SOURSE["m_x" + nam2]) || std::fpclassify(SOURSE["m_x" + nam2]) == FP_SUBNORMAL ||
						std::isnan(SOURSE["E" + nam2]) || std::fpclassify(SOURSE["E" + nam2]) == FP_SUBNORMAL)
					{
						cout << "Error 7847889ergyw453" << endl;
						exit(-1);
					}
				}


			}
		}
		//cout << "B8 " << endl;
		// Теперь притоки
		Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>* A;
		Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>* B;
		if (zone == 1)
		{
			A = &this->phys_param->hydrogen_arise_1;
			B = &this->phys_param->proton_arise_1;
		}
		else if (zone == 2)
		{
			A = &this->phys_param->hydrogen_arise_2;
			B = &this->phys_param->proton_arise_2;
		}
		else if (zone == 3)
		{
			A = &this->phys_param->hydrogen_arise_3;
			B = &this->phys_param->proton_arise_3;
		}
		else
		{
			A = &this->phys_param->hydrogen_arise_4;
			B = &this->phys_param->proton_arise_4;
		}

		for (short int i = 0; i < A->rows(); ++i) 
		{
			for (short int j = 0; j < A->cols(); ++j)
			{
				short iH = (*A)(i, j);
				auto nam2 = this->phys_param->H_name[i];
				auto nam1 = this->phys_param->p_pui_name[j];



				SOURSE["rho" + this->phys_param->H_name[iH - 1]] += 
					Hrho[nam2 + nam1];
				SOURSE["m_x" + this->phys_param->H_name[iH - 1]] +=
					Hm[nam2 + nam1][0];
				SOURSE["m_y" + this->phys_param->H_name[iH - 1]] +=
					Hm[nam2 + nam1][1];
				SOURSE["m_z" + this->phys_param->H_name[iH - 1]] +=
					Hm[nam2 + nam1][2];
				SOURSE["E" + this->phys_param->H_name[iH - 1]] +=
					HE[nam2 + nam1];

				if (this->regim_otladki == true)
				{
					if (std::isnan(Hrho[nam2 + nam1]) || std::fpclassify(Hrho[nam2 + nam1]) == FP_SUBNORMAL ||
						std::isnan(Hm[nam2 + nam1][1]) || std::fpclassify(Hm[nam2 + nam1][1]) == FP_SUBNORMAL ||
						std::isnan(HE[nam2 + nam1]) || std::fpclassify(HE[nam2 + nam1]) == FP_SUBNORMAL)
					{
						cout << "Error 98567345hrytrwe" << endl;
						whach(iH);
						whach(zone);
						whach(i);
						whach(j);
						whach(nam1);
						whach(nam2);
						whach(Hrho[nam2 + nam1]);
						whach(Hm[nam2 + nam1][0]);
						whach(Hm[nam2 + nam1][1]);
						whach(Hm[nam2 + nam1][2]);
						whach(HE[nam2 + nam1]);
						whach(rho[nam1]);
						whach(rho[nam2]);
						whach(T[nam1]);
						whach(T[nam2]);
						whach(V[nam2]);
						whach(V["_p"]);
						whach(U[nam1 + nam2]);
						whach(U[nam2 + nam1]);
						exit(-1);
					}
				}

			}
		}
		//cout << "B9 " << endl;
		// Заполняем источники пикапов -----------------------------------------------------
		pui_n = -1;
		for (const auto& nam2 : this->phys_param->pui_name)
		{
			pui_n++;
			if (pui_n >= 0 && this->phys_param->pui_in_zone(zone_, pui_n) == false) continue;

			SOURSE["rho" + nam2] = 0.0;
			SOURSE["p" + nam2] = 0.0;
		}

		// Сначала потери
		pui_n = -1;
		for (const auto& nam1 : this->phys_param->pui_name)
		{
			pui_n++;
			if (pui_n >= 0 && this->phys_param->pui_in_zone(zone_, pui_n) == false) continue;
			for (const auto& nam2 : this->phys_param->H_name)
			{
				SOURSE["p" + nam1] += (-Hp[nam2 + nam1]);
			}
		}
		//cout << "B10 " << endl;
		// Притоки/потери массы
		//Здесь подход отличается от остальных, хотя можно было одинаковый для всех сделать
		if (this->phys_param->is_PUI == true)
		{
			for (short int i = 0; i < B->rows(); ++i)
			{
				for (short int j = 0; j < B->cols(); ++j)
				{
					short ipui = (*B)(i, j);

					auto nam2 = this->phys_param->H_name[i];
					auto nam1 = this->phys_param->p_pui_name[j];

					if (j == 0)
					{
						if (ipui == 0) continue;
						SOURSE["rho" + this->phys_param->pui_name[ipui - 1]] +=
							Hrho[nam1 + nam2];
					}
					else
					{
						if (ipui != j)
						{
							SOURSE["rho" + nam1] -=
								Hrho[nam1 + nam2];

							if (ipui > 0)
							{
								SOURSE["rho" + this->phys_param->pui_name[ipui - 1]] +=
									Hrho[nam1 + nam2];
							}
						}
					}
				}
			}
		}
		//cout << "B11 " << endl;
		// Притоки энергии
		if (this->phys_param->is_PUI == true)
		{
			for (short int i = 0; i < B->rows(); ++i)
			{
				for (short int j = 0; j < B->cols(); ++j)
				{
					short ipui = (*B)(i, j);

					auto nam2 = this->phys_param->H_name[i];
					auto nam1 = this->phys_param->p_pui_name[j];

					if (ipui == 0) continue;

					SOURSE["p" + this->phys_param->pui_name[ipui - 1]] +=
						Hp[nam1 + nam2];
				}
			}
		}
		//cout << "B12 " << endl;
	}

	//cout << "B13 " << endl;
	double ddp = (this->phys_param->par_n_H_LISM / this->phys_param->par_Kn);
	for (auto& [key, value] : SOURSE) 
	{
		value *= ddp;  
	}

	// Источники водорода не надо умножать на концентрацию водорода, иначе получится не правильно
	for (auto& nam : this->phys_param->H_name)
	{
		SOURSE["rho" + nam] /= this->phys_param->par_n_H_LISM;
		SOURSE["m_x" + nam] /= this->phys_param->par_n_H_LISM;
		SOURSE["m_y" + nam] /= this->phys_param->par_n_H_LISM;
		SOURSE["m_z" + nam] /= this->phys_param->par_n_H_LISM;
		SOURSE["E" + nam] /= this->phys_param->par_n_H_LISM;
	}

	/*for (const auto& [key, value] : SOURSE)
	{
		cout << key << "  " << value << endl;
	}*/
	//cout << "B14 " << endl;
	//exit(-1);
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
	double p_Th;

	unordered_map<string, double> param;

	this->phys_param->Plasma_components(zone, C->parameters[now], param, true);

	rho_Th = param["rho_Th"];
	p_Th = param["p_Th"];


	double S1 = 0.0;
	double S2 = 0.0;


	int i = 0;
	for (auto& nam : this->phys_param->H_name)
	{
		U_M_H[i] = sqrt(kv(C->parameters[now]["Vx"] - C->parameters[now]["Vx" + nam])
			+ kv(C->parameters[now]["Vy"] - C->parameters[now]["Vy" + nam])
			+ kv(C->parameters[now]["Vz"] - C->parameters[now]["Vz" + nam])
			+ (64.0 / (9.0 * const_pi)) *
			(2.0 * p_Th / rho_Th
				+ 2.0 * C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam]));

		U_H[i] = sqrt(kv(C->parameters[now]["Vx"] - C->parameters[now]["Vx" + nam])
			+ kv(C->parameters[now]["Vy"] - C->parameters[now]["Vy" + nam])
			+ kv(C->parameters[now]["Vz"] - C->parameters[now]["Vz" + nam])
			+ (4.0 / const_pi) *
			(2.0 * p_Th / rho_Th
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
			(2.0 * C->parameters[now]["p" + nam] / C->parameters[now]["rho" + nam] - 2.0 * p_Th / rho_Th));
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
				+ (U_H[i] / U_M_H[i]) * (2.0 * p_Th / rho_Th));
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
	// 1, 2, 3, 4
	double rho = C->parameters[now]["rho"];
	double p = C->parameters[now]["p"];
	double u = C->parameters[now]["Vx"];
	double v = C->parameters[now]["Vy"];
	double w = C->parameters[now]["Vz"];
	double M = norm2(u, v, w) / sqrt(this->phys_param->gamma * p / rho);

	double x = C->center[now][0];
	double y = C->center[now][1];
	double z = C->center[now][2];
	double r = norm2(x, y, z);

	// Геометрическое определение
	if (true)
	{
		if (C->type == Type_cell::Zone_1 || C->type == Type_cell::Zone_2)
		{
			if (C->type == Type_cell::Zone_1)
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
			if (true && C->type == Type_cell::Zone_4)
			{
				return 4;
			}
			else
			{
				return 3;
			}
		}
	}
	else
	{
		if (C->parameters[now]["Q"] / rho < 50.0)
			//if (C->type == Type_cell::Zone_1 || C->type == Type_cell::Zone_2)
		{
			//if (M > 3.0 && r < 45)
			if (C->type == Type_cell::Zone_1)
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
			//if (M > 1 && x > 0.0)
			if (true && C->type == Type_cell::Zone_4)
			{
				return 4;
			}
			else
			{
				return 3;
			}
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

	this->Calculating_measure(0);
	this->Calculating_measure(1);


	for (auto& ce : this->All_Cell)
	{
		ce->parameters[1] = ce->parameters[0];
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

	vector<Gran*> gran_list_;
	vector<Cell*> cell_list_;

	cout << "Vibor area" << endl;
	// Если хотим отдельно считать внутреннюю и наружнюю области
	if (true)
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
	else if(true)
	{
		gran_list = &this->All_Gran;
		cell_list = &this->All_Cell;

		//gran_list = &this->Gran_outer_area;
		//cell_list = &this->Cell_outer_area;
	}
	else
	{
		for (auto& gr : this->All_Gran)
		{
			if (gr->cells.size() == 2)
			{
				if (gr->cells[0]->type == Type_cell::Zone_3 || gr->cells[0]->type == Type_cell::Zone_4
					|| gr->cells[1]->type == Type_cell::Zone_3 || gr->cells[1]->type == Type_cell::Zone_4)
				{
					gran_list_.push_back(gr);
				}
			}
			else
			{
				if (gr->cells[0]->type == Type_cell::Zone_3 || gr->cells[0]->type == Type_cell::Zone_4)
				{
					gran_list_.push_back(gr);
				}
			}
		}
		for (auto& ce : this->All_Cell)
		{
			if (ce->type == Type_cell::Zone_3 || ce->type == Type_cell::Zone_4)
			{
				cell_list_.push_back(ce);
			}
		}
		gran_list = &gran_list_;
		cell_list = &cell_list_;
	}
	cout << "END Vibor area" << endl;

	Cell* A, B;

	double time = this->phys_param->prev_step_time;  // Текущий шаг по времени
	double loc_time = this->phys_param->prev_step_time;  // Текущий шаг по времени

	double xc_min = 0.0, yc_min = 0.0, zc_min = 0.0;
	string name_min_time = "___";


	for (unsigned int step = 1; step <= steps; step++)
	{
		if (step % 100 == 0 || step == 10)
		{
			cout << "Global step = " << step << endl;
			whach(time);
			whach(xc_min);
			whach(yc_min);
			whach(zc_min);
			whach(name_min_time);
			whach(this->Cell_Center->parameters[now1]["Vx_H4"]);
			cout << "__________________" << endl;
		}

		now2 = now1;
		now1 = (now1 + 1)%2;

		this->phys_param->prev_step_time = time;
		time = loc_time;
		loc_time = 100000000.0;

		

		//omp_set_num_threads(1); // 32
		
		// Обновляем граничное условие
		if (false)
		{
			this->Init_physics_with_time(this->phys_param->ALL_Time);
		}

		this->phys_param->ALL_Time += time;
		

		// Считаем скорости граней и сразу передвигаем опорные узлы
		//if (is_inner_area == false)
		if(this->phys_param->move_setka == true)
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


		//cout << "A" << endl;
		// Расчитываем потоки через грани
		// в private не добавляются нормально vectora, надо либо обычные массивы делать, либо 
		// создавать их внутри в каждом потоке
		#pragma omp parallel for reduction(min:loc_time) schedule(dynamic)
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
		//cout << "B" << endl;

		bool print_p_less_0 = false;

		// Расчитываем законы сохранения в ячейках
		#pragma omp parallel for schedule(dynamic)
		for (size_t i_step = 0; i_step < cell_list->size(); i_step++)
		{
			auto& cell = (*cell_list)[i_step];
			double Volume = cell->volume[now1];
			double Volume2 = cell->volume[now2];
			int zone;
			double radius = norm2(cell->center[now1][0], cell->center[now1][1], cell->center[now1][2]);

			unordered_map<string, double> POTOK;

			//std::vector<double> POTOK;
			//POTOK.resize(9);


			/*boost::multi_array<double, 2> SOURSE(boost::extents[this->phys_param->H_name.size() + 1][5]);

			for (short int i = 0; i < this->phys_param->H_name.size() + 1; ++i)
			{
				for (short int j = 0; j < 5; ++j) {
					SOURSE[i][j] = 0.0;
				}
			}*/

			for (const auto& names : this->phys_param->param_names)
			{
				POTOK[names] = 0.0;
			}
			POTOK["divB"] = 0.0;
			if (this->phys_param->is_div_V_in_cell == true) POTOK["div_V"] = 0.0;

			for (auto& gran : cell->grans)
			{
				short int sign_potok = 1;
				if (gran->cells[0] != cell) sign_potok = -1;

				for (const auto& names : this->phys_param->param_names)
				{
					if (gran->parameters.find("P" + names) == gran->parameters.end()) continue; // TODO!
					POTOK[names] += sign_potok * gran->parameters["P" + names];
				}

				POTOK["divB"] += sign_potok * gran->parameters["PdivB"];

				if (this->phys_param->is_div_V_in_cell == true)
				{
					if (gran->type2 != Type_Gran_surf::TS && gran->type2 != Type_Gran_surf::BS)
					{
						POTOK["div_V"] += sign_potok * gran->parameters["Pdiv_V"];
					}
					else
					{
						if (gran->cells[0] == cell)
						{
							POTOK["div_V"] += sign_potok * gran->parameters["Pdiv_V_L"];
						}
						else
						{
							POTOK["div_V"] += sign_potok * gran->parameters["Pdiv_V_R"];
						}
					}
				}
			}

			unordered_map<string, double> SOURSE;

			zone = this->determ_zone(cell, now1);
			//this->Calc_sourse_MF(cell, SOURSE, now1, zone);
			//cout << "A1" << endl;
			//cout << "B1" << endl;
			this->Calc_sourse_MF_Bera(cell, SOURSE, now1, zone);
			//cout << "B2" << endl;
			//cout << "A2" << endl;

			double rho3, u3, v3, w3, bx3, by3, bz3, p3, Q3, rho_He3, rho_He;
			double rho, vx, vy, vz, p, bx, by, bz, dsk, Q;

			// Считаем плазму
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

				if (cell->parameters[now1].find("rho_He") != cell->parameters[now1].end())
				{
					rho_He = cell->parameters[now1]["rho_He"];
				}
				else
				{
					rho_He = 0.0;
				}



				vx = cell->parameters[now1]["Vx"];
				vy = cell->parameters[now1]["Vy"];
				vz = cell->parameters[now1]["Vz"];
				p = cell->parameters[now1]["p"];
				bx = cell->parameters[now1]["Bx"];
				by = cell->parameters[now1]["By"];
				bz = cell->parameters[now1]["Bz"];
				dsk = scalarProductFast(vx, vy, vz, bx, by, bz);

				rho3 = rho * Volume / Volume2 - time * POTOK["rho"] / Volume2;

				bool rho_null = false;
				if (rho3 < 1e-7)
				{
					rho_null = true;
					rho3 = 0.05;
					Q3 = Q / rho * rho3;
					rho_He3 = 0.0;
					cout << "Plasma  rho < 0" << endl;
					cout << cell->center[now2][0] << " " <<
						cell->center[now2][1] << " " <<
						cell->center[now2][2] << endl;
				}

				if (cell->parameters[now1].find("Q") != cell->parameters[now1].end())
				{
					Q3 = Q * Volume / Volume2 - time * POTOK["Q"] / Volume2;
				}
				else
				{
					Q3 = 0.0;
				}

				if (cell->parameters[now1].find("rho_He") != cell->parameters[now1].end())
				{
					rho_He3 = rho_He * Volume / Volume2 - time * POTOK["rho_He"] / Volume2;
				}
				else
				{
					rho_He3 = 0.0;
				}


				if (rho_He3 < 0.0)
				{
					rho_He3 = 0.0;
				}


				if (this->phys_param->is_PUI == true)
				{
					double rho_pui_sum = 0.0;
					short int pui_n = -1;
					for (const auto& nam : this->phys_param->pui_name)
					{
						pui_n++;
						if (this->phys_param->pui_in_zone(zone - 1, pui_n) == false) continue;

						if (this->regim_otladki == true)
						{
							if (POTOK.find("rho" + nam) == POTOK.end() || POTOK.find("p" + nam) == POTOK.end())
							{
								cout << "Error 9876341200" << endl;
								exit(-1);
							}
						}

						double rhorho = cell->parameters[now1]["rho" + nam]
							* Volume / Volume2 - time * POTOK["rho" + nam] / Volume2
							+ time * SOURSE["rho" + nam];

						if (rhorho <= 0.0) rhorho = 1e-7;
						cell->parameters[now2]["rho" + nam] = rhorho;
						rho_pui_sum += rhorho;

						double ppp = cell->parameters[now1]["p" + nam];
						double pppp = ppp
							* Volume / Volume2 - time * POTOK["p" + nam] / Volume2
							+ time * 0.95 * SOURSE["p" + nam]
							- time * this->phys_param->g1 * ppp * POTOK["div_V"] / Volume2;

						if (pppp < 0.0) pppp = 1e-8;

						cell->parameters[now2]["p" + nam] = pppp;

						if (this->regim_otladki == true)
						{
							if (std::isnan(pppp) || std::fpclassify(pppp) == FP_SUBNORMAL)
							{
								cout << "Error 9651211156" << endl;
								exit(-1);
							}
						}
					}


					if (rho_pui_sum > (rho3 - rho_He3) * 0.995)
					{
						double kkl = ((rho3 - rho_He3) * 0.995) / rho_pui_sum;

						pui_n = -1;
						for (const auto& nam : this->phys_param->pui_name)
						{
							pui_n++;
							if (this->phys_param->pui_in_zone(zone - 1, pui_n) == false) continue;

							cell->parameters[now2]["rho" + nam] *= kkl;
						}
					}

					
				}

				if (this->regim_otladki == true)
				{
					if (SOURSE.find("m_x") == SOURSE.end() || 
						SOURSE.find("m_y") == SOURSE.end() || 
						SOURSE.find("m_z") == SOURSE.end() || 
						SOURSE.find("E") == SOURSE.end())
					{
						cout << "Error 3412089634" << endl;
						exit(-1);
					}

					if (std::isnan(SOURSE["m_x"]) || std::fpclassify(SOURSE["m_x"]) == FP_SUBNORMAL ||
						std::isnan(SOURSE["m_y"]) || std::fpclassify(SOURSE["m_y"]) == FP_SUBNORMAL ||
						std::isnan(SOURSE["m_z"]) || std::fpclassify(SOURSE["m_z"]) == FP_SUBNORMAL ||
						std::isnan(SOURSE["E"]) || std::fpclassify(SOURSE["E"]) == FP_SUBNORMAL)
					{
						cout << "Error 8634108946" << endl;
						cout << SOURSE["m_x"] << endl;
						cout << SOURSE["m_y"] << endl;
						cout << SOURSE["m_z"] << endl;
						cout << SOURSE["E"] << endl;
						cout << zone << endl;
						cout << cell->center[0][0] << " " << cell->center[0][1] << " " <<
							cell->center[0][2] << endl;
						exit(-1);
					}

				}

				if (rho_null == false)
				{
					u3 = (rho * vx * Volume / Volume2 - time * (POTOK["Vx"] + (bx / cpi4) * POTOK["divB"]) / Volume2
						+ time * SOURSE["m_x"]) / rho3;
					v3 = (rho * vy * Volume / Volume2 - time * (POTOK["Vy"] + (by / cpi4) * POTOK["divB"]) / Volume2
						+ time * SOURSE["m_y"]) / rho3;
					w3 = (rho * vz * Volume / Volume2 - time * (POTOK["Vz"] + (bz / cpi4) * POTOK["divB"]) / Volume2
						+ time * SOURSE["m_z"]) / rho3;

					bx3 = bx * Volume / Volume2 - time * (POTOK["Bx"] + vx * POTOK["divB"]) / Volume2;
					by3 = by * Volume / Volume2 - time * (POTOK["By"] + vy * POTOK["divB"]) / Volume2;
					bz3 = bz * Volume / Volume2 - time * (POTOK["Bz"] + vz * POTOK["divB"]) / Volume2;

					p3 = (((p / this->phys_param->g1 + 0.5 * rho * kvv(vx, vy, vz) + kvv(bx, by, bz) / 25.13274122871834590768) * Volume / Volume2
						- time * (POTOK["p"] + (dsk / cpi4) * POTOK["divB"]) / Volume2 + time * SOURSE["E"]) -
						0.5 * rho3 * kvv(u3, v3, w3) - kvv(bx3, by3, bz3) / 25.13274122871834590768) * this->phys_param->g1;

					// Эффективная магнитная диссипация
					//if (cell->type == Type_cell::Zone_2) // && p3/(kvv(bx3, by3, bz3) / (8.0 * const_pi)) < 1.0)
					if(false)
					{
						double p4, bx4, by4, bz4;
						double tau = 2.5; // 2.5456; // 10.0

						bx4 = bx3 - time * bx3 / tau;
						by4 = by3 - time * by3 / tau;
						bz4 = bz3 - time * bz3 / tau;
						p4 = p3 + (time / tau - kv(time / tau)) * (this->phys_param->g1) * kvv(bx3, by3, bz3) / (4.0 * const_pi);

						bx3 = bx4;
						by3 = by4;
						bz3 = bz4;
						p3 = p4;
					}


				}
				else
				{
					u3 = vx;
					v3 = vy;
					w3 = vz;
					bx3 = bx;
					by3 = by;
					bz3 = bz;
					p3 = rho3/2.0;
				}

				if (std::isnan(rho3) || std::fpclassify(rho3) == FP_SUBNORMAL ||
					std::isnan(Q3) || std::fpclassify(Q3) == FP_SUBNORMAL)
				{
					cout << "Error  9851234578" << endl;
					whach(Volume);
					whach(Volume2);
					whach(rho);
					whach(time);
					whach(Q);
					whach(rho_He);

					for (short unsigned int ik = 0; ik < 9; ik++)
					{
						cout << "ik: " << endl;// << POTOK[ik] << endl;
					}
					exit(-1);
				}


				if (p3 < 0.00000001)
				{
					p3 = 0.000001;

					if (step % 25 == 0 && print_p_less_0 == false)
					{
						#pragma omp critical (firstf) 
						{
							if (print_p_less_0 == false)
							{
								cout << "Plasma  p < 0" << endl;
								cout << cell->center[now2][0] << " " <<
									cell->center[now2][1] << " " <<
									cell->center[now2][2] << endl;
								print_p_less_0 = true;
							}
						}
					}
					

				}


				if (rho_He3 > rho3)
				{
					rho_He3 = rho3 * 0.1;
				}

				cell->parameters[now2]["rho"] = rho3;
				if (cell->parameters[now1].find("Q") != cell->parameters[now1].end())
				{
					cell->parameters[now2]["Q"] = Q3;
				}
				if (cell->parameters[now1].find("rho_He") != cell->parameters[now1].end())
				{
					cell->parameters[now2]["rho_He"] = rho_He3;
				}
				cell->parameters[now2]["Vx"] = u3;
				cell->parameters[now2]["Vy"] = v3;
				cell->parameters[now2]["Vz"] = w3;
				if (cell->is_TVD == true)
				{
					cell->parameters[now2]["Bx"] = bx3;
					cell->parameters[now2]["By"] = by3;
					cell->parameters[now2]["Bz"] = bz3;
				}
				cell->parameters[now2]["p"] = p3;
			}

			//cout << "B3" << endl;
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
			if (this->phys_param->culc_atoms == true)
			{
				for (auto& nam : this->phys_param->H_name)
				{
					if (this->phys_param->Culc_hidrogen[nam] == false) continue;
					rho = cell->parameters[now1]["rho" + nam];
					vx = cell->parameters[now1]["Vx" + nam];
					vy = cell->parameters[now1]["Vy" + nam];
					vz = cell->parameters[now1]["Vz" + nam];
					p = cell->parameters[now1]["p" + nam];

					if (this->regim_otladki == true)
					{
						if (SOURSE.find("m_x" + nam) == SOURSE.end() ||
							SOURSE.find("m_y" + nam) == SOURSE.end() ||
							SOURSE.find("m_z" + nam) == SOURSE.end() ||
							SOURSE.find("E" + nam) == SOURSE.end())
						{
							cout << "Error 6548712112" << endl;
							exit(-1);
						}

						if (std::isnan(SOURSE["m_x" + nam]) || std::fpclassify(SOURSE["m_x" + nam]) == FP_SUBNORMAL ||
							std::isnan(SOURSE["m_y" + nam]) || std::fpclassify(SOURSE["m_y" + nam]) == FP_SUBNORMAL ||
							std::isnan(SOURSE["m_z" + nam]) || std::fpclassify(SOURSE["m_z" + nam]) == FP_SUBNORMAL ||
							std::isnan(SOURSE["E" + nam]) || std::fpclassify(SOURSE["E" + nam]) == FP_SUBNORMAL)
						{
							cout << "Error 726t4366352" << endl;
							exit(-1);
						}

					}

					if (rho < 1e-7 || radius <= 1.0)  // Отключаем источники, так как этой жидкости фактически нет
					{
						SOURSE["rho" + nam] = 0.0;
						SOURSE["m_x" + nam] = 0.0;
						SOURSE["m_y" + nam] = 0.0;
						SOURSE["m_z" + nam] = 0.0;
						SOURSE["E" + nam] = 0.0;

					}

					rho3 = rho * Volume / Volume2 - time * POTOK["rho" + nam] / Volume2 
						+ time * SOURSE["rho" + nam];

					if (rho3 < 1e-8)
					{
						rho3 = 1e-8;
						//cout << "Hidrogen  rho < 0  " << nam << endl;
					}

					u3 = (rho * vx * Volume / Volume2 - time * (POTOK["Vx" + nam]) / Volume2
						+ time * SOURSE["m_x" + nam]) / rho3;
					v3 = (rho * vy * Volume / Volume2 - time * (POTOK["Vy" + nam]) / Volume2
						+ time * SOURSE["m_y" + nam]) / rho3;
					w3 = (rho * vz * Volume / Volume2 - time * (POTOK["Vz" + nam]) / Volume2
						+ time * SOURSE["m_z" + nam]) / rho3;


					p3 = (((p / this->phys_param->g1 + 0.5 * rho * kvv(vx, vy, vz)) * Volume / Volume2
						- time * (POTOK["p" + nam]) / Volume2 + time * SOURSE["E" + nam]) -
						0.5 * rho3 * kvv(u3, v3, w3)) * this->phys_param->g1;

					if (p3 < 1e-8)
					{
						p3 = 1e-8;
						//cout << "Hidrogen  p < 0  " << nam << endl;
						//cout << "Center = " << cell->center[now2][0] << " " <<
						//	cell->center[now2][1] << " " <<
						//	cell->center[now2][2] << endl;
					}

					if (norm2(u3, v3, w3) > 1e3)
					{
						double ddf = norm2(u3, v3, w3);
						u3 = u3 / ddf * 1e3;
						v3 = v3 / ddf * 1e3;
						w3 = w3 / ddf * 1e3;
					}

					cell->parameters[now2]["rho" + nam] = rho3;
					cell->parameters[now2]["Vx" + nam] = u3;
					cell->parameters[now2]["Vy" + nam] = v3;
					cell->parameters[now2]["Vz" + nam] = w3;
					cell->parameters[now2]["p" + nam] = p3;

					if (this->regim_otladki == true)
					{
						if (norm2(u3, v3, w3) > 1e20)
						{
							cout << "Error 9875463498" << endl;
							whach(nam);
							whach(zone);
							whach(cell->center[0][0]);
							whach(cell->center[0][1]);
							whach(cell->center[0][2]);

							whach(u3);
							whach(v3);
							whach(w3);
							whach(vx);
							whach(vy);
							whach(vz);
							whach(rho3);
							whach(rho);
							whach(p);
							whach(p3);

							whach(SOURSE["rho" + nam]);
							whach(POTOK["rho" + nam]);

							whach(SOURSE["m_x" + nam]);
							whach(POTOK["Vx" + nam]);

							whach(SOURSE["m_y" + nam]);
							whach(POTOK["Vy" + nam]);

							whach(SOURSE["m_z" + nam]);
							whach(POTOK["Vz" + nam]);

							whach(SOURSE["E" + nam]);
							whach(POTOK["p" + nam]);

							exit(-1);
						}
					}

					i++;
				}
			}
		}

		#pragma omp barrier

		//cout << "C" << endl;

		// Увеличим давление PUI за ударной волной
		if (is_inner_area == false && this->phys_param->is_PUI == true)
		{
			#pragma omp parallel for
			for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
			{
				auto gr = this->Gran_TS[i_step];
				auto A = gr->cells[0];
				auto B = gr->cells[1];

				B->parameters[now2]["p_Pui_1"] = A->parameters[now2]["p_Pui_1"] *
					kv(B->parameters[now2]["rho"] / A->parameters[now2]["rho"]);
			}
		}

		// Считаем фиктивную центральную ячейку
		if (this->phys_param->culc_atoms == true && is_inner_area == true)
		{
			auto& cell = this->Cell_Center;
			double Volume = 4.0 * const_pi * kyb(this->geo->R0)/3.0;

			unordered_map<string, double> POTOK_F;

			for (const auto& names: this->phys_param->H_name)
			{
				POTOK_F["rho" + names] = 0.0;
				POTOK_F["Vx" + names] = 0.0;
				POTOK_F["Vy" + names] = 0.0;
				POTOK_F["Vz" + names] = 0.0;
				POTOK_F["p" + names] = 0.0;
			}


			for (auto& gran : this->All_boundary_Gran)
			{
				if (gran->type != Type_Gran::Inner_Hard) continue;

				short int sign_potok = -1;

				for (const auto& names : this->phys_param->H_name)
				{
					POTOK_F["rho" + names] += sign_potok * gran->parameters["Prho" + names];
					POTOK_F["Vx" + names] += sign_potok * gran->parameters["PVx" + names];
					POTOK_F["Vy" + names] += sign_potok * gran->parameters["PVy" + names];
					POTOK_F["Vz" + names] += sign_potok * gran->parameters["PVz" + names];
					POTOK_F["p" + names] += sign_potok * gran->parameters["Pp" + names];
				}
			}

			double rho3, u3, v3, w3, p3;
			double rho, vx, vy, vz, p;

			// Теперь считаем для остальных жидкостей
			int i = 1;
			for (auto& nam : this->phys_param->H_name)
			{
				if (this->phys_param->Culc_hidrogen[nam] == false) continue;

				if (nam == "_H1") continue;
				if (nam == "_H5") continue;

				rho = cell->parameters[now1]["rho" + nam];
				vx = cell->parameters[now1]["Vx" + nam];
				vy = cell->parameters[now1]["Vy" + nam];
				vz = cell->parameters[now1]["Vz" + nam];
				p = cell->parameters[now1]["p" + nam];

				rho3 = rho - time * POTOK_F["rho" + nam] / Volume;

				if (rho3 < 1e-8)
				{
					rho3 = 1e-8;
				}

				u3 = (rho * vx - time * (POTOK_F["Vx" + nam]) / Volume) / rho3;
				v3 = (rho * vy - time * (POTOK_F["Vy" + nam]) / Volume) / rho3;
				w3 = (rho * vz - time * (POTOK_F["Vz" + nam]) / Volume) / rho3;


				p3 = (((p / this->phys_param->g1 + 0.5 * rho * kvv(vx, vy, vz))
					- time * (POTOK_F["p" + nam]) / Volume) -
					0.5 * rho3 * kvv(u3, v3, w3)) * this->phys_param->g1;

				if (p3 < 1e-8)
				{
					p3 = 1e-8;
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
	int now2 = (now + 1) % 2;

	double area = (gr->area[now] + gr->area[now2])/2.0;
	name = "______";
	short int metod_ = metod;

	unordered_map<string, double> par_left, par_right;
	double w = 0.0;
	if(this->phys_param->move_setka == true) w = gr->culc_velosity(now, time);
	double loc_time = 10000000.0;

	double dist;
	if (gr->type == Type_Gran::Us)
	{
		// Это для расчёта шага по времени
		auto A = gr->cells[0];
		auto B = gr->cells[1];
		//dist = norm2(A->center[now][0] - B->center[now][0],
		//	A->center[now][1] - B->center[now][1],
		//	A->center[now][2] - B->center[now][2]) / 2.0;
		dist = min(norm2(A->center[now][0] - gr->center[now][0],
			A->center[now][1] - gr->center[now][1],
			A->center[now][2] - gr->center[now][2]), norm2(gr->center[now][0] - B->center[now][0],
				gr->center[now][1] - B->center[now][1],
				gr->center[now][2] - B->center[now][2]));
	}

	if (this->phys_param->culc_plasma == true)
	{
		this->Snos_on_Gran(gr, par_left, par_right, now, true);
	}

	if (this->phys_param->culc_atoms == true)
	{
		this->Snos_on_Gran(gr, par_left, par_right, now, false);
	}

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

		if (par_left.find("rho_He") != par_left.end())
		{
			konvect_left.push_back(par_left["rho_He"]);
			konvect_right.push_back(par_right["rho_He"]);
			konvect.push_back(0.0);
		}

		if (this->phys_param->is_PUI == true)
		{
			for (const auto& nam : this->phys_param->pui_name)
			{
				konvect_left.push_back(par_left["rho" + nam]);
				konvect_right.push_back(par_right["rho" + nam]);
				konvect.push_back(0.0);

				konvect_left.push_back(par_left["p" + nam]);
				konvect_right.push_back(par_right["p" + nam]);
				konvect.push_back(0.0);
			}
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
				n, qqq, dsl, dsr, dsc, w, this->phys_param->contact_hard);

			qqq[5] = qqq[6] = qqq[7] = 0.0;
			konvect[0] = 0.0;
			konvect[1] = 0.0;
		}
		else
		{
			bool left_ydar = false;
			bool contact = false;

			if (gr->type2 == Type_Gran_surf::TS && this->phys_param->TS_hard)
				left_ydar = true;

			if (gr->type2 == Type_Gran_surf::HP && this->phys_param->contact_hard)
				contact = true;

			this->phys_param->chlld(metod_, gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 0, // 1
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option, left_ydar, contact);
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

		gr->parameters["Prho"] = qqq[0] * area;
		gr->parameters["PVx"] = qqq[1] * area;
		gr->parameters["PVy"] = qqq[2] * area;
		gr->parameters["PVz"] = qqq[3] * area;
		gr->parameters["Pp"] = qqq[4] * area;
		gr->parameters["PBx"] = qqq[5] * area;
		gr->parameters["PBy"] = qqq[6] * area;
		gr->parameters["PBz"] = qqq[7] * area;
		gr->parameters["PdivB"] = 0.5 * scalarProductFast(gr->normal[now][0],
			gr->normal[now][1], gr->normal[now][2],
			qqq1[5] + qqq2[5], qqq1[6] + qqq2[6], qqq1[7] + qqq2[7]) * area;

		if (this->phys_param->is_div_V_in_cell == true)
		{
			if (gr->type2 == Type_Gran_surf::Us)
			{
				gr->parameters["Pdiv_V"] = 0.5 * scalarProductFast(gr->normal[now][0],
					gr->normal[now][1], gr->normal[now][2],
					qqq1[1] + qqq2[1], qqq1[2] + qqq2[2], qqq1[3] + qqq2[3]) * area;
			}
			else if (gr->type2 == Type_Gran_surf::TS || gr->type2 == Type_Gran_surf::BS)
			{
				gr->parameters["Pdiv_V_L"] = scalarProductFast(gr->normal[now][0],
					gr->normal[now][1], gr->normal[now][2],
					qqq1[1], qqq1[2], qqq1[3]) * area;
				gr->parameters["Pdiv_V_R"] = scalarProductFast(gr->normal[now][0],
					gr->normal[now][1], gr->normal[now][2],
					qqq2[1], qqq2[2], qqq2[3]) * area;
			}
		}

		// Для контакта поток Bn равен нулю
		if (gr->type2 == Type_Gran_surf::HP && this->phys_param->bn_in_p_on_HP == true)
		{
			gr->parameters["PdivB"] = 0.0;
		}

		// Для контакта поток Vn равен нулю
		if (gr->type2 == Type_Gran_surf::HP && this->phys_param->is_div_V_in_cell == true)
		{
			gr->parameters["Pdiv_V"] = 0.0;
		}



		if (par_left.find("Q") != par_left.end())
		{
			gr->parameters["PQ"] = konvect[0] * area; // Может быть неправльный порядок при других конвективных переменных
		}

		if (par_left.find("rho_He") != par_left.end())
		{
			gr->parameters["Prho_He"] = konvect[1] * area; // Может быть неправльный порядок при других конвективных переменных
		}

		if (this->phys_param->is_PUI == true)
		{
			uint8_t ip = 2;
			for (const auto& nam : this->phys_param->pui_name)
			{
				gr->parameters["Prho" + nam] = konvect[ip] * area;
				ip++;
				gr->parameters["Pp" + nam] = konvect[ip] * area;
				ip++;
			}
		}

	}

	konvect_left.clear();
	konvect_right.clear();
	konvect.clear();

	if (this->phys_param->culc_atoms == true)
	{
		for (auto& nam : this->phys_param->H_name)
		{
			if (this->phys_param->Culc_hidrogen[nam] == false) 
			{
				continue;
			}

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
				w, qqq1, qqq2, qqq, false, 0, //1
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

			gr->parameters["Prho" + nam] = qqq[0] * area;
			gr->parameters["PVx" + nam] = qqq[1] * area;
			gr->parameters["PVy" + nam] = qqq[2] * area;
			gr->parameters["PVz" + nam] = qqq[3] * area;
			gr->parameters["Pp" + nam] = qqq[4] * area;

		}
	}

	return loc_time;
}

void Setka::Save_for_interpolate(string filename, bool razriv)
{
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		cout << "Error 097564537  Can not  open file to writing: " + filename << endl;
		exit(-1);
	}

	// Записываем есть ли особенная интерполяция на разрывах (или всё сплошным образом)
	out.write(reinterpret_cast<const char*>(&razriv), sizeof(bool));

	// Записываем до какого расстояния слева выделяется HP
	out.write(reinterpret_cast<const char*>(&this->geo->L6), sizeof(double));

	// Добавляем ещё переменные для вывода  "BB/8pi"
	if (true)
	{
		this->phys_param->param_names.push_back("BB/8pi");
		for (const auto& Cel : this->All_Cell)
		{
			Cel->parameters[0]["BB/8pi"] = kvv(Cel->parameters[0]["Bx"],
				Cel->parameters[0]["By"], Cel->parameters[0]["Bz"]) / (8.0 * const_pi);
		}
		this->Cell_Center->parameters[0]["BB/8pi"] = kvv(this->Cell_Center->parameters[0]["Bx"],
			this->Cell_Center->parameters[0]["By"], this->Cell_Center->parameters[0]["Bz"]) / (8.0 * const_pi);
	}


	// Записываем количество строк
	size_t size = this->phys_param->param_names.size() + 1;
	out.write(reinterpret_cast<const char*>(&size), sizeof(size));

	// Записываем каждую строку
	for (const auto& str : this->phys_param->param_names) {
		// Сначала записываем длину строки
		size_t str_size = str.size();
		out.write(reinterpret_cast<const char*>(&str_size), sizeof(str_size));
		// Затем саму строку
		out.write(str.data(), str_size);
	}

	cout << "All parameters (send): " << endl;
	for (const auto& i : this->phys_param->param_names)
	{
		cout << i << "  ";
	}
	cout << endl;

	// Добавляем геометрическую зону
	if (true)
	{
		string str = "zone_geo";
		size_t str_size = str.size();
		out.write(reinterpret_cast<const char*>(&str_size), sizeof(str_size));
		out.write(str.data(), str_size);
	}

	// Считаем сколько дополнительных ячеек будет на внешней границе
	unsigned int gr_b = 0;

	if (razriv == true)
	{
		//  Записываем первую зону
		if (true)
		{
			for (const auto& Cel : this->All_Cell)
			{
				Cel->is_need = static_cast<short int>(Cel->type);
				if (Cel->is_need == 1) gr_b++;
			}

			// Записываем количество ячеек
			size = gr_b + 1 + 2 * this->Gran_TS.size(); // Центр + доп точки на границе и за ней
			out.write(reinterpret_cast<const char*>(&size), sizeof(size));

			for (const auto& Cel : this->All_Cell)
			{
				if (Cel->is_need != 1) continue;
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


				double zzz = 1.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

			unordered_map<string, double> par_left, par_right;
			for (const auto& gr : this->Gran_TS)
			{
				auto C1 = gr->cells[0];
				auto C2 = gr->cells[1];

				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				this->Snos_on_Gran(gr, par_left, par_right, 0, false);



				Eigen::Vector3d A1, A2, A3;
				A1 << C1->center[0][0], C1->center[0][1], C1->center[0][2];

				double aa = gr->center[0][0];
				double bb = gr->center[0][1];
				double cc = gr->center[0][2];
				A2 << aa, bb, cc;
				out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = 0.0;

					if (par_left.find(i) != par_left.end())
					{
						aa = par_left[i];
					}

					out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				}

				double zzz = 1.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));

				A3 = A2 + (A2 - A1);
				out.write(reinterpret_cast<const char*>(&A3[0]), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&A3[1]), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&A3[2]), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = par_left[i];
					bb = C1->parameters[0][i];
					cc = aa + (aa - bb);

					out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));
				}

				zzz = 1.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
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

				//cout << "================  " << this->Cell_Center->parameters[0]["rho"] << endl;

				double zzz = 1.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
				//cout << "================  " << zzz << endl;
			}
		}

		// Записываем вторую зону
		if (true)
		{
			gr_b = 0;
			for (const auto& Cel : this->All_Cell)
			{
				Cel->is_need = static_cast<short int>(Cel->type);
				if ((Cel->is_need == 2 || Cel->is_need == 3) && Cel->center[0][0] >= this->geo->L6 - 50.0) gr_b++;
			}

			// Записываем количество ячеек
			size = gr_b + 2 * this->Gran_TS.size() + 0 * this->Gran_HP.size(); //доп точки на TS и HP
			out.write(reinterpret_cast<const char*>(&size), sizeof(size));

			for (const auto& Cel : this->All_Cell)
			{
				if (Cel->is_need == 1 || Cel->is_need == 4 || Cel->center[0][0] < this->geo->L6 - 50.0) continue;
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


				double zzz = static_cast<short int>(Cel->type);
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

			unordered_map<string, double> par_left, par_right;
			for (const auto& gr : this->Gran_TS)
			{
				auto C1 = gr->cells[0];
				auto C2 = gr->cells[1];
				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				this->Snos_on_Gran(gr, par_left, par_right, 0, false);


				Eigen::Vector3d A1, A2, A3;
				A1 << C2->center[0][0], C2->center[0][1], C2->center[0][2];

				double aa = gr->center[0][0];
				double bb = gr->center[0][1];
				double cc = gr->center[0][2];
				A2 << aa, bb, cc;
				out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = 0.0;

					if (par_right.find(i) != par_right.end())
					{
						aa = par_right[i];
					}

					out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				}

				double zzz = 2.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));

				A3 = A2 + (A2 - A1);
				out.write(reinterpret_cast<const char*>(&A3[0]), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&A3[1]), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&A3[2]), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = par_right[i];
					bb = C2->parameters[0][i];
					cc = aa + (aa - bb);

					out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));
				}

				zzz = 2.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

			if (false)
			{
				for (const auto& gr : this->Gran_HP)
				{
					auto C1 = gr->cells[0];
					auto C2 = gr->cells[1];

					this->Snos_on_Gran(gr, par_left, par_right, 0, true);
					this->Snos_on_Gran(gr, par_left, par_right, 0, false);


					Eigen::Vector3d A1, A2, A3, nn;
					A1 << C2->center[0][0], C2->center[0][1], C2->center[0][2];

					double aa = gr->center[0][0];
					double bb = gr->center[0][1];
					double cc = gr->center[0][2];
					A2 << aa, bb, cc;
					nn << gr->normal[0][0], gr->normal[0][1], gr->normal[0][2];
					double zzz;

					/*out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
					out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
					out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

					for (const auto& i : this->phys_param->param_names)
					{
						aa = 0.0;

						if (par_left.find(i) != par_left.end())
						{
							aa = par_left[i];
						}

						out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
					}

					zzz = 2.0;
					out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));*/

					for (short int io = 1; io <= 1; io++)
					{

						A3 = A1;

						out.write(reinterpret_cast<const char*>(&A3[0]), sizeof(double));
						out.write(reinterpret_cast<const char*>(&A3[1]), sizeof(double));
						out.write(reinterpret_cast<const char*>(&A3[2]), sizeof(double));

						for (const auto& i : this->phys_param->param_names)
						{
							aa = par_left[i];
							bb = C1->parameters[0][i];
							cc = bb; //bb + io * (aa - bb);

							out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));
						}

						zzz = 2.0;
						out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
					}
				}
			}
		}

		// Записываем третью зону
		if (true)
		{
			gr_b = 0;
			for (const auto& Cel : this->All_Cell)
			{
				Cel->is_need = static_cast<short int>(Cel->type);
				if (Cel->is_need >= 2 && Cel->center[0][0] >= -20.0) gr_b++;
			}

			// Записываем количество ячеек
			size = gr_b + 0 * this->Gran_HP.size() + 2 * this->Gran_BS.size(); //доп точки на TS и HP
			out.write(reinterpret_cast<const char*>(&size), sizeof(size));

			for (const auto& Cel : this->All_Cell)
			{
				if (Cel->is_need < 2 || Cel->center[0][0] < -20.0) continue;
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


				double zzz = static_cast<short int>(Cel->type);
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

			unordered_map<string, double> par_left, par_right;

			if (false)
			{
				for (const auto& gr : this->Gran_HP)
				{
					auto C1 = gr->cells[0];
					auto C2 = gr->cells[1];

					this->Snos_on_Gran(gr, par_left, par_right, 0, true);
					this->Snos_on_Gran(gr, par_left, par_right, 0, false);


					Eigen::Vector3d A1, A2, A3;
					A1 << C2->center[0][0], C2->center[0][1], C2->center[0][2];

					double aa = gr->center[0][0];
					double bb = gr->center[0][1];
					double cc = gr->center[0][2];
					A2 << aa, bb, cc;
					out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
					out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
					out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

					for (const auto& i : this->phys_param->param_names)
					{
						aa = 0.0;

						if (par_right.find(i) != par_right.end())
						{
							aa = par_right[i];
						}

						out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
					}

					double zzz = 3.0;
					out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));

					for (short int io = 1; io <= 5; io++)
					{
						A3 = A2 + io * (A2 - A1);
						out.write(reinterpret_cast<const char*>(&A3[0]), sizeof(aa));
						out.write(reinterpret_cast<const char*>(&A3[1]), sizeof(bb));
						out.write(reinterpret_cast<const char*>(&A3[2]), sizeof(cc));

						for (const auto& i : this->phys_param->param_names)
						{
							aa = par_right[i];
							bb = C2->parameters[0][i];
							cc = aa;// + io * (aa - bb);

							out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));
						}

						zzz = 3.0;
						out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
					}
				}
			}

			for (const auto& gr : this->Gran_BS)
			{
				auto C1 = gr->cells[0];
				auto C2 = gr->cells[1];
				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				this->Snos_on_Gran(gr, par_left, par_right, 0, false);

				Eigen::Vector3d A1, A2, A3;
				A1 << C1->center[0][0], C1->center[0][1], C1->center[0][2];

				double aa = gr->center[0][0];
				double bb = gr->center[0][1];
				double cc = gr->center[0][2];
				A2 << aa, bb, cc;
				out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = 0.0;

					if (par_left.find(i) != par_left.end())
					{
						aa = par_left[i];
					}

					out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				}

				double zzz = 3.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));

				A3 = A2 + (A2 - A1);
				out.write(reinterpret_cast<const char*>(&A3[0]), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&A3[1]), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&A3[2]), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = par_left[i];
					bb = C1->parameters[0][i];
					cc = aa + (aa - bb);

					out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));
				}

				zzz = 3.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

		}

		// Записываем четвёртую зону
		if (true)
		{
			gr_b = 0;
			for (const auto& Cel : this->All_Cell)
			{
				Cel->is_need = static_cast<short int>(Cel->type);
				if ((Cel->is_need == 4 || Cel->is_need == 3) && Cel->center[0][0] >= -100.0) gr_b++;
			}

			// Записываем количество ячеек
			size = gr_b + 3 * this->Gran_BS.size(); //доп точки на TS и HP
			out.write(reinterpret_cast<const char*>(&size), sizeof(size));

			for (const auto& Cel : this->All_Cell)
			{
				if ((Cel->is_need != 4 && Cel->is_need != 3) || Cel->center[0][0] < -100.0) continue;
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


				double zzz = static_cast<short int>(Cel->type);
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

			unordered_map<string, double> par_left, par_right;
			for (const auto& gr : this->Gran_BS)
			{
				auto C1 = gr->cells[0];
				auto C2 = gr->cells[1];
				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				this->Snos_on_Gran(gr, par_left, par_right, 0, false);


				Eigen::Vector3d A1, A2, A3;
				A1 << C2->center[0][0], C2->center[0][1], C2->center[0][2];

				double aa = gr->center[0][0];
				double bb = gr->center[0][1];
				double cc = gr->center[0][2];
				A2 << aa, bb, cc;
				out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = 0.0;

					if (par_right.find(i) != par_right.end())
					{
						aa = par_right[i];
					}

					out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				}

				double zzz = 4.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));

				A3 = A2 + (A2 - A1);
				out.write(reinterpret_cast<const char*>(&A3[0]), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&A3[1]), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&A3[2]), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = par_right[i];
					bb = C2->parameters[0][i];
					cc = aa + (aa - bb);

					out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));
				}

				zzz = 4.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));

				A3 = A2 + 2 * (A2 - A1);
				out.write(reinterpret_cast<const char*>(&A3[0]), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&A3[1]), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&A3[2]), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = par_right[i];
					bb = C2->parameters[0][i];
					cc = aa + 2 * (aa - bb);

					out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));
				}

				zzz = 4.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

		}

		// Записываем пятую зону
		if (true)
		{
			gr_b = 0;
			for (const auto& Cel : this->All_Cell)
			{
				Cel->is_need = static_cast<short int>(Cel->type);
				if ((Cel->is_need >= 2) && Cel->center[0][0] >= this->geo->L6 - 50.0
					&& Cel->center[0][0] <= 60.0) gr_b++;
			}

			// Записываем количество ячеек
			size = gr_b + 2 * this->Gran_HP.size(); //доп точки на TS и HP
			out.write(reinterpret_cast<const char*>(&size), sizeof(size));

			for (const auto& Cel : this->All_Cell)
			{
				if ((Cel->is_need == 1) || Cel->center[0][0] < this->geo->L6 - 50.0
					|| Cel->center[0][0] > 60.0) continue;
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


				double zzz = static_cast<short int>(Cel->type);
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

			unordered_map<string, double> par_left, par_right;
			for (const auto& gr : this->Gran_HP)
			{
				auto C1 = gr->cells[0];
				auto C2 = gr->cells[1];
				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				this->Snos_on_Gran(gr, par_left, par_right, 0, false);

				Eigen::Vector3d A1, A2, A3;
				A1 << C2->center[0][0], C2->center[0][1], C2->center[0][2];

				double aa = gr->center[0][0];
				double bb = gr->center[0][1];
				double cc = gr->center[0][2];
				A2 << aa, bb, cc;
				out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&bb), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = 0.0;

					if (par_right.find(i) != par_right.end())
					{
						aa = par_right[i];
					}

					out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
				}

				double zzz = 3.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));

				A3 = A2 + (A2 - A1);
				out.write(reinterpret_cast<const char*>(&A3[0]), sizeof(aa));
				out.write(reinterpret_cast<const char*>(&A3[1]), sizeof(bb));
				out.write(reinterpret_cast<const char*>(&A3[2]), sizeof(cc));

				for (const auto& i : this->phys_param->param_names)
				{
					aa = par_right[i];
					bb = C2->parameters[0][i];
					cc = aa + (aa - bb);

					out.write(reinterpret_cast<const char*>(&cc), sizeof(cc));
				}

				zzz = 3.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}



		}

		// Записываем шестую зону
		if (true)
	{
		gr_b = 0;
		for (const auto& Cel : this->All_Cell)
		{
			Cel->is_need = static_cast<short int>(Cel->type);
			if (Cel->center[0][0] <= this->geo->L6 + 50.0) gr_b++;
		}

		// Записываем количество ячеек
		size = gr_b; //доп точки на TS и HP
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));

		for (const auto& Cel : this->All_Cell)
		{
			if (Cel->center[0][0] > this->geo->L6 + 50.0) continue;
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


			double zzz = static_cast<short int>(Cel->type);
			out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
		}

	}
	}
	else
	{
		if (true)
		{
			for (const auto& Cel : this->All_Cell)
			{
				Cel->is_need = static_cast<short int>(Cel->type);
			}

			// Записываем количество ячеек
			size = this->All_Cell.size(); // + Центр
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


				double zzz = static_cast<double>(Cel->is_need);
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

			// Записываем центральную точку
			if (false)
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

				//cout << "================  " << this->Cell_Center->parameters[0]["rho"] << endl;

				double zzz = 1.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
				//cout << "================  " << zzz << endl;
			}
		}
	}

	// Записываем дополнительные точки (на небольшом удалении от внешней границы, чтоб 
	// убрать артефакты в интерполяции
	if (false)
	{
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
				double zzz = static_cast<double>(A->type);
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));

			}
		}
	}


	unordered_map<string, double> par_left, par_right;
	// Запишем координаты поверхностей (на самом деле центров граней)

	int test_i = 121;
	out.write(reinterpret_cast<const char*>(&test_i), sizeof(int));

	// TS
	if (true)
	{
		size = this->Gran_TS.size() * 7;
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));
		std::ofstream outfile("TS_interpol.txt");

		for (const auto& gr : this->Gran_TS)
		{
			double aa = gr->center[0][0];
			double bb = gr->center[0][1];
			double cc = gr->center[0][2];

			double r_1, the_1, phi_1;

			r_1 = sqrt(aa * aa + bb * bb + cc * cc);
			the_1 = acos(aa / r_1);
			phi_1 = polar_angle(bb, cc);

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия phi
			phi_1 = phi_1 + 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия phi
			phi_1 = phi_1 - 2 * const_pi;
			phi_1 = phi_1 - 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия - theta
			phi_1 = phi_1 + 2 * const_pi;  // вернул
			the_1 = -the_1;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия - theta  phi
			phi_1 = phi_1 + 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия - theta  phi
			phi_1 = phi_1 - 2 * const_pi; // вернул
			phi_1 = phi_1 - 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия - theta
			phi_1 = phi_1 + 2 * const_pi; // вернул
			the_1 = -the_1;
			the_1 = 2 * const_pi - the_1;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}
		}

		outfile.close();
	}

	test_i = 122;
	out.write(reinterpret_cast<const char*>(&test_i), sizeof(int));

	// HP
	if (true)
	{
		// сначала радиальная запись HP
		size = 0;
		for (const auto& gr : this->Gran_HP)
		{
			double aa = gr->center[0][0];
			if (aa >= -5.0) size++;
		}
		size = size * 4 + 101 + 120;
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));
		std::ofstream outfile("HP_interpol.txt");

		for (const auto& gr : this->Gran_HP)
		{
			double aa = gr->center[0][0];
			double bb = gr->center[0][1];
			double cc = gr->center[0][2];

			if (aa < -5.0) continue;

			double r_1, the_1, phi_1;

			r_1 = sqrt(aa * aa + bb * bb + cc * cc);
			the_1 = acos(aa / r_1);
			phi_1 = polar_angle(bb, cc);

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);
			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия phi
			phi_1 = phi_1 + 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия phi
			phi_1 = phi_1 - 2 * const_pi;
			phi_1 = phi_1 - 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия - theta
			phi_1 = phi_1 + 2 * const_pi;
			the_1 = -the_1;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

		}

		// добавим нулевую точку 101 раз
		if(true)
		{
			double time;
			bool bb = false;
			Gran* G = nullptr;
			Eigen::Vector3d orig, Vel;

			for (const auto& gr : this->Gran_HP)
			{
				orig << 0.0, 0.0, 0.0;
				Vel << 1.0, 0.0, 0.0;
				bool aa = gr->Luch_crossing(orig, Vel, time);
				if (aa == true && time > 0)
				{
					G = gr;
					bb = true;
					break;
				}
			}

			if (bb == false)
			{
				cout << "Error 7657367y5h4w56346  " << endl;
				exit(-1);
			}

			Eigen::Vector3d CC;
			CC = orig + Vel * time;

			double r_1, the_1, phi_1;

			r_1 = CC.norm();
			the_1 = acos(CC[0] / r_1);
			phi_1 = -0.06;

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			Gran* gr = G;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);
			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));

			for(size_t ii = 0; ii < 100; ii++)
			{
				phi_1 = -0.05 + ii * (2.2 * const_pi) / 100.0;
				outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

				out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));
			}

		}

		// добавим нулевую точку 101 раз
		if (true)
		{
			double time;
			bool bb = false;
			Gran* G = nullptr;
			Eigen::Vector3d orig, Vel;
			for (size_t ii = 0; ii < 120; ii++)
			{
				double phi_1 = -const_pi / 9.176784 + ii * (2 * const_pi / (120 - 10));
				for (const auto& gr : this->Gran_HP)
				{
					orig << 0.0, 0.0, 0.0;
					Vel << 0.0, cos(phi_1), sin(phi_1);
					bool aa = gr->Luch_crossing(orig, Vel, time);
					if (aa == true && time > 0)
					{
						G = gr;
						bb = true;
						break;
					}
				}

				if (bb == false)
				{
					cout << "Error 7657367y5h4w56346  " << endl;
					exit(-1);
				}

				Eigen::Vector3d CC;
				CC = orig + Vel * time;

				double r_1, the_1;

				r_1 = CC.norm();
				the_1 = acos(CC[0] / r_1);

				outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

				Gran* gr = G;

				out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));
			}

		}

	
		outfile.close();

		// Теперь цилиндрическая запись HP

		size_t size_x = 100;
		size_t size_phi = 120;

		out.write(reinterpret_cast<const char*>(&size_x), sizeof(size_x));
		out.write(reinterpret_cast<const char*>(&size_phi), sizeof(size_phi));
		outfile.open("HP_2_interpol.txt");

		for (size_t i = 0; i < size_phi; ++i)
		{
			for (size_t j = 0; j < size_x; ++j)
			{
				double phi_ = -const_pi / 9.176784 + i * (2 * const_pi / (size_phi - 10));
				double x_ = this->geo->L6 * 0.9999 + j * (-this->geo->L6 + 5) / size_x;
				double time;
				bool bb = false;
				Gran* G = nullptr;
				Eigen::Vector3d orig, Vel;

				for (const auto& gr : this->Gran_HP)
				{
					orig << x_, 0.0, 0.0;
					Vel << 0.0, cos(phi_), sin(phi_);
					bool aa = gr->Luch_crossing(orig, Vel, time);
					if (aa == true && time > 0)
					{
						G = gr;
						bb = true;
						break;
					}
				}

				if (bb == false)
				{
					cout << "Error 8765ugeugg346  " << endl;
					cout << x_ << " " << phi_ << endl;
					exit(-1);
				}

				Eigen::Vector3d CC;
				CC = orig + Vel * time;

				double r_1 = sqrt(kv(CC[1]) + kv(CC[2]));

				outfile << x_ << " " << phi_ << " " << r_1 << endl;

				Gran* gr = G;

				out.write(reinterpret_cast<const char*>(&x_), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));

			}
		}
		outfile.close();
	}

	test_i = 123;
	out.write(reinterpret_cast<const char*>(&test_i), sizeof(int));

	// BS
	if (true)
	{
		// сначала радиальная запись BS
		size = 0;
		for (const auto& gr : this->Gran_BS)
		{
			double aa = gr->center[0][0];
			size++;
		}
		size = size * 4 + 101 + 120;
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));
		std::ofstream outfile("BS_interpol.txt");

		for (const auto& gr : this->Gran_BS)
		{
			double aa = gr->center[0][0];
			double bb = gr->center[0][1];
			double cc = gr->center[0][2];

			double r_1, the_1, phi_1;

			r_1 = sqrt(aa * aa + bb * bb + cc * cc);
			the_1 = acos(aa / r_1);
			phi_1 = polar_angle(bb, cc);

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);
			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия phi
			phi_1 = phi_1 + 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия phi
			phi_1 = phi_1 - 2 * const_pi;
			phi_1 = phi_1 - 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия - theta
			phi_1 = phi_1 + 2 * const_pi;
			the_1 = -the_1;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

		}

		// добавим нулевую точку 101 раз
		if (true)
		{
			double time;
			bool bb = false;
			Gran* G = nullptr;
			Eigen::Vector3d orig, Vel;

			for (const auto& gr : this->Gran_BS)
			{
				orig << 0.0, 0.0, 0.0;
				Vel << 1.0, 0.0, 0.0;
				bool aa = gr->Luch_crossing(orig, Vel, time);
				if (aa == true && time > 0)
				{
					G = gr;
					bb = true;
					break;
				}
			}

			if (bb == false)
			{
				cout << "Error 7657367y5h4w56346  " << endl;
				exit(-1);
			}

			Eigen::Vector3d CC;
			CC = orig + Vel * time;

			double r_1, the_1, phi_1;

			r_1 = CC.norm();
			the_1 = acos(CC[0] / r_1);
			phi_1 = -0.06;

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			Gran* gr = G;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);
			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));

			for (size_t ii = 0; ii < 100; ii++)
			{
				phi_1 = -0.05 + ii * (2.2 * const_pi) / 100.0;
				outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

				out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));
			}

		}

		// добавим доп точку 120 раз
		if (true)
		{
			double time;
			bool bb = false;
			Gran* G = nullptr;
			Eigen::Vector3d orig, Vel;
			for (size_t ii = 0; ii < 120; ii++)
			{
				double phi_1 = -const_pi / 9.176784 + ii * (2 * const_pi / (120 - 10));
				for (const auto& gr : this->Gran_BS)
				{
					orig << 0.0001, 0.0, 0.0;
					Vel << 0.0, cos(phi_1), sin(phi_1);
					bool aa = gr->Luch_crossing(orig, Vel, time);
					if (aa == true && time > 0)
					{
						G = gr;
						bb = true;
						break;
					}
				}

				if (bb == false)
				{
					cout << "Error yurthrtegdk76  " << endl;
					exit(-1);
				}

				Eigen::Vector3d CC;
				CC = orig + Vel * time;

				double r_1, the_1;

				r_1 = CC.norm();
				the_1 = acos(CC[0] / r_1);

				outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

				Gran* gr = G;

				out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));
			}

		}


		outfile.close();
	}

	test_i = 124;
	out.write(reinterpret_cast<const char*>(&test_i), sizeof(int));


	for (size_t i = 0; i < 999; i++)
	{
		bool aa = false;
		out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
	}

	out.close();
}

void Setka::Save_for_interpolate_one_zone_only(string filename, Type_cell ZONA)
{
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		cout << "Error 097564537  Can not  open file to writing: " + filename << endl;
		exit(-1);
	}

	bool razriv = false;
	// Записываем есть ли особенная интерполяция на разрывах (или всё сплошным образом)
	out.write(reinterpret_cast<const char*>(&razriv), sizeof(bool));

	// Записываем до какого расстояния слева выделяется HP
	out.write(reinterpret_cast<const char*>(&this->geo->L6), sizeof(double));

	// Добавляем ещё переменные для вывода  "BB/8pi"
	if (true)
	{
		this->phys_param->param_names.push_back("BB/8pi");
		for (const auto& Cel : this->All_Cell)
		{
			Cel->parameters[0]["BB/8pi"] = kvv(Cel->parameters[0]["Bx"],
				Cel->parameters[0]["By"], Cel->parameters[0]["Bz"]) / (8.0 * const_pi);
		}
		this->Cell_Center->parameters[0]["BB/8pi"] = kvv(this->Cell_Center->parameters[0]["Bx"],
			this->Cell_Center->parameters[0]["By"], this->Cell_Center->parameters[0]["Bz"]) / (8.0 * const_pi);
	}


	// Записываем количество строк
	size_t size = this->phys_param->param_names.size() + 1;
	out.write(reinterpret_cast<const char*>(&size), sizeof(size));

	// Записываем каждую строку
	for (const auto& str : this->phys_param->param_names) {
		// Сначала записываем длину строки
		size_t str_size = str.size();
		out.write(reinterpret_cast<const char*>(&str_size), sizeof(str_size));
		// Затем саму строку
		out.write(str.data(), str_size);
	}

	cout << "All parameters (send): " << endl;
	for (const auto& i : this->phys_param->param_names)
	{
		cout << i << "  ";
	}
	cout << endl;

	// Добавляем геометрическую зону
	if (true)
	{
		string str = "zone_geo";
		size_t str_size = str.size();
		out.write(reinterpret_cast<const char*>(&str_size), sizeof(str_size));
		out.write(str.data(), str_size);
	}

	// Считаем сколько дополнительных ячеек будет на внешней границе
	unsigned int gr_b = 0;

	if (true)
	{
		if (true)
		{
			unsigned int N_cell = 0;
			for (const auto& Cel : this->All_Cell)
			{
				Cel->is_need = static_cast<short int>(Cel->type);

				if (Cel->type == ZONA) N_cell++;
			}

			if (ZONA == Type_cell::Zone_1) N_cell++;

			// Записываем количество ячеек
			size = N_cell; // + Центр?
			out.write(reinterpret_cast<const char*>(&size), sizeof(size));

			for (const auto& Cel : this->All_Cell)
			{
				if (Cel->type != ZONA) continue;
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


				double zzz = static_cast<double>(Cel->is_need);
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
			}

			// Записываем центральную точку
			if (ZONA == Type_cell::Zone_1)
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

				//cout << "================  " << this->Cell_Center->parameters[0]["rho"] << endl;

				double zzz = 1.0;
				out.write(reinterpret_cast<const char*>(&zzz), sizeof(zzz));
				//cout << "================  " << zzz << endl;
			}
		}
	}

	unordered_map<string, double> par_left, par_right;
	// Запишем координаты поверхностей (на самом деле центров граней)

	int test_i = 121;
	out.write(reinterpret_cast<const char*>(&test_i), sizeof(int));

	// TS
	if (true)
	{
		size = this->Gran_TS.size() * 7;
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));
		std::ofstream outfile("TS_interpol.txt");

		for (const auto& gr : this->Gran_TS)
		{
			double aa = gr->center[0][0];
			double bb = gr->center[0][1];
			double cc = gr->center[0][2];

			double r_1, the_1, phi_1;

			r_1 = sqrt(aa * aa + bb * bb + cc * cc);
			the_1 = acos(aa / r_1);
			phi_1 = polar_angle(bb, cc);

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия phi
			phi_1 = phi_1 + 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия phi
			phi_1 = phi_1 - 2 * const_pi;
			phi_1 = phi_1 - 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия - theta
			phi_1 = phi_1 + 2 * const_pi;  // вернул
			the_1 = -the_1;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия - theta  phi
			phi_1 = phi_1 + 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия - theta  phi
			phi_1 = phi_1 - 2 * const_pi; // вернул
			phi_1 = phi_1 - 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}

			// Симметрия - theta
			phi_1 = phi_1 + 2 * const_pi; // вернул
			the_1 = -the_1;
			the_1 = 2 * const_pi - the_1;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			for (const auto& str : this->phys_param->param_names)
			{
				out.write(reinterpret_cast<const char*>(&par_left[str]), sizeof(cc));
				out.write(reinterpret_cast<const char*>(&par_right[str]), sizeof(cc));
			}
		}

		outfile.close();
	}

	test_i = 122;
	out.write(reinterpret_cast<const char*>(&test_i), sizeof(int));

	// HP
	if (true)
	{
		// сначала радиальная запись HP
		size = 0;
		for (const auto& gr : this->Gran_HP)
		{
			double aa = gr->center[0][0];
			if (aa >= -5.0) size++;
		}
		size = size * 4 + 101 + 120;
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));
		std::ofstream outfile("HP_interpol.txt");

		for (const auto& gr : this->Gran_HP)
		{
			double aa = gr->center[0][0];
			double bb = gr->center[0][1];
			double cc = gr->center[0][2];

			if (aa < -5.0) continue;

			double r_1, the_1, phi_1;

			r_1 = sqrt(aa * aa + bb * bb + cc * cc);
			the_1 = acos(aa / r_1);
			phi_1 = polar_angle(bb, cc);

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);
			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия phi
			phi_1 = phi_1 + 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия phi
			phi_1 = phi_1 - 2 * const_pi;
			phi_1 = phi_1 - 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия - theta
			phi_1 = phi_1 + 2 * const_pi;
			the_1 = -the_1;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

		}

		// добавим нулевую точку 101 раз
		if (true)
		{
			double time;
			bool bb = false;
			Gran* G = nullptr;
			Eigen::Vector3d orig, Vel;

			for (const auto& gr : this->Gran_HP)
			{
				orig << 0.0, 0.0, 0.0;
				Vel << 1.0, 0.0, 0.0;
				bool aa = gr->Luch_crossing(orig, Vel, time);
				if (aa == true && time > 0)
				{
					G = gr;
					bb = true;
					break;
				}
			}

			if (bb == false)
			{
				cout << "Error 7657367y5h4w56346  " << endl;
				exit(-1);
			}

			Eigen::Vector3d CC;
			CC = orig + Vel * time;

			double r_1, the_1, phi_1;

			r_1 = CC.norm();
			the_1 = acos(CC[0] / r_1);
			phi_1 = -0.06;

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			Gran* gr = G;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);
			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));

			for (size_t ii = 0; ii < 100; ii++)
			{
				phi_1 = -0.05 + ii * (2.2 * const_pi) / 100.0;
				outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

				out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));
			}

		}

		// добавим нулевую точку 101 раз
		if (true)
		{
			double time;
			bool bb = false;
			Gran* G = nullptr;
			Eigen::Vector3d orig, Vel;
			for (size_t ii = 0; ii < 120; ii++)
			{
				double phi_1 = -const_pi / 9.176784 + ii * (2 * const_pi / (120 - 10));
				for (const auto& gr : this->Gran_HP)
				{
					orig << 0.0, 0.0, 0.0;
					Vel << 0.0, cos(phi_1), sin(phi_1);
					bool aa = gr->Luch_crossing(orig, Vel, time);
					if (aa == true && time > 0)
					{
						G = gr;
						bb = true;
						break;
					}
				}

				if (bb == false)
				{
					cout << "Error 7657367y5h4w56346  " << endl;
					exit(-1);
				}

				Eigen::Vector3d CC;
				CC = orig + Vel * time;

				double r_1, the_1;

				r_1 = CC.norm();
				the_1 = acos(CC[0] / r_1);

				outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

				Gran* gr = G;

				out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));
			}

		}


		outfile.close();

		// Теперь цилиндрическая запись HP

		size_t size_x = 100;
		size_t size_phi = 120;

		out.write(reinterpret_cast<const char*>(&size_x), sizeof(size_x));
		out.write(reinterpret_cast<const char*>(&size_phi), sizeof(size_phi));
		outfile.open("HP_2_interpol.txt");

		for (size_t i = 0; i < size_phi; ++i)
		{
			for (size_t j = 0; j < size_x; ++j)
			{
				double phi_ = -const_pi / 9.176784 + i * (2 * const_pi / (size_phi - 10));
				double x_ = this->geo->L6 * 0.9999 + j * (-this->geo->L6 + 5) / size_x;
				double time;
				bool bb = false;
				Gran* G = nullptr;
				Eigen::Vector3d orig, Vel;

				for (const auto& gr : this->Gran_HP)
				{
					orig << x_, 0.0, 0.0;
					Vel << 0.0, cos(phi_), sin(phi_);
					bool aa = gr->Luch_crossing(orig, Vel, time);
					if (aa == true && time > 0)
					{
						G = gr;
						bb = true;
						break;
					}
				}

				if (bb == false)
				{
					cout << "Error 8765ugeugg346  " << endl;
					cout << x_ << " " << phi_ << endl;
					exit(-1);
				}

				Eigen::Vector3d CC;
				CC = orig + Vel * time;

				double r_1 = sqrt(kv(CC[1]) + kv(CC[2]));

				outfile << x_ << " " << phi_ << " " << r_1 << endl;

				Gran* gr = G;

				out.write(reinterpret_cast<const char*>(&x_), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));

			}
		}
		outfile.close();
	}

	test_i = 123;
	out.write(reinterpret_cast<const char*>(&test_i), sizeof(int));

	// BS
	if (true)
	{
		// сначала радиальная запись BS
		size = 0;
		for (const auto& gr : this->Gran_BS)
		{
			double aa = gr->center[0][0];
			size++;
		}
		size = size * 4 + 101 + 120;
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));
		std::ofstream outfile("BS_interpol.txt");

		for (const auto& gr : this->Gran_BS)
		{
			double aa = gr->center[0][0];
			double bb = gr->center[0][1];
			double cc = gr->center[0][2];

			double r_1, the_1, phi_1;

			r_1 = sqrt(aa * aa + bb * bb + cc * cc);
			the_1 = acos(aa / r_1);
			phi_1 = polar_angle(bb, cc);

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);
			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия phi
			phi_1 = phi_1 + 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия phi
			phi_1 = phi_1 - 2 * const_pi;
			phi_1 = phi_1 - 2 * const_pi;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

			// Симметрия - theta
			phi_1 = phi_1 + 2 * const_pi;
			the_1 = -the_1;
			out.write(reinterpret_cast<const char*>(&the_1), sizeof(bb));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(aa));
			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(cc));

			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(cc));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(cc));

		}

		// добавим нулевую точку 101 раз
		if (true)
		{
			double time;
			bool bb = false;
			Gran* G = nullptr;
			Eigen::Vector3d orig, Vel;

			for (const auto& gr : this->Gran_BS)
			{
				orig << 0.0, 0.0, 0.0;
				Vel << 1.0, 0.0, 0.0;
				bool aa = gr->Luch_crossing(orig, Vel, time);
				if (aa == true && time > 0)
				{
					G = gr;
					bb = true;
					break;
				}
			}

			if (bb == false)
			{
				cout << "Error 7657367y5h4w56346  " << endl;
				exit(-1);
			}

			Eigen::Vector3d CC;
			CC = orig + Vel * time;

			double r_1, the_1, phi_1;

			r_1 = CC.norm();
			the_1 = acos(CC[0] / r_1);
			phi_1 = -0.06;

			outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

			Gran* gr = G;

			out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
			out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
			out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

			out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

			this->Snos_on_Gran(gr, par_left, par_right, 0, true);
			out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
			out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));

			for (size_t ii = 0; ii < 100; ii++)
			{
				phi_1 = -0.05 + ii * (2.2 * const_pi) / 100.0;
				outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

				out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));
			}

		}

		// добавим доп точку 120 раз
		if (true)
		{
			double time;
			bool bb = false;
			Gran* G = nullptr;
			Eigen::Vector3d orig, Vel;
			for (size_t ii = 0; ii < 120; ii++)
			{
				double phi_1 = -const_pi / 9.176784 + ii * (2 * const_pi / (120 - 10));
				for (const auto& gr : this->Gran_BS)
				{
					orig << 0.0001, 0.0, 0.0;
					Vel << 0.0, cos(phi_1), sin(phi_1);
					bool aa = gr->Luch_crossing(orig, Vel, time);
					if (aa == true && time > 0)
					{
						G = gr;
						bb = true;
						break;
					}
				}

				if (bb == false)
				{
					cout << "Error yurthrtegdk76  " << endl;
					exit(-1);
				}

				Eigen::Vector3d CC;
				CC = orig + Vel * time;

				double r_1, the_1;

				r_1 = CC.norm();
				the_1 = acos(CC[0] / r_1);

				outfile << the_1 << " " << phi_1 << " " << r_1 << endl;

				Gran* gr = G;

				out.write(reinterpret_cast<const char*>(&the_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&phi_1), sizeof(double));
				out.write(reinterpret_cast<const char*>(&r_1), sizeof(double));

				out.write(reinterpret_cast<const char*>(&gr->normal[0][0]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][1]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&gr->normal[0][2]), sizeof(double));

				this->Snos_on_Gran(gr, par_left, par_right, 0, true);
				out.write(reinterpret_cast<const char*>(&par_left["rho"]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&par_right["rho"]), sizeof(double));
			}

		}


		outfile.close();
	}

	test_i = 124;
	out.write(reinterpret_cast<const char*>(&test_i), sizeof(int));


	for (size_t i = 0; i < 999; i++)
	{
		bool aa = false;
		out.write(reinterpret_cast<const char*>(&aa), sizeof(aa));
	}

	out.close();
}


void Setka::Save_cell_pui_parameters(string filename)
{
	std::ofstream out("PUI_" + filename, std::ios::binary);
	if (!out) {
		cout << "Error 097564537  Can not open file to writing: " + filename << endl;
		exit(-1);
	}

	double d;
	for (auto& i : this->All_Cell)
	{
		short int zone = this->determ_zone(i, 0);

		if (zone != 2)
		{
			d = i->parameters[0]["rho_Pui_1"] / i->parameters[0]["rho"];
			out.write(reinterpret_cast<const char*>(&d), sizeof(double));
			d = i->parameters[0]["p_Pui_1"] / i->parameters[0]["p"];
			out.write(reinterpret_cast<const char*>(&d), sizeof(double));
		}
		else
		{
			d = i->parameters[0]["rho_Pui_1"] / i->parameters[0]["rho"];
			out.write(reinterpret_cast<const char*>(&d), sizeof(double));
			d = i->parameters[0]["rho_Pui_2"] / i->parameters[0]["rho"];
			out.write(reinterpret_cast<const char*>(&d), sizeof(double));
			d = i->parameters[0]["p_Pui_1"] / i->parameters[0]["p"];
			out.write(reinterpret_cast<const char*>(&d), sizeof(double));
			d = i->parameters[0]["p_Pui_2"] / i->parameters[0]["p"];
			out.write(reinterpret_cast<const char*>(&d), sizeof(double));
		}
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

	out.close();
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
				//cout << "Error 8765509090" << endl;
				//cout << pair << endl;
				//cout << "-------------------" << endl;
				i->parameters[0][pair] = 0.0;
			}

			// записываем значение
			out.write(reinterpret_cast<const char*>(&i->parameters[0][pair]), sizeof(double));
		}
	}

	bool bb;

	// Записываем для будущих считываний
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

	cout << "Start Download_cell_MK_parameters  for  setka " << this->name << endl;

	vector<string> param_file;
	//cout << "MK_parameters: ";
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

		//cout << key << " ";
		param_file.push_back(key);
	}
	//cout << endl;


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

void Setka::PereInterpolate(string filename, bool move, bool MK_only)
{
	cout << "PereInterpolate: step 1/4" << endl;
	Interpol SS = Interpol(filename);

	this->PereInterpolate(&SS, move, MK_only);
}

void Setka::PereInterpolate(Interpol* SS, bool move, bool MK_only)
{
	cout << "PereInterpolate: step 2/4" << endl;
	// Сначала двигаем все поверхности
	if (move)
	{
		double x, y, z, r, rr, the, phi, R_BS;
		std::unordered_map<string, double> param2;

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


					SS->Get_TS(x, y, z, param2);
					rr = param2["r"];
					//rr = Surf->Get_TS(phi, the);

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

					SS->Get_HP(x, y, z, param2);
					rr = param2["r"];
					//rr = Surf->Get_HP(phi, the, 0);

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

					if (x < 0.01) x = 0.01;
					SS->Get_BS(x, y, z, param2);
					rr = param2["r"];
					//rr = Surf->Get_BS(phi, the);

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


				SS->Get_TS(x, y, z, param2);
				rr = param2["r"];
				//rr = Surf->Get_TS(phi, the);

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

				SS->Get_HP(x, y, z, param2);
				rr = param2["r"];
				//rr = Surf->Get_HP(phi, the, 0);

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

				SS->Get_BS(x, y, z, param2);
				rr = param2["r"];
				//rr = Surf->Get_BS(phi, the);

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

					SS->Get_TS(x, y, z, param2);
					rr = param2["r"];
					//rr = Surf->Get_TS(phi, the);

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

					SS->Get_HP(x, y, z, param2);
					rr = param2["r"];
					//rr = Surf->Get_HP(phi, x, 1);
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


					SS->Get_TS(x, y, z, param2);
					rr = param2["r"];
					//rr = Surf->Get_TS(phi, the);

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


				SS->Get_TS(x, y, z, param2);
				rr = param2["r"];
				//rr = Surf->Get_TS(phi, the);

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
				for (auto& j : i)
				{
					x = j->Yzels_opor[1]->coord[0][0];
					y = j->Yzels_opor[1]->coord[0][1];
					z = j->Yzels_opor[1]->coord[0][2];
					r = sqrt(kvv(0.0, y, z));

					SS->Get_HP(x, y, z, param2);
					rr = param2["r"];
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

					SS->Get_HP(x, y, z, param2);
					rr = param2["r"];
					//rr = Surf->Get_HP(phi, x, 1);

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

	}

	cout << "PereInterpolate: step 3/4" << endl;

	// Теперь переинтерполируем все значения в ячейках

	std::unordered_map<string, double> param;
	std::array<Cell_handle, 6> prev_cell;
	std::array<Cell_handle, 6> next_cell;
	for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();

	for (auto& cel : this->All_Cell)
	{
		double x, y, z;
		x = cel->center[0][0];
		y = cel->center[0][1];
		z = cel->center[0][2];
		bool bb;
		bb = SS->Get_param(x, y, z, param, prev_cell, next_cell);     // Интерполируем переменные
		while(bb == false)
		{	
			x = x * 0.99;
			y = y * 0.99;
			z = z * 0.99;
			bb = SS->Get_param(x, y, z, param, prev_cell, next_cell);     // Интерполируем переменные
			if (bb == false)
			{
				cout << "ifuhueiruoweiohgf83478  Warning! pereinterpolate x,y,z " << x << " " << y << " " << z << endl;
			}
		}

		for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i]; // Обновляем предыдущую ячейку
		for (const auto& [key, value] : param)
		{
			if (MK_only == false || std::find(this->phys_param->MK_param.begin(), this->phys_param->MK_param.end(), key) != this->phys_param->MK_param.end())
			{
				cel->parameters[0][key] = value;
			}
		}

		cel->parameters[1] = cel->parameters[0];
	}
	cout << "PereInterpolate: step 4/4" << endl;
}

void Setka::Culc_rotors_in_cell(void)
{
	// Это не обычный ротер, а делённый на B^2
	cout << "Start: Culc_rotor_in_cell" << endl;

	this->phys_param->param_names.push_back("rotB/b2_x");
	this->phys_param->param_names.push_back("rotB/b2_y");
	this->phys_param->param_names.push_back("rotB/b2_z");
	// Добавили переменную для интерполяции
	unsigned int k1 = 0;

#pragma omp parallel for schedule(dynamic)
	//for (auto& cell : this->All_Cell)
	for (size_t idx = 0; idx < this->All_Cell.size(); ++idx) 
	{
		auto& cell = this->All_Cell[idx];
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
							this->Snos_on_Gran(gr_, par_left, par_right, 0, true);

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

void Setka::Culc_usual_rotors_in_cell(void)
{
	cout << "Start: Culc_rotor_in_cell" << endl;

	this->phys_param->param_names.push_back("rotB_x");
	this->phys_param->param_names.push_back("rotB_y");
	this->phys_param->param_names.push_back("rotB_z");
	// Добавили переменную для интерполяции
	unsigned int k1 = 0;

#pragma omp parallel for schedule(dynamic)
	//for (auto& cell : this->All_Cell)
	for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
	{
		auto& cell = this->All_Cell[idx];
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
							this->Snos_on_Gran(gr_, par_left, par_right, 0, true);

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

		cell->parameters[0]["rotB_x"] = vvv[0];
		cell->parameters[0]["rotB_y"] = vvv[1];
		cell->parameters[0]["rotB_z"] = vvv[2];

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
	this->phys_param->param_names.push_back("DIVk_x");
	this->phys_param->param_names.push_back("DIVk_y");
	this->phys_param->param_names.push_back("DIVk_z");
	//this->phys_param->param_names.push_back("div_et");
	//this->phys_param->param_names.push_back("gradK_x");
	//this->phys_param->param_names.push_back("gradK_y");
	//this->phys_param->param_names.push_back("gradK_z");
	//this->phys_param->param_names.push_back("et_x");
	//this->phys_param->param_names.push_back("et_y");
	//this->phys_param->param_names.push_back("et_z");
	// Добавили переменную для интерполяции

#pragma omp parallel for
	for (size_t i_step = 0; i_step < this->All_Cell.size(); i_step++)
	{
		unordered_map<string, double> par_left, par_right;
		Eigen::Vector3d normal, Vel;
		Eigen::Vector3d eB;

		auto& cell = this->All_Cell[i_step];
		double divV = 0.0;
		double DIVk_x = 0.0;
		double DIVk_y = 0.0;
		double DIVk_z = 0.0;


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
				this->Snos_on_Gran(gr, par_left, par_right, 0, true);

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
			Eigen::Vector3d D1, D2, D3;
			if (norm > 1e-10)
			{
				double k_paral = 1.0 / norm;
				double k_perpend = 0.05 * k_paral;

				D1 << k_perpend + (k_paral - k_perpend) * eB[0] * eB[0] / kv(norm),
					(k_paral - k_perpend)* eB[0] * eB[1] / kv(norm),
					(k_paral - k_perpend)* eB[0] * eB[2] / kv(norm);
				D2 << (k_paral - k_perpend) * eB[1] * eB[0] / kv(norm),
					k_perpend + (k_paral - k_perpend) * eB[1] * eB[1] / kv(norm),
					(k_paral - k_perpend)* eB[1] * eB[2] / kv(norm);
				D3 << (k_paral - k_perpend) * eB[2] * eB[0] / kv(norm),
					(k_paral - k_perpend)* eB[2] * eB[1] / kv(norm),
					k_perpend + (k_paral - k_perpend) * eB[2] * eB[2] / kv(norm);
			}
			else
			{
				D1 << 0.0, 0.0, 0.0;
				D2 << 0.0, 0.0, 0.0;
				D3 << 0.0, 0.0, 0.0;
			}

			divV += Vel.dot(normal) * gr->area[0];
			DIVk_x += D1.dot(normal) * gr->area[0];
			DIVk_y += D2.dot(normal) * gr->area[0];
			DIVk_z += D3.dot(normal) * gr->area[0];

		}

		cell->parameters[0]["divV"] = divV/cell->volume[0];
		cell->parameters[0]["DIVk_x"] = DIVk_x / cell->volume[0];
		cell->parameters[0]["DIVk_y"] = DIVk_y / cell->volume[0];
		cell->parameters[0]["DIVk_z"] = DIVk_z / cell->volume[0];
	}

	cout << "End: Culc_divergence_in_cell" << endl;
}


void Setka::Culc_gradient_in_cell(void)
{
	cout << "Start: Culc_grad_in_cell" << endl;

	this->phys_param->param_names.push_back("gradBB_x");
	this->phys_param->param_names.push_back("gradBB_y");
	this->phys_param->param_names.push_back("gradBB_z");
	//this->phys_param->param_names.push_back("div_et");
	//this->phys_param->param_names.push_back("gradK_x");
	//this->phys_param->param_names.push_back("gradK_y");
	//this->phys_param->param_names.push_back("gradK_z");
	//this->phys_param->param_names.push_back("et_x");
	//this->phys_param->param_names.push_back("et_y");
	//this->phys_param->param_names.push_back("et_z");
	// Добавили переменную для интерполяции

#pragma omp parallel for
	for (size_t i_step = 0; i_step < this->All_Cell.size(); i_step++)
	{
		unordered_map<string, double> par_left, par_right;
		Eigen::Vector3d normal;
		Eigen::Vector3d eB;

		auto& cell = this->All_Cell[i_step];
		double gr_x = 0.0;
		double gr_y = 0.0;
		double gr_z = 0.0;


		// пробегаемся по всем граням 
		for (const auto& gr : cell->grans)
		{
			eB << 0.0, 0.0, 0.0;


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
				this->Snos_on_Gran(gr, par_left, par_right, 0, true);

				if (gr->type2 == Type_Gran_surf::Us)
				{
					eB << (par_left["Bx"] + par_right["Bx"]) / 2.0,
						(par_left["By"] + par_right["By"]) / 2.0,
						(par_left["Bz"] + par_right["Bz"]) / 2.0;
				}
				else
				{
					if (gr->cells[0] == cell)
					{
						eB << par_left["Bx"],
							par_left["By"],
							par_left["Bz"];
					}
					else
					{
						eB << par_right["Bx"],
							par_right["By"],
							par_right["Bz"];
					}
				}
			}
			else
			{
				eB << cell->parameters[0]["Bx"], cell->parameters[0]["By"], cell->parameters[0]["Bz"];
			}

			double BB = kv(eB.norm());

			gr_x += BB * normal[0] * gr->area[0];
			gr_y += BB * normal[1] * gr->area[0];
			gr_z += BB * normal[2] * gr->area[0];

		}

		cell->parameters[0]["gradBB_x"] = gr_x / cell->volume[0];
		cell->parameters[0]["gradBB_y"] = gr_y / cell->volume[0];
		cell->parameters[0]["gradBB_z"] = gr_z / cell->volume[0];
	}

	cout << "End: Culc_grad_in_cell" << endl;
}
