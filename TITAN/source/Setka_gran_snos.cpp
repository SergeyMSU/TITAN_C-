#include "Setka.h"

void Setka::Snos_on_Gran(Gran* gr, unordered_map<string, double>& par_left,
	unordered_map<string, double>& par_right, short int now)
{
	// Для граничных граней не надо делать ТВД
	if (gr->type != Type_Gran::Us)
	{
		if (gr->type == Type_Gran::Inner_Hard || gr->type == Type_Gran::Outer_Hard)
		{
			// Нужно ли считать плазменные поля
			if (this->phys_param->culc_plasma == true)
			{
				par_left["rho"] = gr->parameters["rho"];
				par_left["n_He"] = gr->parameters["n_He"];
				par_left["Q"] = gr->parameters["Q"];
				par_left["p"] = gr->parameters["p"];
				par_left["Vx"] = gr->parameters["Vx"];
				par_left["Vy"] = gr->parameters["Vy"];
				par_left["Vz"] = gr->parameters["Vz"];
				par_left["Bx"] = gr->parameters["Bx"];
				par_left["By"] = gr->parameters["By"];
				par_left["Bz"] = gr->parameters["Bz"];

				par_right["rho"] = par_left["rho"];
				par_right["n_He"] = par_left["n_He"];
				par_right["Q"] = par_left["Q"];
				par_right["p"] = par_left["p"];
				par_right["Vx"] = par_left["Vx"];
				par_right["Vy"] = par_left["Vy"];
				par_right["Vz"] = par_left["Vz"];
				par_right["Bx"] = par_left["Bx"];
				par_right["By"] = par_left["By"];
				par_right["Bz"] = par_left["Bz"];
			}

			for (auto& nam : this->phys_param->H_name)
			{
				if (gr->type == Type_Gran::Outer_Hard || (nam == "_H1"))
				{
					par_left["rho" + nam] = gr->parameters["rho" + nam];
					par_left["p" + nam] = gr->parameters["p" + nam];
					par_left["Vx" + nam] = gr->parameters["Vx" + nam];
					par_left["Vy" + nam] = gr->parameters["Vy" + nam];
					par_left["Vz" + nam] = gr->parameters["Vz" + nam];

					par_right["rho" + nam] = par_left["rho" + nam];
					par_right["p" + nam] = par_left["p" + nam];
					par_right["Vx" + nam] = par_left["Vx" + nam];
					par_right["Vy" + nam] = par_left["Vy" + nam];
					par_right["Vz" + nam] = par_left["Vz" + nam];
				}
				/*if (gr->type == Type_Gran::Inner_Hard && (nam == "_H2"))
				{
					par_left["rho" + nam] = 0.00001;
					par_left["p" + nam] = 0.00001;
					par_left["Vx" + nam] = 2.0 * fabs(this->phys_param->Velosity_inf);
					par_left["Vy" + nam] = 2.0 * fabs(this->phys_param->Velosity_inf);
					par_left["Vz" + nam] = 2.0 * fabs(this->phys_param->Velosity_inf);

					par_right["rho" + nam] = par_left["rho" + nam];
					par_right["p" + nam] = par_left["p" + nam];
					par_right["Vx" + nam] = par_left["Vx" + nam];
					par_right["Vy" + nam] = par_left["Vy" + nam];
					par_right["Vz" + nam] = par_left["Vz" + nam];
				}*/
				else // Для  H2  H3  H4    и  Inner_Hard  -  делаем через центральную фиктивную ячейку
				{
					auto A = gr->cells[0];
					auto B = this->Cell_Center;
					par_left["rho" + nam] = A->parameters[now]["rho" + nam];
					par_left["p" + nam] = A->parameters[now]["p" + nam];
					par_left["Vx" + nam] = A->parameters[now]["Vx" + nam];
					par_left["Vy" + nam] = A->parameters[now]["Vy" + nam];
					par_left["Vz" + nam] = A->parameters[now]["Vz" + nam];

					par_right["rho" + nam] = B->parameters[now]["rho" + nam];
					par_right["p" + nam] = B->parameters[now]["p" + nam];
					par_right["Vx" + nam] = B->parameters[now]["Vx" + nam];
					par_right["Vy" + nam] = B->parameters[now]["Vy" + nam];
					par_right["Vz" + nam] = B->parameters[now]["Vz" + nam];
				}
			}
		}
		else if (gr->type == Type_Gran::Outer_Soft)
		{
			auto C = gr->cells[0];
			if (this->phys_param->culc_plasma == true)
			{
				par_left["rho"] = C->parameters[now]["rho"];
				par_left["n_He"] = C->parameters[now]["n_He"];
				par_left["Q"] = C->parameters[now]["Q"];
				par_left["p"] = C->parameters[now]["p"];
				par_left["Vx"] = C->parameters[now]["Vx"];
				par_left["Vy"] = C->parameters[now]["Vy"];
				par_left["Vz"] = C->parameters[now]["Vz"];
				par_left["Bx"] = C->parameters[now]["Bx"];
				par_left["By"] = C->parameters[now]["By"];
				par_left["Bz"] = C->parameters[now]["Bz"];

				par_right["rho"] = par_left["rho"];
				par_right["n_He"] = par_left["n_He"];
				par_right["Q"] = par_left["Q"];
				par_right["p"] = par_left["p"];
				par_right["Vx"] = par_left["Vx"];
				par_right["Vy"] = par_left["Vy"];
				par_right["Vz"] = par_left["Vz"];
				par_right["Bx"] = par_left["Bx"];
				par_right["By"] = par_left["By"];
				par_right["Bz"] = par_left["Bz"];

				// Запрещаем затекание жидкости через мягкие граничные условия
				if (par_left["Vx"] * gr->normal[now][0] + 
					par_left["Vy"] * gr->normal[now][1] +
					par_left["Vz"] * gr->normal[now][2] < 0.0)
				{
					par_left["Vx"] = 0.0;
					par_left["Vy"] = 0.0;
					par_left["Vz"] = 0.0;
				}

				if (gr->normal[now][0] < -0.9) // Для задней границы
				{
					// Отсос
					if (par_left["Vx"] > this->phys_param->Velosity_inf / 7.0)
					{
						par_left["Vx"] = this->phys_param->Velosity_inf / 5.0;
					}
				}

			}

			for (auto& nam : this->phys_param->H_name)
			{
				par_left["rho" + nam] = C->parameters[now]["rho" + nam];
				par_left["p" + nam] = C->parameters[now]["p" + nam];
				par_left["Vx" + nam] = C->parameters[now]["Vx" + nam];
				par_left["Vy" + nam] = C->parameters[now]["Vy" + nam];
				par_left["Vz" + nam] = C->parameters[now]["Vz" + nam];

				par_right["rho" + nam] = par_left["rho" + nam];
				par_right["p" + nam] = par_left["p" + nam];
				par_right["Vx" + nam] = par_left["Vx" + nam];
				par_right["Vy" + nam] = par_left["Vy" + nam];
				par_right["Vz" + nam] = par_left["Vz" + nam];
			}
		}
		else
		{
			cout << "Error 5473190855" << endl;
			exit(-1);
		}
	}
	else
	{
		// Делаем ли ТВД?
		if (this->phys_param->TVD == true)
		{
			Eigen::Vector3d Ac, Bc, AAc, BBc, G, vec;
			auto A = gr->cells[0];
			auto B = gr->cells[1];
			auto AA = gr->cells_TVD[0];
			auto BB = gr->cells_TVD[1];
			if (AA == nullptr || BB == nullptr) goto a1; // Для таких не надо делать TVD
		
			G << gr->center[now][0], gr->center[now][1], gr->center[now][2];
			Ac << A->center[now][0], A->center[now][1], A->center[now][2];
			Bc << B->center[now][0], B->center[now][1], B->center[now][2];
			AAc << AA->center[now][0], AA->center[now][1], AA->center[now][2];
			BBc << BB->center[now][0], BB->center[now][1], BB->center[now][2];
			double d1, d2, dd1, dd2;
			vec = Ac - G;
			d1 = vec.norm();
			vec = Bc - G;
			d2 = vec.norm();
			vec = Ac - AAc;
			dd1 = vec.norm();
			vec = Bc - BBc;
			dd2 = vec.norm();

			for (auto& nam: this->phys_param->param_names)
			{
				par_left[nam]  = linear(-dd1 - d1, AA->parameters[now][nam],
					-d1, A->parameters[now][nam], 
					d2, B->parameters[now][nam], 0.0);
				par_right[nam] = linear(-d1, A->parameters[now][nam],
					d2, B->parameters[now][nam], 
					d2 + dd2, BB->parameters[now][nam], 0.0);

				if (nam == "rho" || nam == "p" || nam == "Vx" || nam == "Vy"
					|| nam == "Vz")
				{
					for (auto& nam2 : this->phys_param->H_name)
					{
						par_left[nam + nam2] = linear(-dd1 - d1, AA->parameters[now][nam + nam2],
							-d1, A->parameters[now][nam + nam2],
							d2, B->parameters[now][nam + nam2], 0.0);
						par_right[nam + nam2] = linear(-d1, A->parameters[now][nam + nam2],
							d2, B->parameters[now][nam + nam2],
							d2 + dd2, BB->parameters[now][nam + nam2], 0.0);
					}
				}
			}

			if (par_left["rho"] < 0.0000001) par_left["rho"] = A->parameters[now]["rho"];
			if (par_left["Q"] < 0.0000001) par_left["Q"] = A->parameters[now]["Q"];
			if (par_left["n_He"] < 0.0000001) par_left["n_He"] = A->parameters[now]["n_He"];
			if (par_left["p"] < 0.0000001) par_left["p"] = A->parameters[now]["p"];
			if (par_right["rho"] < 0.0000001) par_right["rho"] = B->parameters[now]["rho"];
			if (par_right["Q"] < 0.0000001) par_right["Q"] = B->parameters[now]["Q"];
			if (par_right["n_He"] < 0.0000001) par_right["n_He"] = B->parameters[now]["n_He"];
			if (par_right["p"] < 0.0000001) par_right["p"] = B->parameters[now]["p"];
			
			for (auto& nam2 : this->phys_param->H_name)
			{
				if (par_left["rho" + nam2] < 0.0000001) par_left["rho" + nam2] = A->parameters[now]["rho" + nam2];
				if (par_left["p" + nam2] < 0.0000001) par_left["p" + nam2] = A->parameters[now]["p" + nam2];
				if (par_right["rho" + nam2] < 0.0000001) par_right["rho" + nam2] = B->parameters[now]["rho" + nam2];
				if (par_right["p" + nam2] < 0.0000001) par_right["p" + nam2] = B->parameters[now]["p" + nam2];
			}
		}
		else
		{
			a1:
			auto A = gr->cells[0];
			auto B = gr->cells[1];

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

			if (this->phys_param->culc_plasma == true)
			{
				par_left["rho"] = A->parameters[now]["rho"];
				par_left["Q"] = A->parameters[now]["Q"];
				par_left["n_He"] = A->parameters[now]["n_He"];
				par_left["p"] = A->parameters[now]["p"];
				par_left["Vx"] = A->parameters[now]["Vx"];
				par_left["Vy"] = A->parameters[now]["Vy"];
				par_left["Vz"] = A->parameters[now]["Vz"];
				par_left["Bx"] = A->parameters[now]["Bx"];
				par_left["By"] = A->parameters[now]["By"];
				par_left["Bz"] = A->parameters[now]["Bz"];

				par_right["rho"] = B->parameters[now]["rho"];
				par_right["Q"] = B->parameters[now]["Q"];
				par_right["n_He"] = B->parameters[now]["n_He"];
				par_right["p"] = B->parameters[now]["p"];
				par_right["Vx"] = B->parameters[now]["Vx"];
				par_right["Vy"] = B->parameters[now]["Vy"];
				par_right["Vz"] = B->parameters[now]["Vz"];
				par_right["Bx"] = B->parameters[now]["Bx"];
				par_right["By"] = B->parameters[now]["By"];
				par_right["Bz"] = B->parameters[now]["Bz"];
			}

			for (auto& nam : this->phys_param->H_name)
			{
				par_left["rho" + nam] = A->parameters[now]["rho" + nam];
				par_left["p" + nam] = A->parameters[now]["p" + nam];
				par_left["Vx" + nam] = A->parameters[now]["Vx" + nam];
				par_left["Vy" + nam] = A->parameters[now]["Vy" + nam];
				par_left["Vz" + nam] = A->parameters[now]["Vz" + nam];

				par_right["rho" + nam] = B->parameters[now]["rho" + nam];
				par_right["p" + nam] = B->parameters[now]["p" + nam];
				par_right["Vx" + nam] = B->parameters[now]["Vx" + nam];
				par_right["Vy" + nam] = B->parameters[now]["Vy" + nam];
				par_right["Vz" + nam] = B->parameters[now]["Vz" + nam];
			}


		}
	}
}