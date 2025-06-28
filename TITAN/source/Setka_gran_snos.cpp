#include "Setka.h"

void Setka::Snos_on_Gran(Gran* gr, unordered_map<string, double>& par_left,
	unordered_map<string, double>& par_right, short int now)
{
	// Для граничных граней не надо делать ТВД
	if (gr->type != Type_Gran::Us)
	{
		if (gr->type == Type_Gran::Inner_Hard || gr->type == Type_Gran::Outer_Hard)
		{
			auto C = gr->cells[0];
			// Нужно ли считать плазменные поля
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

				par_right["rho"] = gr->parameters["rho"];
				par_right["n_He"] = gr->parameters["n_He"];
				par_right["Q"] = gr->parameters["Q"];
				par_right["p"] = gr->parameters["p"];
				par_right["Vx"] = gr->parameters["Vx"];
				par_right["Vy"] = gr->parameters["Vy"];
				par_right["Vz"] = gr->parameters["Vz"];
				par_right["Bx"] = gr->parameters["Bx"];
				par_right["By"] = gr->parameters["By"];
				par_right["Bz"] = gr->parameters["Bz"];

			}

			for (auto& nam : this->phys_param->H_name)
			{
				// Для четвёртого сорта жёсткие гран условия
				if (gr->type == Type_Gran::Outer_Hard && nam == "_H4")
				{
					par_left["rho" + nam] = C->parameters[now]["rho" + nam];
					par_left["p" + nam] = C->parameters[now]["p" + nam];
					par_left["Vx" + nam] = C->parameters[now]["Vx" + nam];
					par_left["Vy" + nam] = C->parameters[now]["Vy" + nam];
					par_left["Vz" + nam] = C->parameters[now]["Vz" + nam];

					par_right["rho" + nam] = gr->parameters["rho" + nam];
					par_right["p" + nam] = gr->parameters["p" + nam];
					par_right["Vx" + nam] = gr->parameters["Vx" + nam];
					par_right["Vy" + nam] = gr->parameters["Vy" + nam];
					par_right["Vz" + nam] = gr->parameters["Vz" + nam];
				}
				else if (gr->type == Type_Gran::Outer_Hard || nam == "_H1")
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
				else if (gr->type == Type_Gran::Inner_Hard)// Для  H2  H3  H4    и  Inner_Hard  -  делаем через центральную фиктивную ячейку
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
				else
				{
					cout << "Error 8653496143 " << endl;
					cout << nam << endl;
					exit(-1);
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

				// Запрещаем затекание жидкости через мягкие граничные условия
				if (par_left["Vx"] * gr->normal[now][0] + 
					par_left["Vy"] * gr->normal[now][1] +
					par_left["Vz"] * gr->normal[now][2] < 0.0)
				{
					par_left["Vx"] = 0.1 * gr->normal[now][0];
					par_left["Vy"] = 0.1 * gr->normal[now][1];
					par_left["Vz"] = 0.1 * gr->normal[now][2];
				}

				if (gr->normal[now][0] < -0.9) // Для задней границы
				{
					// Отсос
					if (par_left["Vx"] > this->phys_param->Velosity_inf / 7.0)
					{
						par_left["Vx"] = this->phys_param->Velosity_inf / 5.0;
					}
				}

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
				par_left["rho" + nam] = C->parameters[now]["rho" + nam];
				par_left["p" + nam] = C->parameters[now]["p" + nam];
				par_left["Vx" + nam] = C->parameters[now]["Vx" + nam];
				par_left["Vy" + nam] = C->parameters[now]["Vy" + nam];
				par_left["Vz" + nam] = C->parameters[now]["Vz" + nam];

				// Запрещаем затекание жидкости через мягкие граничные условия
				if (par_left["Vx" + nam] * gr->normal[now][0] +
					par_left["Vy" + nam] * gr->normal[now][1] +
					par_left["Vz" + nam] * gr->normal[now][2] < 0.0)
				{
					par_left["Vx" + nam] = 0.1 * gr->normal[now][0];
					par_left["Vy" + nam] = 0.1 * gr->normal[now][1];
					par_left["Vz" + nam] = 0.1 * gr->normal[now][2];
				}



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
			// AA  -  A  -|-  B  -  BB
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

			// AA  -  A  -|-  B  -  BB
			// r3,    r1, rr, r2,   r4
			double r3, r1, rr, r2, r4;

			if (AA->type == Type_cell::Zone_1 && A->type == Type_cell::Zone_1 &&
				B->type == Type_cell::Zone_1 && BB->type == Type_cell::Zone_1)
			{
				Eigen::Vector3d VAA, VA, VB, VBB, Vleft, Vright;
				r3 = AAc.norm();
				r1 = Ac.norm();
				rr = G.norm();
				r2 = Bc.norm();
				r4 = BBc.norm();

				// Преобразуем нужные векторные величины
				if (true)
				{
					std::array<string, 3> V1;
					std::array<string, 3> V2;
					std::array<string, 3> V3;

					V1[0] = "Vx";
					V2[0] = "Vy";
					V3[0] = "Vz";

				    V1[1] = "Bx";
					V2[1] = "By";
					V3[1] = "Bz";

					V1[2] = "Vx_H1";
					V2[2] = "Vy_H1";
					V3[2] = "Vz_H1";

					for (short int ik = 0; ik < 3; ik++)
					{
						/*cout << AA->parameters[now][V1[ik]] << " " <<
							AA->parameters[now][V2[ik]] << " " <<
							AA->parameters[now][V3[ik]] << endl;*/

						// Переводим скорости в сферическую с.к.
						spherical_skorost(AAc[2], AAc[0], AAc[1],
							AA->parameters[now][V3[ik]], AA->parameters[now][V1[ik]],
							AA->parameters[now][V2[ik]], VAA[0], VAA[1], VAA[2]);

						spherical_skorost(Ac[2], Ac[0], Ac[1],
							A->parameters[now][V3[ik]], A->parameters[now][V1[ik]],
							A->parameters[now][V2[ik]], VA[0], VA[1], VA[2]);

						spherical_skorost(Bc[2], Bc[0], Bc[1],
							B->parameters[now][V3[ik]], B->parameters[now][V1[ik]],
							B->parameters[now][V2[ik]], VB[0], VB[1], VB[2]);

						spherical_skorost(BBc[2], BBc[0], BBc[1],
							BB->parameters[now][V3[ik]], BB->parameters[now][V1[ik]],
							BB->parameters[now][V2[ik]], VBB[0], VBB[1], VBB[2]);

						// Интерполируем скорости (в сферической с.к.)
						for (short int i = 0; i < 3; i++)
						{
							Vleft[i] = linear(-dd1 - d1, VAA[i], -d1, VA[i], d2, VB[i], 0.0);
							Vright[i] = linear(-d1, VA[i], d2, VB[i], d2 + dd2, VBB[i], 0.0);
						}

						dekard_skorost(G[2], G[0], G[1],
							Vleft[0], Vleft[1], Vleft[2],
							par_left[V3[ik]], par_left[V1[ik]], par_left[V2[ik]]);

						dekard_skorost(G[2], G[0], G[1],
							Vright[0], Vright[1], Vright[2],
							par_right[V3[ik]], par_right[V1[ik]], par_right[V2[ik]]);

						/*cout << " = " << par_left[V1[ik]] << " " <<
							par_left[V2[ik]] << " " <<
							par_left[V3[ik]] << endl;*/
					}
				}


			
				for (auto& nam : this->phys_param->param_names)
				{
					if (nam == "rho" || nam == "n_He" || nam == "Q")
					{
						par_left[nam] = linear(-dd1 - d1, AA->parameters[now][nam] * kv(r3),
							-d1, A->parameters[now][nam] * kv(r1),
							d2, B->parameters[now][nam] * kv(r2), 0.0)/kv(rr);
						par_right[nam] = linear(-d1, A->parameters[now][nam] * kv(r1),
							d2, B->parameters[now][nam] * kv(r2),
							d2 + dd2, BB->parameters[now][nam] * kv(r4), 0.0) / kv(rr);
					}
					else if(nam == "p")
					{
						par_left[nam] = linear(-dd1 - d1, AA->parameters[now][nam] * kvg(r3),
							-d1, A->parameters[now][nam] * kvg(r1),
							d2, B->parameters[now][nam] * kvg(r2), 0.0) / kvg(rr);
						par_right[nam] = linear(-d1, A->parameters[now][nam] * kvg(r1),
							d2, B->parameters[now][nam] * kvg(r2),
							d2 + dd2, BB->parameters[now][nam] * kvg(r4), 0.0) / kvg(rr);
					}
					else if (nam == "Vx" || nam == "Vy" || nam == "Vz" 
						|| nam == "Bx" || nam == "By" || nam == "Bz")
					{
						// PASS
						;
					}
					else
					{
						par_left[nam] = linear(-dd1 - d1, AA->parameters[now][nam],
							-d1, A->parameters[now][nam],
							d2, B->parameters[now][nam], 0.0);
						par_right[nam] = linear(-d1, A->parameters[now][nam],
							d2, B->parameters[now][nam],
							d2 + dd2, BB->parameters[now][nam], 0.0);
					}

					if (nam == "rho" || nam == "p" || nam == "Vx" || nam == "Vy"
						|| nam == "Vz")
					{
						for (auto& nam2 : this->phys_param->H_name)
						{
							if (nam2 == "_H1" && (nam == "Vx" || nam == "Vy" || nam == "Vz")) continue;
							par_left[nam + nam2] = linear(-dd1 - d1, AA->parameters[now][nam + nam2],
								-d1, A->parameters[now][nam + nam2],
								d2, B->parameters[now][nam + nam2], 0.0);
							par_right[nam + nam2] = linear(-d1, A->parameters[now][nam + nam2],
								d2, B->parameters[now][nam + nam2],
								d2 + dd2, BB->parameters[now][nam + nam2], 0.0);
						}
					}
				}

			}
			else
			{
				for (auto& nam : this->phys_param->param_names)
				{
					par_left[nam] = linear(-dd1 - d1, AA->parameters[now][nam],
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
				Eigen::Vector3d VAA, VA, VB, VBB, Vleft, Vright;

				if (A->type == Type_cell::Zone_1)
				{
					Eigen::Vector3d Ac, G;
					Ac << A->center[now][0], A->center[now][1], A->center[now][2];
					G << gr->center[now][0], gr->center[now][1], gr->center[now][2];
					double r1 = Ac.norm();
					double rr = G.norm();

					// Преобразуем нужные векторные величины
					if (true)
					{
						std::array<string, 3> V1;
						std::array<string, 3> V2;
						std::array<string, 3> V3;

						V1[0] = "Vx";
						V2[0] = "Vy";
						V3[0] = "Vz";

						V1[1] = "Bx";
						V2[1] = "By";
						V3[1] = "Bz";

						V1[2] = "Vx_H1";
						V2[2] = "Vy_H1";
						V3[2] = "Vz_H1";

						for (short int ik = 0; ik < 3; ik++)
						{
							// Переводим скорости в сферическую с.к.
							spherical_skorost(Ac[2], Ac[0], Ac[1],
								A->parameters[now][V3[ik]], A->parameters[now][V1[ik]],
								A->parameters[now][V2[ik]], VA[0], VA[1], VA[2]);

							dekard_skorost(G[2], G[0], G[1],
								VA[0], VA[1], VA[2],
								par_left[V3[ik]], par_left[V1[ik]], par_left[V2[ik]]);
						}
					}

					par_left["rho"] = A->parameters[now]["rho"] * kv(r1) / kv(rr);
					par_left["Q"] = A->parameters[now]["Q"] * kv(r1) / kv(rr);
					par_left["n_He"] = A->parameters[now]["n_He"] * kv(r1) / kv(rr);
					par_left["p"] = A->parameters[now]["p"] * kvg(r1) / kvg(rr);

					for (auto& nam : this->phys_param->H_name)
					{
						par_left["rho" + nam] = A->parameters[now]["rho" + nam];
						par_left["p" + nam] = A->parameters[now]["p" + nam];
						if (nam == "_H1") continue;
						par_left["Vx" + nam] = A->parameters[now]["Vx" + nam];
						par_left["Vy" + nam] = A->parameters[now]["Vy" + nam];
						par_left["Vz" + nam] = A->parameters[now]["Vz" + nam];
					}
				}
				else
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

					for (auto& nam : this->phys_param->H_name)
					{
						par_left["rho" + nam] = A->parameters[now]["rho" + nam];
						par_left["p" + nam] = A->parameters[now]["p" + nam];
						par_left["Vx" + nam] = A->parameters[now]["Vx" + nam];
						par_left["Vy" + nam] = A->parameters[now]["Vy" + nam];
						par_left["Vz" + nam] = A->parameters[now]["Vz" + nam];
					}
				}

				if (B->type == Type_cell::Zone_1)
				{
					Eigen::Vector3d Ac, G;
					Ac << B->center[now][0], B->center[now][1], B->center[now][2];
					G << gr->center[now][0], gr->center[now][1], gr->center[now][2];
					double r1 = Ac.norm();
					double rr = G.norm();

					// Преобразуем нужные векторные величины
					if (true)
					{
						std::array<string, 3> V1;
						std::array<string, 3> V2;
						std::array<string, 3> V3;

						V1[0] = "Vx";
						V2[0] = "Vy";
						V3[0] = "Vz";

						V1[1] = "Bx";
						V2[1] = "By";
						V3[1] = "Bz";

						V1[2] = "Vx_H1";
						V2[2] = "Vy_H1";
						V3[2] = "Vz_H1";

						for (short int ik = 0; ik < 3; ik++)
						{
							// Переводим скорости в сферическую с.к.
							spherical_skorost(Ac[2], Ac[0], Ac[1],
								B->parameters[now][V3[ik]], B->parameters[now][V1[ik]],
								B->parameters[now][V2[ik]], VA[0], VA[1], VA[2]);

							dekard_skorost(G[2], G[0], G[1],
								VA[0], VA[1], VA[2],
								par_right[V3[ik]], par_right[V1[ik]], par_right[V2[ik]]);
						}
					}

					par_right["rho"] = B->parameters[now]["rho"] * kv(r1) / kv(rr);
					par_right["Q"] = B->parameters[now]["Q"] * kv(r1) / kv(rr);
					par_right["n_He"] = B->parameters[now]["n_He"] * kv(r1) / kv(rr);
					par_right["p"] = B->parameters[now]["p"] * kvg(r1) / kvg(rr);
					
					for (auto& nam : this->phys_param->H_name)
					{
						par_right["rho" + nam] = B->parameters[now]["rho" + nam];
						par_right["p" + nam] = B->parameters[now]["p" + nam];
						if (nam == "_H1") continue;
						par_right["Vx" + nam] = B->parameters[now]["Vx" + nam];
						par_right["Vy" + nam] = B->parameters[now]["Vy" + nam];
						par_right["Vz" + nam] = B->parameters[now]["Vz" + nam];
					}
				}
				else
				{
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

					for (auto& nam : this->phys_param->H_name)
					{
						par_right["rho" + nam] = B->parameters[now]["rho" + nam];
						par_right["p" + nam] = B->parameters[now]["p" + nam];
						par_right["Vx" + nam] = B->parameters[now]["Vx" + nam];
						par_right["Vy" + nam] = B->parameters[now]["Vy" + nam];
						par_right["Vz" + nam] = B->parameters[now]["Vz" + nam];
					}
				}
			}
		}
	}
}