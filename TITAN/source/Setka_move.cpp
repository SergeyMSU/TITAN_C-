#include "Setka.h"

void Setka::Culc_Velocity_surface(short int now, const double& time, short int metod)
{
	double dsr, dsc, dsl;
	
	int now2 = (now + 1) % 2;

	for (auto& yz : this->All_Yzel)
	{
		yz->num_velocity = 0;
		yz->velocity[0] = 0.0;
		yz->velocity[1] = 0.0;
		yz->velocity[2] = 0.0;
	}

	// ������� �������� TS
#pragma omp parallel for private(dsr, dsc, dsl)
	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		auto A = gr->cells[0];
		auto B = gr->cells[1];

		std::vector<double> qqq, qqq1, qqq2;
		qqq.resize(8);
		qqq1.resize(8);
		qqq2.resize(8);
		std::vector<double> konvect_left, konvect_right, konvect;
		PrintOptions Option = PrintOptions{};


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

		double w = 0.0;

		this->phys_param->chlld(metod, gr->normal[now][0], gr->normal[now][1],
			gr->normal[now][2],
			w, qqq1, qqq2, qqq, false, 3,
			konvect_left, konvect_right, konvect, dsr, dsc, dsl,
			Option);

		for (auto& yz : gr->yzels)
		{
#pragma omp critical (s1) 
			{
				yz->velocity[0] += 0.1 * dsl * gr->normal[now][0];
			}
#pragma omp critical (s2) 
			{
				yz->velocity[1] += 0.1 * dsl * gr->normal[now][1];
			}
#pragma omp critical (s3) 
			{
				yz->velocity[2] += 0.1 * dsl * gr->normal[now][2];
			}
#pragma omp critical (s4) 
			{
				yz->num_velocity++;
			}
		}

	}

	// ������� �������� HP
	if (this->phys_param->move_HP == true)
	{
#pragma omp parallel for private(dsr, dsc, dsl)
		for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
		{
			auto gr = this->Gran_HP[i_step];
			auto A = gr->cells[0];
			auto B = gr->cells[1];

			std::vector<double> qqq, qqq1, qqq2;
			qqq.resize(8);
			qqq1.resize(8);
			qqq2.resize(8);
			std::vector<double> konvect_left, konvect_right, konvect;
			PrintOptions Option = PrintOptions{};


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

			double w = 0.0;

			this->phys_param->chlld(metod, gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 3,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option);

			for (auto& yz : gr->yzels)
			{
#pragma omp critical (s1) 
				{
					yz->velocity[0] += 0.1 * dsc * gr->normal[now][0];
				}
#pragma omp critical (s2) 
				{
					yz->velocity[1] += 0.1 * dsc * gr->normal[now][1];
				}
#pragma omp critical (s3) 
				{
					yz->velocity[2] += 0.1 * dsc * gr->normal[now][2];
				}
#pragma omp critical (s4) 
				{
					yz->num_velocity++;
				}
			}

		}
	}

	// ��������� ����� ���������� ����� �� �����������
	for (auto& yz : this->All_Yzel)
	{
		if (yz->num_velocity > 0)
		{
			yz->velocity[0] /= yz->num_velocity;
			yz->velocity[1] /= yz->num_velocity;
			yz->velocity[2] /= yz->num_velocity;
		}

		// ��������� ����� ���������� ����
		if (yz->type == Type_yzel::TS)
		{
			Eigen::Vector3d A, B;
			A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];
			Eigen::Vector3d V;
			V << yz->velocity[0], yz->velocity[1], yz->velocity[2];
			B = A;
			B.normalize();
			A = A + V.dot(B) * time * B;
			yz->coord[now2][0] = A[0];
			yz->coord[now2][1] = A[1];
			yz->coord[now2][2] = A[2];
		}
		else if (yz->type == Type_yzel::HP && this->phys_param->move_HP == true)
		{
			Eigen::Vector3d A, B;
			A << yz->coord[now][0], yz->coord[now][1], yz->coord[now][2];
			Eigen::Vector3d V;
			V << yz->velocity[0], yz->velocity[1], yz->velocity[2];
			if (A(0) >= 0.0)
			{
				B = A;
			}
			else
			{
				B << 0.0, A(1), A(2);
			}
			B.normalize();
			A = A + V.dot(B) * time * B;
			yz->coord[now2][0] = A[0];
			yz->coord[now2][1] = A[1];
			yz->coord[now2][2] = A[2];
		}
	}


	// ����������� ����������� (��� ��������� ����� ���������� ����� �� �����������)



	// ��������� ���� �� HP (������������ �����) ���� ���������
	if (this->phys_param->move_HP == true)
		//if(false)
	{
		short int NN = this->D_Luch[0].size() - 1;
		for (auto& L : this->D_Luch)
		{
			double h1 = norm2(0.0, L[this->geo->N4]->Yzels_opor[1]->coord[now2][1],
				L[this->geo->N4]->Yzels_opor[1]->coord[now2][2]);
			double h2 = norm2(0.0, L[NN]->Yzels_opor[1]->coord[now2][1],
				L[NN]->Yzels_opor[1]->coord[now2][2]);

			for (short int i = this->geo->N4; i <= NN; i++)
			{
				double h = h1 + (i - this->geo->N4) * (h2 - h1) / (NN - this->geo->N4);
				auto yz = L[i]->Yzels_opor[1];
				double hh = norm2(0.0, yz->coord[now2][1], yz->coord[now2][2]);
				if (hh < 0.0000001 || h < 0.0000001 || std::isnan(hh) || std::isnan(h) ||
					std::fpclassify(h) == FP_SUBNORMAL || std::fpclassify(hh) == FP_SUBNORMAL)
					// || fabs(1.0 - h / hh) > 0.0001)
				{
					cout << "Error 8563529613" << endl;
					whach(h);
					whach(hh);
					whach(h1);
					whach(h2);
				}


				yz->coord[now2][1] = yz->coord[now2][1] * h / hh;
				yz->coord[now2][2] = yz->coord[now2][2] * h / hh;
			}
		}
	}
	
}