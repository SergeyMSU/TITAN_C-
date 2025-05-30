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

	// Считаем движение TS
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
			yz->velocity[0] += 0.1 * dsl * gr->normal[now][0];
			yz->velocity[1] += 0.1 * dsl * gr->normal[now][1];
			yz->velocity[2] += 0.1 * dsl * gr->normal[now][2];
			yz->num_velocity++;
		}

	}


	for (auto& yz : this->All_Yzel)
	{
		if (yz->num_velocity > 0)
		{
			yz->velocity[0] /= yz->num_velocity;
			yz->velocity[1] /= yz->num_velocity;
			yz->velocity[2] /= yz->num_velocity;
		}

		// Вычисляем новую координату узла
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
	}


}