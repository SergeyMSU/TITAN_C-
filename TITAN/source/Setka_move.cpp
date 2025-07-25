#include "Setka.h"

// ������ ��������������� ����������� HP
#define  macros1() p[0] = AA->coord[now][0]; \
p[1] = norm2(0.0, AA->coord[now][1], AA->coord[now][2]);\
\
p1[0] = AA1->coord[now][0];\
p1[1] = norm2(0.0, AA1->coord[now][1], AA1->coord[now][2]);\
\
p11[0] = AA11->coord[now][0];\
p11[1] = norm2(0.0, AA11->coord[now][1], AA11->coord[now][2]);\
\
p3[0] = AA3->coord[now][0];\
p3[1] = norm2(0.0, AA3->coord[now][1], AA3->coord[now][2]);\
\
p33[0] = AA33->coord[now][0];\
p33[1] = norm2(0.0, AA33->coord[now][1], AA33->coord[now][2]);\
\
solveQuadraticEquation(p3[0], p3[1],\
	p1[0], p1[1], p11[0], p11[1], a, b, c);\
r1 = a * kv(p[0]) + b * p[0] + c;\
\
solveQuadraticEquation(p33[0], p33[1],\
	p3[0], p3[1], p1[0], p1[1], a, b, c);\
r2 = a * kv(p[0]) + b * p[0] + c;\
\
p2[0] = polar_angle(AA2->coord[now][1], AA2->coord[now][2]);\
p2[1] = norm2(0.0, AA2->coord[now][1], AA2->coord[now][2]);\
\
p[0] = polar_angle(AA->coord[now][1], AA->coord[now][2]);\
\
p4[0] = polar_angle(AA4->coord[now][1], AA4->coord[now][2]);\
p4[1] = norm2(0.0, AA4->coord[now][1], AA4->coord[now][2]);\
\
p22[0] = polar_angle(AA22->coord[now][1], AA22->coord[now][2]);\
p22[1] = norm2(0.0, AA22->coord[now][1], AA22->coord[now][2]);\
\
p44[0] = polar_angle(AA44->coord[now][1], AA44->coord[now][2]);\
p44[1] = norm2(0.0, AA44->coord[now][1], AA44->coord[now][2]);\
\
solveQuadraticEquation(p2[0], p2[1],\
	p4[0], p4[1], p44[0], p44[1], a, b, c);\
r3 = a * kv(p[0]) + b * p[0] + c;\
\
solveQuadraticEquation(p4[0], p4[1],\
	p2[0], p2[1], p22[0], p22[1], a, b, c);\
r4 = a * kv(p[0]) + b * p[0] + c;\
\
A << AA->coord[now][0], AA->coord[now][1], AA->coord[now][2];\
\
B = A;\
B[1] *= r1 / p[1];\
B[2] *= r1 / p[1];\
\
V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_along * this->phys_param->sglag_HP_k * (B - A) / time;\
\
AA->mut.lock();\
AA->velocity[0] += V[0];\
AA->velocity[1] += V[1];\
AA->velocity[2] += V[2];\
AA->num_velocity++;\
AA->mut.unlock();\
\
B = A;\
B[1] *= r2 / p[1];\
B[2] *= r2 / p[1];\
\
V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_along * this->phys_param->sglag_HP_k * (B - A) / time;\
\
AA->mut.lock();\
AA->velocity[0] += V[0];\
AA->velocity[1] += V[1];\
AA->velocity[2] += V[2];\
AA->num_velocity++;\
AA->mut.unlock();\
\
B = A;\
B[1] *= r3 / p[1];\
B[2] *= r3 / p[1];\
\
V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_angle * this->phys_param->sglag_HP_k * (B - A) / time;\
\
AA->mut.lock();\
AA->velocity[0] += V[0];\
AA->velocity[1] += V[1];\
AA->velocity[2] += V[2];\
AA->num_velocity++;\
AA->mut.unlock();\
\
B = A;\
B[1] *= r4 / p[1];\
B[2] *= r4 / p[1];\
\
V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_angle * this->phys_param->sglag_HP_k * (B - A) / time;\
\
AA->mut.lock();\
AA->velocity[0] += V[0];\
AA->velocity[1] += V[1];\
AA->velocity[2] += V[2];\
AA->num_velocity++;\
AA->mut.unlock();


#define  macros2(p, AA) p[0] = polar_angle(AA->coord[now][0],\
	norm2(0.0, AA->coord[now][1], AA->coord[now][2]));\
	if (p[0] > const_pi) p[0] = p[0] - 2 * const_pi;\
	p[1] = norm2(AA->coord[now][0], AA->coord[now][1], AA->coord[now][2]);



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
	if (this->phys_param->move_TS == true)
	{
		#pragma omp parallel for private(dsr, dsc, dsl)
		for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
		{
			auto gr = this->Gran_TS[i_step];
			auto A = gr->cells[0];
			auto B = gr->cells[1];


			unordered_map<string, double> par_left;
			unordered_map<string, double> par_right;

			Snos_on_Gran(gr, par_left, par_right, now);

			std::vector<double> qqq, qqq1, qqq2;
			qqq.resize(8);
			qqq1.resize(8);
			qqq2.resize(8);
			std::vector<double> konvect_left, konvect_right, konvect;
			PrintOptions Option = PrintOptions{};
			Option.fluid = "plasma_TS";

			//A->parameters[now]
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

			double w = 0.0;


			this->phys_param->chlld(gr->Get_method(), gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 3,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option);

			for (auto& yz : gr->yzels)
			{
				yz->mut.lock();
				yz->velocity[0] += this->phys_param->velocity_TS * dsl * gr->normal[now][0];
				yz->velocity[1] += this->phys_param->velocity_TS * dsl * gr->normal[now][1];
				yz->velocity[2] += this->phys_param->velocity_TS * dsl * gr->normal[now][2];
				yz->num_velocity++;
				yz->mut.unlock();
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
			short int metod_ = gr->Get_method();

			if (gr->type2 != Type_Gran_surf::HP) cout << "Error 7823467276345679264978234" << endl;


			if (this->phys_param->null_bn_on_HP == true)
			{
				// �������� bn � �������
				Eigen::Vector3d nnn;
				nnn << gr->normal[now][0], gr->normal[now][1], gr->normal[now][2];
				double dtpr = scalarProductFast(nnn[0], nnn[1], nnn[2],
					A->parameters[now]["Bx"], A->parameters[now]["By"],
					A->parameters[now]["Bz"]);
				A->parameters[now]["Bx"] -= dtpr * nnn[0];
				A->parameters[now]["By"] -= dtpr * nnn[1];
				A->parameters[now]["Bz"] -= dtpr * nnn[2];

				dtpr = scalarProductFast(nnn[0], nnn[1], nnn[2],
					B->parameters[now]["Bx"], B->parameters[now]["By"],
					B->parameters[now]["Bz"]);
				B->parameters[now]["Bx"] -= dtpr * nnn[0];
				B->parameters[now]["By"] -= dtpr * nnn[1];
				B->parameters[now]["Bz"] -= dtpr * nnn[2];
			}

			std::vector<double> qqq, qqq1, qqq2;
			qqq.resize(8);
			qqq1.resize(8);
			qqq2.resize(8);
			std::vector<double> konvect_left, konvect_right, konvect;
			PrintOptions Option = PrintOptions{};
			Option.fluid = "plasma_HP";

			// �������� ��������� ��� �������� ��������
			unordered_map<string, double> par_left, par_right;
			this->Snos_on_Gran(gr, par_left, par_right, now);

			if (false) // ��������� � ������� �����
			{
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
			else
			{
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

				// ��� �������� �������� ���� ����� �������� Bn
				if (this->phys_param->null_bn_on_HP == true)
				{
					// �������� bn � �������
					Eigen::Vector3d nnn;
					nnn << gr->normal[now][0], gr->normal[now][1], gr->normal[now][2];
					double dtpr = scalarProductFast(nnn[0], nnn[1], nnn[2],
						qqq1[5], qqq1[6],
						qqq1[7]);
					qqq1[5] -= dtpr * nnn[0];
					qqq1[6] -= dtpr * nnn[1];
					qqq1[7] -= dtpr * nnn[2];

					dtpr = scalarProductFast(nnn[0], nnn[1], nnn[2],
						qqq2[5], qqq2[6], qqq2[7]);
					qqq2[5] -= dtpr * nnn[0];
					qqq2[6] -= dtpr * nnn[1];
					qqq2[7] -= dtpr * nnn[2];
				}
			}

			double w = 0.0;

			// ���������� ��������� �������� � �������
			// � ������� ��������� ����
			if (this->phys_param->bn_in_p_on_HP == true)
			{
				metod_ = 2;

				qqq1[4] += kvv(qqq1[5], qqq1[6], qqq1[7]) / (8.0 * const_pi);
				qqq1[5] = qqq1[6] = qqq1[7] = 0.0;

				qqq2[4] += kvv(qqq2[5], qqq2[6], qqq2[7]) / (8.0 * const_pi);
				qqq2[5] = qqq2[6] = qqq2[7] = 0.0;
			}

			if (this->phys_param->bn_in_p_on_HP == true)
			{
				//cout << "A" << endl;
				std::vector<double> n(3);
				n[0] = gr->normal[now][0];
				n[1] = gr->normal[now][1];
				n[2] = gr->normal[now][2];

				this->phys_param->Godunov_Solver_Alexashov(qqq1, qqq2,//
					n, qqq, dsl, dsr, dsc, w);
				qqq[5] = qqq[6] = qqq[7] = 0.0;
			}
			else
			{
				//cout << "B" << endl;
				this->phys_param->chlld(metod_, gr->normal[now][0], gr->normal[now][1],
					gr->normal[now][2],
					w, qqq1, qqq2, qqq, false, 3,
					konvect_left, konvect_right, konvect, dsr, dsc, dsl,
					Option);
			}

			if (std::isnan(dsc) || std::fpclassify(dsc) == FP_SUBNORMAL)
			{
				cout << "Error  5436542867" << endl;
				whach(dsc);
				whach(dsl);
				whach(dsr);
				whach(qqq1[0]);
				whach(qqq2[0]);
				whach(qqq[0]);
				whach(qqq[1]);
				whach(qqq[2]);
				whach(qqq[3]);
				whach(qqq[4]);
				exit(-1);
			}

			for (auto& yz : gr->yzels)
			{
				yz->mut.lock();
				yz->velocity[0] += this->phys_param->velocity_HP * dsc * gr->normal[now][0];
				yz->velocity[1] += this->phys_param->velocity_HP * dsc * gr->normal[now][1];
				yz->velocity[2] += this->phys_param->velocity_HP * dsc * gr->normal[now][2];
				yz->num_velocity++;
				yz->mut.unlock();
			}

		}
	}

	// ������� �������� BS
	if (this->phys_param->move_BS == true)
	{
#pragma omp parallel for private(dsr, dsc, dsl)
		for (int i_step = 0; i_step < this->Gran_BS.size(); i_step++)
		{
			auto gr = this->Gran_BS[i_step];
			auto A = gr->cells[0];
			auto B = gr->cells[1];


			std::vector<double> qqq, qqq1, qqq2;
			qqq.resize(8);
			qqq1.resize(8);
			qqq2.resize(8);
			std::vector<double> konvect_left, konvect_right, konvect;
			PrintOptions Option = PrintOptions{};
			Option.fluid = "plasma_BS";

			// �������� ��������� ��� �������� ��������
			unordered_map<string, double> par_left, par_right;
			this->Snos_on_Gran(gr, par_left, par_right, now);

			if (false) // ��������� � ������� �����
			{
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
			else
			{
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
			}


			double w = 0.0;

			this->phys_param->chlld(gr->Get_method(), gr->normal[now][0], gr->normal[now][1],
				gr->normal[now][2],
				w, qqq1, qqq2, qqq, false, 3,
				konvect_left, konvect_right, konvect, dsr, dsc, dsl,
				Option);

			for (auto& yz : gr->yzels)
			{
				yz->mut.lock();
				yz->velocity[0] += 0.1 * dsr * gr->normal[now][0];
				yz->velocity[1] += 0.1 * dsr * gr->normal[now][1];
				yz->velocity[2] += 0.1 * dsr * gr->normal[now][2];
				yz->num_velocity++;
				yz->mut.unlock();
			}

		}
	}


	// ����������� ������������ (��� ��������� ����� ���������� ����� �� �����������)

	if (this->phys_param->sglag_TS == true)
	{
#pragma omp parallel for
		for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
		{
			Eigen::Vector3d A, B, V;
			auto gr = this->Gran_TS[i_step];
			A << gr->center[now][0], gr->center[now][1], gr->center[now][2];
			double rr = A.norm();
			double r = 0.0;
			for (auto& j : gr->grans_surf)
			{
				r += norm2(j->center[now][0], j->center[now][1], j->center[now][2]);
			}
			r /= gr->grans_surf.size();

			B = A * r / rr;
			
			double pk = this->phys_param->sglag_TS_k;
			double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
			if (phi > this->geo->tetta1)
			{
				pk = this->phys_param->sglag_TS_k_sphere_tail;
			}
			else if (phi < this->geo->tetta0)
			{
				pk = this->phys_param->sglag_TS_k_sphere_head;
			}

			V = this->phys_param->velocity_TS * pk * (B - A) / time;

			for (auto& yz : gr->yzels)
			{
				yz->mut.lock();
				yz->velocity[0] += V(0);
				yz->velocity[1] += V(1);
				yz->velocity[2] += V(2);
				yz->num_velocity++;
				yz->mut.unlock();
			}
		}

		
	}


	if (this->phys_param->sglag_HP == true)
	{
		double phi_a = polar_angle(this->A_Luch[0][1]->Yzels[0]->coord[0][0], 
			norm2(0.0, this->A_Luch[0][1]->Yzels[0]->coord[0][1], 
				this->A_Luch[0][1]->Yzels[0]->coord[0][2]));

		Eigen::Vector3d A, B, V;

		// ����������� � �������� ������� � > 0
		// ������ ������ ������ � ����������� ��. ��� ����������
		// ������ � ���������� �������� ����� ����� � "����������" �����������
		if (true) // ��� ������ ������� �����������, ������ �������� ������
		{
			
			if (true)
			{   // ������ � �����������
				#pragma omp parallel for private(A, B, V)
				for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
				{
					auto gr = this->Gran_HP[i_step];
					A << gr->center[now][0], gr->center[now][1], gr->center[now][2];

					double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));

					// �������� ���������-����������� ����������� ������ � �������� ����
					if (phi > phi_a + const_pi / 180.0) continue;


					double rr = A.norm();
					double r = 0.0;
					for (auto& j : gr->grans_surf)
					{
						r += norm2(j->center[now][0], j->center[now][1], j->center[now][2]);
					}
					r /= gr->grans_surf.size();

					B = A * r / rr;

					V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_k_sphere * (B - A) / time;

					for (auto& yz : gr->yzels)
					{
						yz->mut.lock();
						yz->velocity[0] += V(0);
						yz->velocity[1] += V(1);
						yz->velocity[2] += V(2);
						yz->num_velocity++;
						yz->mut.unlock();
					}
				}
			}
			if (true)
			{
				// ������ � ����������
				#pragma omp parallel for private(A, B, V)
				for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
				{
					auto gr = this->Gran_HP[i_step];
					A << gr->center[now][0], gr->center[now][1], gr->center[now][2];
					
					double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));

					// �������� ���������-����������� ����������� ������ � �������� ����
					if (phi > phi_a + const_pi / 180.0) continue;
					B = 2 * A;

					for (auto& j : gr->grans_surf)
					{
						B[0] += j->center[now][0];
						B[1] += j->center[now][1];
						B[2] += j->center[now][2];
					}
					B /= (gr->grans_surf.size() + 2.0);


					V = this->phys_param->velocity_HP * this->phys_param->sglag_HP_k_sphere * (B - A) / time;

					for (auto& yz : gr->yzels)
					{
						yz->mut.lock();
						yz->velocity[0] += V(0);
						yz->velocity[1] += V(1);
						yz->velocity[2] += V(2);
						yz->num_velocity++;
						yz->mut.unlock();
					}
				}
			}

		}

		Yzel* AA, * AA1, * AA11, * AA2, * AA22, * AA3, * AA33, * AA4, * AA44;
		// AA33  AA3  AA  AA1  AA11 - ����� HP  ����������� � ������
		// A22  AA2  AA  AA4  AA44  - �������� HP
		std::array<double, 2> p, p1, p11, p3, p33, p2, p22, p4, p44;
		double a, b, c;
		double r1, r2, r3, r4;

		// ����������� � �������� ������� (����� - ���������)
		if (false)
		{
#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
			//for (auto& yz : this->Yzels_HP_sglag)
			for (size_t idx = 0; idx < this->Yzels_HP_sglag.size(); ++idx)
			{
				auto& yz = this->Yzels_HP_sglag[idx];
				AA = yz;
				AA1 = yz->Yzel_sosed_sglag["AA1"];
				AA11 = yz->Yzel_sosed_sglag["AA11"];
				AA2 = yz->Yzel_sosed_sglag["AA3"];
				AA22 = yz->Yzel_sosed_sglag["AA33"];
				AA3 = yz->Yzel_sosed_sglag["AA2"];
				AA33 = yz->Yzel_sosed_sglag["AA22"];
				AA4 = yz->Yzel_sosed_sglag["AA4"];
				AA44 = yz->Yzel_sosed_sglag["AA44"];


				macros2(p, AA);
				macros2(p1, AA1);
				macros2(p11, AA11);
				macros2(p3, AA3);
				macros2(p33, AA33);

				solveQuadraticEquation(p3[0], p3[1],
					p1[0], p1[1], p11[0], p11[1], a, b, c);
				r1 = a * kv(p[0]) + b * p[0] + c;

				solveQuadraticEquation(p33[0], p33[1],
					p3[0], p3[1], p1[0], p1[1], a, b, c);
				r2 = a * kv(p[0]) + b * p[0] + c;


				p2[0] = polar_angle(AA2->coord[now][1], AA2->coord[now][2]);
				p2[1] = norm2(AA2->coord[now][0], AA2->coord[now][1], AA2->coord[now][2]);

				p[0] = polar_angle(AA->coord[now][1], AA->coord[now][2]);

				p4[0] = polar_angle(AA4->coord[now][1], AA4->coord[now][2]);
				p4[1] = norm2(AA4->coord[now][0], AA4->coord[now][1], AA4->coord[now][2]);

				p22[0] = polar_angle(AA22->coord[now][1], AA22->coord[now][2]);
				p22[1] = norm2(AA22->coord[now][0], AA22->coord[now][1], AA22->coord[now][2]);

				p44[0] = polar_angle(AA44->coord[now][1], AA44->coord[now][2]);
				p44[1] = norm2(AA44->coord[now][0], AA44->coord[now][1], AA44->coord[now][2]);

				solveQuadraticEquation(p2[0], p2[1],
					p4[0], p4[1], p44[0], p44[1], a, b, c);
				r3 = a * kv(p[0]) + b * p[0] + c;

				solveQuadraticEquation(p4[0], p4[1],
					p2[0], p2[1], p22[0], p22[1], a, b, c);
				r4 = a * kv(p[0]) + b * p[0] + c;

				A << AA->coord[now][0], AA->coord[now][1], AA->coord[now][2];

				B = A;
				B[0] *= r1 / p[1];
				B[1] *= r1 / p[1];
				B[2] *= r1 / p[1];

				V = this->phys_param->velocity_HP *
					this->phys_param->sglag_HP_sphere *
					this->phys_param->sglag_HP_k_sphere * (B - A) / time;

				AA->mut.lock();
				AA->velocity[0] += V[0];
				AA->velocity[1] += V[1];
				AA->velocity[2] += V[2];
				AA->num_velocity++;
				AA->mut.unlock();

				B = A;
				B[0] *= r2 / p[1];
				B[1] *= r2 / p[1];
				B[2] *= r2 / p[1];

				V = this->phys_param->velocity_HP *
					this->phys_param->sglag_HP_sphere *
					this->phys_param->sglag_HP_k_sphere * (B - A) / time; 

				AA->mut.lock();
				AA->velocity[0] += V[0];
				AA->velocity[1] += V[1];
				AA->velocity[2] += V[2];
				AA->num_velocity++;
				AA->mut.unlock();

				B = A;
				B[0] *= r3 / p[1];
				B[1] *= r3 / p[1];
				B[2] *= r3 / p[1];

				V = this->phys_param->velocity_HP *
					this->phys_param->sglag_HP_sphere *
					this->phys_param->sglag_HP_k_sphere * (B - A) / time; 

				AA->mut.lock();
				AA->velocity[0] += V[0];
				AA->velocity[1] += V[1];
				AA->velocity[2] += V[2];
				AA->num_velocity++;
				AA->mut.unlock();

				B = A;
				B[0] *= r4 / p[1];
				B[1] *= r4 / p[1];
				B[2] *= r4 / p[1];

				V = this->phys_param->velocity_HP *
					this->phys_param->sglag_HP_sphere *
					this->phys_param->sglag_HP_k_sphere * (B - A) / time; 

				AA->mut.lock();
				AA->velocity[0] += V[0];
				AA->velocity[1] += V[1];
				AA->velocity[2] += V[2];
				AA->num_velocity++;
				AA->mut.unlock();

			}
		}


		// �������������� ����������� �� A �����
		int nn = this->A_Luch.size();
		int mm = this->A_Luch[0].size();

#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
		for (size_t i = 0; i < nn; i++)
		{
			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 2; j < mm; j++)
			{
				AA = this->A_Luch[i][j]->Yzels_opor[2];
				AA2 = this->A_Luch[im][j]->Yzels_opor[2];
				AA22 = this->A_Luch[imm][j]->Yzels_opor[2];
				AA4 = this->A_Luch[ip][j]->Yzels_opor[2];
				AA44 = this->A_Luch[ipp][j]->Yzels_opor[2];

				if (j == mm - 1)
				{
					AA1 = this->A_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->B_Luch[i][0]->Yzels_opor[2];
					AA33 = this->B_Luch[i][1]->Yzels_opor[2];
				}
				else if (j == mm - 2)
				{
					AA1 = this->A_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->A_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->B_Luch[i][0]->Yzels_opor[2];
				}
				else
				{
					AA1 = this->A_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->A_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->A_Luch[i][j + 2]->Yzels_opor[2];
				}

				macros1();
			}
		}


		// �������������� ����������� �� B �����
		nn = this->B_Luch.size();
		mm = this->B_Luch[0].size();
#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
		for (size_t i = 0; i < nn; i++)
		{
			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 0; j < mm; j++)
			{
				AA = this->B_Luch[i][j]->Yzels_opor[2];
				AA2 = this->B_Luch[im][j]->Yzels_opor[2];
				AA22 = this->B_Luch[imm][j]->Yzels_opor[2];
				AA4 = this->B_Luch[ip][j]->Yzels_opor[2];
				AA44 = this->B_Luch[ipp][j]->Yzels_opor[2];

				if (j == 0)
				{
					int nnn = this->A_Luch[0].size();
					AA1 = this->A_Luch[i][nnn - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][nnn - 2]->Yzels_opor[2];
					AA3 = this->B_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->B_Luch[i][j + 2]->Yzels_opor[2];
				}
				else if (j == 1)
				{
					int nnn = this->A_Luch[0].size();
					AA1 = this->B_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->A_Luch[i][nnn - 1]->Yzels_opor[2];
					AA3 = this->B_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->B_Luch[i][j + 2]->Yzels_opor[2];
				}
				else if (j == mm - 1)
				{
					AA1 = this->B_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->B_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->E_Luch[i][0]->Yzels_opor[1];
					AA33 = this->E_Luch[i][1]->Yzels_opor[1];
				}
				else if (j == mm - 2)
				{
					AA1 = this->B_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->B_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->B_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->E_Luch[i][0]->Yzels_opor[1];
				}
				else
				{
					AA1 = this->B_Luch[i][j - 1]->Yzels_opor[2];
					AA11 = this->B_Luch[i][j - 2]->Yzels_opor[2];
					AA3 = this->B_Luch[i][j + 1]->Yzels_opor[2];
					AA33 = this->B_Luch[i][j + 2]->Yzels_opor[2];
				}

				macros1();
			}
		}
	
		// �������������� ����������� �� E �����
		nn = this->E_Luch.size();
		mm = this->E_Luch[0].size();

#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
		for (size_t i = 0; i < nn; i++)
		{
			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 0; j < mm; j++)
			{
				AA = this->E_Luch[i][j]->Yzels_opor[1];
				AA2 = this->E_Luch[im][j]->Yzels_opor[1];
				AA22 = this->E_Luch[imm][j]->Yzels_opor[1];
				AA4 = this->E_Luch[ip][j]->Yzels_opor[1];
				AA44 = this->E_Luch[ipp][j]->Yzels_opor[1];

				if (j == 0)
				{
					int nnn = this->B_Luch[0].size();
					AA1 = this->B_Luch[i][nnn - 1]->Yzels_opor[2];
					AA11 = this->B_Luch[i][nnn - 2]->Yzels_opor[2];
					AA3 = this->E_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->E_Luch[i][j + 2]->Yzels_opor[1];
				}
				else if (j == 1)
				{
					int nnn = this->B_Luch[0].size();
					AA1 = this->E_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->B_Luch[i][nnn - 1]->Yzels_opor[2];
					AA3 = this->E_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->E_Luch[i][j + 2]->Yzels_opor[1];
				}
				else if (j == mm - 1)
				{
					AA1 = this->E_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][j - 2]->Yzels_opor[1];
					AA3 = this->D_Luch[i][0]->Yzels_opor[1];
					AA33 = this->D_Luch[i][1]->Yzels_opor[1];
				}
				else if (j == mm - 2)
				{
					AA1 = this->E_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][j - 2]->Yzels_opor[1];
					AA3 = this->E_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->D_Luch[i][0]->Yzels_opor[1];
				}
				else
				{
					AA1 = this->E_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][j - 2]->Yzels_opor[1];
					AA3 = this->E_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->E_Luch[i][j + 2]->Yzels_opor[1];
				}
				macros1();
			}
		}

		// �������������� ����������� �� D �����
		nn = this->D_Luch.size();
		mm = this->geo->N4 - 1; // this->D_Luch[0].size();

#pragma omp parallel for private(A, B, V, a, b, c, r1, r2, r3, r4, p, p1, p11, p3, p33, p2, p22, p4, p44, AA, AA1, AA11, AA2, AA22, AA3, AA33, AA4, AA44)
		for (size_t i = 0; i < nn; i++)
		{
			short int ip = i + 1;
			short int ipp = i + 2;
			short int im = i - 1;
			short int imm = i - 2;

			if (ip == nn) ip = 0;
			if (ipp == nn) ipp = 0;
			if (ipp == nn + 1) ipp = 1;

			if (im == -1) im = nn - 1;
			if (imm == -1) imm = nn - 1;
			if (imm == -2) imm = nn - 2;

			for (size_t j = 0; j < mm; j++)
			{
				AA = this->D_Luch[i][j]->Yzels_opor[1];
				AA2 = this->D_Luch[im][j]->Yzels_opor[1];
				AA22 = this->D_Luch[imm][j]->Yzels_opor[1];
				AA4 = this->D_Luch[ip][j]->Yzels_opor[1];
				AA44 = this->D_Luch[ipp][j]->Yzels_opor[1];

				if (j == 0)
				{
					int nnn = this->E_Luch[0].size();
					AA1 = this->E_Luch[i][nnn - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][nnn - 2]->Yzels_opor[1];
					AA3 = this->D_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->D_Luch[i][j + 2]->Yzels_opor[1];
				}
				else if (j == 1)
				{
					int nnn = this->E_Luch[0].size();
					AA1 = this->D_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->E_Luch[i][nnn - 1]->Yzels_opor[1];
					AA3 = this->D_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->D_Luch[i][j + 2]->Yzels_opor[1];
				}
				else
				{
					AA1 = this->D_Luch[i][j - 1]->Yzels_opor[1];
					AA11 = this->D_Luch[i][j - 2]->Yzels_opor[1];
					AA3 = this->D_Luch[i][j + 1]->Yzels_opor[1];
					AA33 = this->D_Luch[i][j + 2]->Yzels_opor[1];
				}
				macros1();
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
		else if (yz->type == Type_yzel::BS && this->phys_param->move_BS == true)
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

	// ��������� ���� �� HP (������������ �����) ���� ���������
	if (this->phys_param->move_HP == true)
	{
		short int NN = this->D_Luch[0].size() - 1;
		for (auto& L : this->D_Luch)
		{

			double h1 = norm2(0.0, L[this->geo->N4 - 4]->Yzels_opor[1]->coord[now2][1],
				L[this->geo->N4 - 4]->Yzels_opor[1]->coord[now2][2]);

			// ������ �������� � ������ (������� ����������)
			//double h2 = norm2(0.0, L[NN]->Yzels_opor[1]->coord[now2][1],
			//	L[NN]->Yzels_opor[1]->coord[now2][2]);

			for (short int i = this->geo->N4 - 3; i <= NN; i++)
			{
				//double h = h1 + (i - this->geo->N4 + 2) * (h2 - h1) / (NN - this->geo->N4 + 2);
				double h = h1;
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
					//whach(h2);
				}


				yz->coord[now2][1] = yz->coord[now2][1] * h / hh;
				yz->coord[now2][2] = yz->coord[now2][2] * h / hh;
			}
		}
	}
	
	// ���� ��������� ����, ������� ���������� BS
	if (this->phys_param->move_BS == true)
	{
		short int NN = this->A_Luch[0].size() - 1;
		for (short int i = 0; i < this->A_Luch.size(); i++)
		{
			auto yz = this->A_Luch[i][NN]->Yzels_opor[3];
			double H = norm2(0.0, yz->coord[now2][1], yz->coord[now2][2]);
			for (auto& L : this->B_Luch[i])
			{
				auto yyz = L->Yzels_opor[3];
				double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
				yyz->coord[now2][1] *= H / h2;
				yyz->coord[now2][2] *= H / h2;
			}
			for (auto& L : this->E_Luch[i])
			{
				auto yyz = L->Yzels_opor[2];
				double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
				yyz->coord[now2][1] *= H / h2;
				yyz->coord[now2][2] *= H / h2;
			}
			for (auto& L : this->D_Luch[i])
			{
				auto yyz = L->Yzels_opor[2];
				double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
				yyz->coord[now2][1] *= H / h2;
				yyz->coord[now2][2] *= H / h2;
			}
		}


	}
}

void Setka::Smooth_head_TS(void)
{
	cout << "Start Smooth_head_TS" << endl;

	std::vector<double> a; // ������������ a1, a2, a3
	std::vector<double> b; // ������������ b1, b2, b3
	std::vector<double> c; // ������������ b1, b2, b3
	std::vector<double> d; // ������������ b1, b2, b3
	std::vector<double> e; // ������������ b1, b2, b3
	std::vector<double> f; // ������������ b1, b2, b3
	std::vector<double> w; // ����
	Eigen::Vector3d A, B, V;

	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			a.push_back(kv(A[0]));
			b.push_back(kv(A[1]));
			c.push_back(kv(A[2]));
			d.push_back(A[0] * A[1]);
			e.push_back(A[0] * A[2]);
			f.push_back(A[2] * A[1]);
			w.push_back(0.05 + phi);
		}
	}

	unsigned int n = a.size();
	Eigen::MatrixXd AA(n, 6);
	for (int i = 0; i < n; ++i) {
		AA(i, 0) = a[i] * w[i];  // ����������� ��� x
		AA(i, 1) = b[i] * w[i];  // ����������� ��� x
		AA(i, 2) = c[i] * w[i];  // ����������� ��� x
		AA(i, 3) = d[i] * w[i];  // ����������� ��� x
		AA(i, 4) = e[i] * w[i];  // ����������� ��� x
		AA(i, 5) = f[i] * w[i];  // ����������� ��� x
	}

	Eigen::VectorXd b_vec(n);
	for (int i = 0; i < n; ++i) {
		b_vec(i) = 1.0 * w[i];
	}

	// ������ ������� ������� ���������� ���������
	Eigen::VectorXd solution = AA.colPivHouseholderQr().solve(b_vec);

	// ��������� �������
	double x = solution(0);
	double y = solution(1);
	double z = solution(2);
	double zz = solution(3);
	double zzz = solution(4);
	double zzzz = solution(5);

	double ak = sqrt(1.0 / x);
	double bk = sqrt(1.0 / y);
	double ck = sqrt(1.0 / z);
	double dk = zz;
	double ek = zzz;
	double fk = zzzz;
	
	cout << "a = " << ak << endl;
	cout << "b = " << bk << endl;
	cout << "c = " << ck << endl;
	cout << "d = " << dk << endl;
	cout << "e = " << ek << endl;
	cout << "f = " << fk << endl;

	// ������ ���������� �����������

	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		for (auto& yz : gr->yzels)
		{
			A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			//double rr = A.norm();

			double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
			if (phi < this->geo->tetta0 + 0.17)
			{
				double rr = sqrt(kv(A[0]) / kv(ak) + kv(A[1]) / kv(bk)
					+ kv(A[2]) / kv(ck) + A[0] * A[1] * dk + 
					A[0] * A[2] * ek + A[2] * A[1] * fk);
				double aaa = 1.0/rr;
				yz->coord[0][0] *= aaa;
				yz->coord[0][1] *= aaa;
				yz->coord[0][2] *= aaa;
				yz->coord[1][0] = yz->coord[0][0];
				yz->coord[1][1] = yz->coord[0][1];
				yz->coord[1][2] = yz->coord[0][2];
			}
		}
	}


	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(0);
	}
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}


	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "End Smooth_head_TS" << endl;
}

void Setka::Smooth_head_TS2(void)
{
	cout << "Start Smooth_head_TS2" << endl;
	
	vector<unsigned int> nomera;
	vector<double> RRR;

	this->Renumerate();

	Eigen::Vector3d A, B, V;

	// ��� ������� r
	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			for (auto& yz : gr->yzels)
			{
				B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
				double r_yz = B.norm();
				double r_yz_new = 0.0;
				unsigned short int i_yz_new = 0;
				for (auto& ggr : yz->grans)
				{
					for (auto& yyz : ggr->yzels)
					{
						if (yyz->number == yz->number) continue;
						V << yyz->coord[0][0], yyz->coord[0][1], yyz->coord[0][2];
						double r_yyz = V.norm();
						r_yz_new += r_yyz;
						i_yz_new++;
						if (r_yz < r_yyz) goto b1;
					}
				}
				// ���� ����� �� ����, ������ ���� ���� ������������� ������������
				nomera.push_back(yz->number);
				RRR.push_back(r_yz_new/ i_yz_new);
			b1:
				r_yz_new = 0.0;
			}
		}
	}

	for (unsigned int j = 0; j < nomera.size(); j++)
	{
		unsigned int i = nomera[j];
		auto& yz = this->All_Yzel[i - 1];
		double rr = RRR[j];
		A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
		double r = A.norm();

		yz->coord[0][0] *= (rr/r);
		yz->coord[0][1] *= (rr / r);
		yz->coord[0][2] *= (rr / r);
		yz->coord[1][0] = yz->coord[0][0];
		yz->coord[1][1] = yz->coord[0][1];
		yz->coord[1][2] = yz->coord[0][2];
	}

	nomera.clear();
	RRR.clear();

	// ��� ����� R
	for (int i_step = 0; i_step < this->Gran_TS.size(); i_step++)
	{
		auto gr = this->Gran_TS[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			for (auto& yz : gr->yzels)
			{
				B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
				double r_yz = B.norm();
				double r_yz_new = 0.0;
				unsigned short int i_yz_new = 0;
				for (auto& ggr : yz->grans)
				{
					for (auto& yyz : ggr->yzels)
					{
						if (yyz->number == yz->number) continue;
						V << yyz->coord[0][0], yyz->coord[0][1], yyz->coord[0][2];
						double r_yyz = V.norm();
						r_yz_new += r_yyz;
						i_yz_new++;
						if (r_yz > r_yyz) goto bb1;
					}
				}
				// ���� ����� �� ����, ������ ���� ���� ������������� ������������
				nomera.push_back(yz->number);
				RRR.push_back(r_yz_new / i_yz_new);
			bb1:
				r_yz_new = 0.0;
			}
		}
	}

	for (unsigned int j = 0; j < nomera.size(); j++)
	{
		unsigned int i = nomera[j];
		auto& yz = this->All_Yzel[i - 1];
		double rr = RRR[j];
		A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
		double r = A.norm();

		yz->coord[0][0] *= (rr / r);
		yz->coord[0][1] *= (rr / r);
		yz->coord[0][2] *= (rr / r);
		yz->coord[1][0] = yz->coord[0][0];
		yz->coord[1][1] = yz->coord[0][1];
		yz->coord[1][2] = yz->coord[0][2];
	}



	// ---------------

	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(0);
	}
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}


	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "End Smooth_head_TS2" << endl;
}

void Setka::Smooth_head_HP2(void)
{
	cout << "Start Smooth_head_HP2" << endl;

	vector<unsigned int> nomera;
	vector<double> RRR;

	this->Renumerate();

	Eigen::Vector3d A, B, V;

	for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
	{
		auto gr = this->Gran_HP[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			for (auto& yz : gr->yzels)
			{
				B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
				double r_yz = B.norm();
				double r_yz_new = 0.0;
				unsigned short int i_yz_new = 0;
				for (auto& ggr : yz->grans)
				{
					for (auto& yyz : ggr->yzels)
					{
						if (yyz->number == yz->number) continue;
						V << yyz->coord[0][0], yyz->coord[0][1], yyz->coord[0][2];
						double r_yyz = V.norm();
						r_yz_new += r_yyz;
						i_yz_new++;
						if (r_yz < r_yyz) goto b1;
					}
				}
				// ���� ����� �� ����, ������ ���� ���� ������������� ������������
				nomera.push_back(yz->number);
				RRR.push_back(r_yz_new / i_yz_new);
			b1:
				r_yz_new = 0.0;
			}
		}
	}

	for (unsigned int j = 0; j < nomera.size(); j++)
	{
		unsigned int i = nomera[j];
		auto& yz = this->All_Yzel[i - 1];
		double rr = RRR[j];
		A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
		double r = A.norm();

		yz->coord[0][0] *= (rr / r);
		yz->coord[0][1] *= (rr / r);
		yz->coord[0][2] *= (rr / r);
		yz->coord[1][0] = yz->coord[0][0];
		yz->coord[1][1] = yz->coord[0][1];
		yz->coord[1][2] = yz->coord[0][2];
	}

	// -----------------------------------
	nomera.clear();
	RRR.clear();

	// ��� ����� R
	for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
	{
		auto gr = this->Gran_HP[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			for (auto& yz : gr->yzels)
			{
				B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
				double r_yz = B.norm();
				double r_yz_new = 0.0;
				unsigned short int i_yz_new = 0;
				for (auto& ggr : yz->grans)
				{
					for (auto& yyz : ggr->yzels)
					{
						if (yyz->number == yz->number) continue;
						V << yyz->coord[0][0], yyz->coord[0][1], yyz->coord[0][2];
						double r_yyz = V.norm();
						r_yz_new += r_yyz;
						i_yz_new++;
						if (r_yz > r_yyz) goto bb1;
					}
				}
				// ���� ����� �� ����, ������ ���� ���� ������������� ������������
				nomera.push_back(yz->number);
				RRR.push_back(r_yz_new / i_yz_new);
			bb1:
				r_yz_new = 0.0;
			}
		}
	}

	for (unsigned int j = 0; j < nomera.size(); j++)
	{
		unsigned int i = nomera[j];
		auto& yz = this->All_Yzel[i - 1];
		double rr = RRR[j];
		A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
		double r = A.norm();

		yz->coord[0][0] *= (rr / r);
		yz->coord[0][1] *= (rr / r);
		yz->coord[0][2] *= (rr / r);
		yz->coord[1][0] = yz->coord[0][0];
		yz->coord[1][1] = yz->coord[0][1];
		yz->coord[1][2] = yz->coord[0][2];
	}
	// ------------------------------------
	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(0);
	}
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}


	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "End Smooth_head_HP2" << endl;
}

void Setka::Smooth_head_HP(void)
{
	cout << "Start Smooth_head_HP" << endl;

	std::vector<double> a; // ������������ a1, a2, a3
	std::vector<double> b; // ������������ b1, b2, b3
	std::vector<double> c; // ������������ b1, b2, b3
	std::vector<double> d; // ������������ b1, b2, b3
	std::vector<double> e; // ������������ b1, b2, b3
	std::vector<double> f; // ������������ b1, b2, b3
	std::vector<double> w; // ����
	Eigen::Vector3d A, B, V;

	for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
	{
		auto gr = this->Gran_HP[i_step];
		A << gr->center[0][0], gr->center[0][1], gr->center[0][2];
		double rr = A.norm();

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
		if (phi < this->geo->tetta0 + 0.17)
		{
			a.push_back(kv(A[0]));
			b.push_back(kv(A[1]));
			c.push_back(kv(A[2]));
			d.push_back(A[0] * A[1]);
			e.push_back(A[0] * A[2]);
			f.push_back(A[2] * A[1]);
			w.push_back(0.05 + phi);
		}
	}

	unsigned int n = a.size();
	Eigen::MatrixXd AA(n, 6);
	for (int i = 0; i < n; ++i) {
		AA(i, 0) = a[i] * w[i];  // ����������� ��� x
		AA(i, 1) = b[i] * w[i];  // ����������� ��� x
		AA(i, 2) = c[i] * w[i];  // ����������� ��� x
		AA(i, 3) = d[i] * w[i];  // ����������� ��� x
		AA(i, 4) = e[i] * w[i];  // ����������� ��� x
		AA(i, 5) = f[i] * w[i];  // ����������� ��� x
	}

	Eigen::VectorXd b_vec(n);
	for (int i = 0; i < n; ++i) {
		b_vec(i) = 1.0 * w[i];
	}

	// ������ ������� ������� ���������� ���������
	Eigen::VectorXd solution = AA.colPivHouseholderQr().solve(b_vec);

	// ��������� �������
	double x = solution(0);
	double y = solution(1);
	double z = solution(2);
	double zz = solution(3);
	double zzz = solution(4);
	double zzzz = solution(5);

	double ak = sqrt(1.0 / x);
	double bk = sqrt(1.0 / y);
	double ck = sqrt(1.0 / z);
	double dk = zz;
	double ek = zzz;
	double fk = zzzz;

	cout << "a = " << ak << endl;
	cout << "b = " << bk << endl;
	cout << "c = " << ck << endl;
	cout << "d = " << dk << endl;
	cout << "e = " << ek << endl;
	cout << "f = " << fk << endl;

	// ������ ���������� �����������

	for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
	{
		auto gr = this->Gran_HP[i_step];
		for (auto& yz : gr->yzels)
		{
			A << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			//double rr = A.norm();

			double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));
			if (phi < this->geo->tetta0 + 0.17)
			{
				double rr = sqrt(kv(A[0]) / kv(ak) + kv(A[1]) / kv(bk)
					+ kv(A[2]) / kv(ck) + A[0] * A[1] * dk +
					A[0] * A[2] * ek + A[2] * A[1] * fk);
				double aaa = 1.0 / rr;
				yz->coord[0][0] *= aaa;
				yz->coord[0][1] *= aaa;
				yz->coord[0][2] *= aaa;
				yz->coord[1][0] = yz->coord[0][0];
				yz->coord[1][1] = yz->coord[0][1];
				yz->coord[1][2] = yz->coord[0][2];
			}
		}
	}


	for (int i_step = 0; i_step < this->All_Luch.size(); i_step++)
	{
		auto lu = this->All_Luch[i_step];
		lu->dvigenie(0);
	}
	for (auto& i : this->All_Yzel)
	{
		for (unsigned short int j = 0; j < 3; j++)
		{
			i->coord[1][j] = i->coord[0][j];
		}
	}


	this->Calculating_measure(0);
	this->Calculating_measure(1);

	cout << "End Smooth_head_HP" << endl;
}

void Setka::Find_Yzel_Sosed_for_sglag(void)
{
	cout << "Start: Find_Yzel_Sosed_for_sglag" << endl;
	this->Renumerate();
	short int now = 0;
	Eigen::Vector3d A;
	Eigen::Vector3d B, V1, V2, VV;
	double phi_a = polar_angle(this->A_Luch[0][1]->Yzels[0]->coord[0][0],
		norm2(0.0, this->A_Luch[0][1]->Yzels[0]->coord[0][1],
			this->A_Luch[0][1]->Yzels[0]->coord[0][2]));

	// ������� ������� ��� HP
	for (int i_step = 0; i_step < this->Gran_HP.size(); i_step++)
	{
		auto gr = this->Gran_HP[i_step];
		A << gr->center[now][0], gr->center[now][1], gr->center[now][2];

		double phi = polar_angle(A[0], norm2(0.0, A[1], A[2]));

		// �������� ���������-����������� ����������� ������ � �������� ����
		if (phi > phi_a + const_pi / 180.0) continue;

		//  ��� ��������� ����� ����� ������� �����������
		for (auto& yz : gr->yzels)
		{
			if (yz->Yzel_sosed_sglag.size() > 0)
			{
				continue;
			}
			Yzel* AA3, * AA1;
			Yzel* AA2, * AA4;
			AA1 = nullptr;
			AA2 = nullptr;
			AA3 = nullptr;
			AA4 = nullptr;
			// AA33  AA3  (AA)  AA1  AA11 - ����� HP  ����������� � ������
			// A22  AA2  (AA)  AA4  AA44  - �������� HP
			double d1 = -1.0;
			double d2 = 1.0;
			double d3 = -1.0;
			double d4 = 1.0;

			B << yz->coord[0][0], yz->coord[0][1], yz->coord[0][2];
			dekard_skorost(B[0], B[1], B[2],
				0.0, 0.0, 1.0,
				V1[0], V1[1], V1[2]);
			dekard_skorost(B[0], B[1], B[2],
				0.0, 1.0, 0.0,
				V2[0], V2[1], V2[2]);
			V1.normalize();
			V2.normalize();
			for (auto& grr : gr->grans_surf)
			{
				int num = yz->number;
				if (std::find_if(
					grr->yzels.begin(),
					grr->yzels.end(),
					[num](Yzel* yzz) { return yzz->number == num; })
					== grr->yzels.end())
				{
					continue;
				}

				for (auto& yzz : grr->yzels)
				{
					if (yzz->number == yz->number) continue;
					VV << yzz->coord[0][0] - B[0], yzz->coord[0][1] - B[1],
						yzz->coord[0][2] - B[2];
					VV.normalize();
					double sk1 = VV.dot(V1);
					if (sk1 > d1)
					{
						d1 = sk1;
						AA1 = yzz;
					}
					if (sk1 < d2)
					{
						d2 = sk1;
						AA2 = yzz;
					}

					double sk2 = VV.dot(V2);
					if (sk2 > d3)
					{
						d3 = sk2;
						AA3 = yzz;
					}
					if (sk2 < d4)
					{
						d4 = sk2;
						AA4 = yzz;
					}
				}
			}

			if (AA1 != nullptr && AA2 != nullptr && AA3 != nullptr && AA4 != nullptr)
			{
				yz->Yzel_sosed_sglag["AA1"] = AA1;
				yz->Yzel_sosed_sglag["AA2"] = AA2;
				yz->Yzel_sosed_sglag["AA3"] = AA3;
				yz->Yzel_sosed_sglag["AA4"] = AA4;
				this->Yzels_HP_sglag.push_back(yz);
			}
		}
	}
	
	for (auto& yz : this->Yzels_HP_sglag)
	{
		auto s1 = yz->Yzel_sosed_sglag["AA1"];
		auto s2 = yz->Yzel_sosed_sglag["AA2"];
		auto s3 = yz->Yzel_sosed_sglag["AA3"];
		auto s4 = yz->Yzel_sosed_sglag["AA4"];


		if (s1->Yzel_sosed_sglag.find("AA1") != s1->Yzel_sosed_sglag.end() &&
			s2->Yzel_sosed_sglag.find("AA2") != s2->Yzel_sosed_sglag.end() &&
			s3->Yzel_sosed_sglag.find("AA3") != s3->Yzel_sosed_sglag.end() &&
			s4->Yzel_sosed_sglag.find("AA4") != s4->Yzel_sosed_sglag.end())
		{
			yz->Yzel_sosed_sglag["AA11"] = s1->Yzel_sosed_sglag["AA1"];
			yz->Yzel_sosed_sglag["AA22"] = s2->Yzel_sosed_sglag["AA2"];
			yz->Yzel_sosed_sglag["AA33"] = s3->Yzel_sosed_sglag["AA3"];
			yz->Yzel_sosed_sglag["AA44"] = s4->Yzel_sosed_sglag["AA4"];
		}
		else
		{
			yz->Yzel_sosed_sglag.clear();
		}

	}
	
	// ������� ����, � ������� �� ������� ������ �������
	this->Yzels_HP_sglag.erase(
		std::remove_if(
			this->Yzels_HP_sglag.begin(),
			this->Yzels_HP_sglag.end(),
			[](Yzel* yz) { return yz->Yzel_sosed_sglag.size() == 0; }
		),
		this->Yzels_HP_sglag.end()
	);


	cout << "End: Find_Yzel_Sosed_for_sglag" << endl;
}
