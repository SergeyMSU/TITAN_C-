#include "Setka.h"

void Setka::Init_boundary_grans(void)
{
	this->Calculating_measure(0);

	for (auto& i : this->All_Gran)
	{
		i->type = Type_Gran::Us;

		if (i->cells.size() != 1) continue;

		if (i->func_R(0) <= this->geo->R0 + 0.00001)
		{
			i->type = Type_Gran::Inner_Hard;
			this->All_boundary_Gran.push_back(i);
		}
		else if(i->center[0][0] > this->geo->L7 + 0.0001)
		{
			i->type = Type_Gran::Outer_Hard;
			this->All_boundary_Gran.push_back(i);
		}
		else
		{
			i->type = Type_Gran::Outer_Soft;
			this->All_boundary_Gran.push_back(i);
		}
	}
}

void Setka::Init_physics(void)
{
	this->Calculating_measure(0);
	double x, y, z, r, the;
	Eigen::Vector3d vec, cc, vv;
	double BR, BPHI, V1, V2, V3, mV;

	// ����� ��������� ������� �� �����
	if (true)
	{
		for (auto& i : this->All_Cell)
		{
			x = i->center[0][0];
			y = i->center[0][1];
			z = i->center[0][2];
			r = i->func_R(0);

			if (false)//(r < this->geo->R0 * 30.0)
			{

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

	// ����� ��������� ������� (�� ��������� ������)
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

				the = -the + const_pi / 2.0;   // �.�.� ������ �� �� �� 1 �.�.���� �� - 90 �� 90 � ����������
				cc = this->phys_param->Matr * vv;

				mV = this->phys_param->Get_v_0(the);

				i->parameters["rho"] = this->phys_param->Get_rho_0(the);
				i->parameters["p"] = this->phys_param->p_0;
				i->parameters["Vx"] = mV * vec(0);
				i->parameters["Vy"] = mV * vec(1);
				i->parameters["Vz"] = mV * vec(2);
				i->parameters["Bx"] = cc(0);
				i->parameters["By"] = cc(1);
				i->parameters["Bz"] = cc(2);
				i->parameters["Q"] = i->parameters["rho"];
			}
		}
	}
}

void Setka::Go(void)
{
	this->Test_geometr();
	// ��� ��������� ������� ����������� �� ����� Setter.txt
	// ������� ������������� ������ ��� �������� �����
	unsigned int steps = 30;

	unsigned short int now1 = 1;
	unsigned short int now2 = 0;

	Cell* A, B;
	std::vector<double> qqq, qqq1, qqq2;
	std::vector<double> POTOK;
	qqq.resize(8);
	qqq1.resize(8);
	qqq2.resize(8);
	POTOK.resize(9);
	std::vector<double> konvect_left, konvect_right, konvect;
	PrintOptions Option = PrintOptions{};
	double dsr, dsc, dsl;


	double time = 0.000001;  // ������� ��� �� �������
	double loc_time = 0.000001;  // ������� ��� �� �������


	for (unsigned int step = 1; step <= steps; step++)
	{
		cout << "Global step = " << step << endl;
		whach(time);

		now2 = now1;
		now1 = (now1 + 1)%2;

		time = loc_time;
		loc_time = 100000000.0;


		// ����������� ������ ����� �����
		for (auto& gran : this->All_Gran)
		{
			double area = gran->area[now1];
			if (gran->type == Type_Gran::Us) // ������� �����
			{
				auto A = gran->cells[0];
				auto B = gran->cells[1];

				if (A->parameters[now1].find("rho") == A->parameters[now1].end() ||
					B->parameters[now1].find("rho") == B->parameters[now1].end())
				{
					cout << "Error  0956453978" << endl;
					exit(-1);
				}

				if (std::isnan(A->parameters[now1]["rho"]) || std::fpclassify(A->parameters[now1]["rho"]) == FP_SUBNORMAL)
				{
					cout << "Error 0907675453" << endl;
					for (const auto& [key, value] : A->parameters[now1]) {
						std::cout << key << ": " << value << std::endl;
					}
					exit(-1);
				}

				qqq1[0] = A->parameters[now1]["rho"];
				qqq1[1] = A->parameters[now1]["Vx"];
				qqq1[2] = A->parameters[now1]["Vy"];
				qqq1[3] = A->parameters[now1]["Vz"];
				qqq1[4] = A->parameters[now1]["p"];
				qqq1[5] = A->parameters[now1]["Bx"];
				qqq1[6] = A->parameters[now1]["By"];
				qqq1[7] = A->parameters[now1]["Bz"];

				qqq2[0] = B->parameters[now1]["rho"];
				qqq2[1] = B->parameters[now1]["Vx"];
				qqq2[2] = B->parameters[now1]["Vy"];
				qqq2[3] = B->parameters[now1]["Vz"];
				qqq2[4] = B->parameters[now1]["p"];
				qqq2[5] = B->parameters[now1]["Bx"];
				qqq2[6] = B->parameters[now1]["By"];
				qqq2[7] = B->parameters[now1]["Bz"];

				double w = 0.0;

				this->phys_param->chlld(1, gran->normal[now1][0], gran->normal[now1][1], gran->normal[now1][2],
					w, qqq1, qqq2, qqq, false, 1,
					konvect_left, konvect_right, konvect, dsr, dsc, dsl,
					Option);

				double dist = norm2(A->center[now1][0] - B->center[now1][0],
					A->center[now1][1] - B->center[now1][1],
					A->center[now1][2] - B->center[now1][2])/2.0;

				loc_time = min(loc_time, 0.9 * dist / (max(fabs(dsl), fabs(dsr)) + fabs(w)));

				gran->parameters["Pm"] = qqq[0] * area;
				gran->parameters["PVx"] = qqq[1] * area;
				gran->parameters["PVy"] = qqq[2] * area;
				gran->parameters["PVz"] = qqq[3] * area;
				gran->parameters["Pe"] = qqq[4] * area;
				gran->parameters["PBx"] = qqq[5] * area;
				gran->parameters["PBy"] = qqq[6] * area;
				gran->parameters["PBz"] = qqq[7] * area;
				gran->parameters["PdivB"] = 0.5 * scalarProductFast(gran->normal[now1][0], 
					gran->normal[now1][1], gran->normal[now1][2], 
					qqq1[5] + qqq2[5], qqq1[6] + qqq2[6], qqq1[7] + qqq2[7]) * area;
			}
			else if (gran->type == Type_Gran::Inner_Hard || gran->type == Type_Gran::Outer_Hard)
			{
				this->phys_param->Get_Potok(gran->parameters["rho"], gran->parameters["p"],
					gran->parameters["Vx"], gran->parameters["Vy"], gran->parameters["Vz"],
					gran->parameters["Bx"], gran->parameters["By"], gran->parameters["Bz"],
					gran->normal[now1][0], gran->normal[now1][1], gran->normal[now1][2],
					qqq);

				gran->parameters["Pm"] = qqq[0] * area;
				gran->parameters["PVx"] = qqq[1] * area;
				gran->parameters["PVy"] = qqq[2] * area;
				gran->parameters["PVz"] = qqq[3] * area;
				gran->parameters["Pe"] = qqq[7] * area;
				gran->parameters["PBx"] = qqq[4] * area;
				gran->parameters["PBy"] = qqq[5] * area;
				gran->parameters["PBz"] = qqq[6] * area;

				gran->parameters["PdivB"] = scalarProductFast(gran->normal[now1][0],
					gran->normal[now1][1], gran->normal[now1][2],
					gran->parameters["Bx"], gran->parameters["By"], 
					gran->parameters["Bz"]) * area;
			}
			else if (gran->type == Type_Gran::Outer_Soft)
			{
				auto C = gran->cells[0];
				gran->parameters["rho"] = C->parameters[now1]["rho"];
				gran->parameters["p"] = C->parameters[now1]["p"];
				gran->parameters["Vx"] = C->parameters[now1]["Vx"];
				gran->parameters["Vy"] = C->parameters[now1]["Vy"];
				gran->parameters["Vz"] = C->parameters[now1]["Vz"];
				gran->parameters["Bx"] = C->parameters[now1]["Bx"];
				gran->parameters["By"] = C->parameters[now1]["By"];
				gran->parameters["Bz"] = C->parameters[now1]["Bz"];
				this->phys_param->Get_Potok(gran->parameters["rho"], gran->parameters["p"],
					gran->parameters["Vx"], gran->parameters["Vy"], gran->parameters["Vz"],
					gran->parameters["Bx"], gran->parameters["By"], gran->parameters["Bz"],
					gran->normal[now1][0], gran->normal[now1][1], gran->normal[now1][2],
					qqq);

				gran->parameters["Pm"] = qqq[0] * area;
				gran->parameters["PVx"] = qqq[1] * area;
				gran->parameters["PVy"] = qqq[2] * area;
				gran->parameters["PVz"] = qqq[3] * area;
				gran->parameters["Pe"] = qqq[7] * area;
				gran->parameters["PBx"] = qqq[4] * area;
				gran->parameters["PBy"] = qqq[5] * area;
				gran->parameters["PBz"] = qqq[6] * area;
				gran->parameters["PdivB"] = scalarProductFast(gran->normal[now1][0],
					gran->normal[now1][1], gran->normal[now1][2],
					gran->parameters["Bx"], gran->parameters["By"],
					gran->parameters["Bz"]) * area;
			}
			else
			{
				cout << "Error 6323145345" << endl;
				exit(-1);
			}
		}

		// ����������� ������ ���������� � �������
		for (auto& cell : this->All_Cell)
		{
			double Volume = cell->volume[now1];
			double Volume2 = cell->volume[now2];

			POTOK[0] = POTOK[1] = POTOK[2] = POTOK[3] = POTOK[4] = POTOK[5] = 
				POTOK[6] = POTOK[7] = POTOK[8] = 0.0;

			for (auto& gran : cell->grans)
			{
				short int sign_potok = 1;
				if(gran->cells[0] != cell) sign_potok = -1;

				POTOK[0] += sign_potok * gran->parameters["Pm"];
				POTOK[1] += sign_potok * gran->parameters["PVx"];
				POTOK[2] += sign_potok * gran->parameters["PVy"];
				POTOK[3] += sign_potok * gran->parameters["PVz"];
				POTOK[4] += sign_potok * gran->parameters["Pe"];
				POTOK[5] += sign_potok * gran->parameters["PBx"];
				POTOK[6] += sign_potok * gran->parameters["PBy"];
				POTOK[7] += sign_potok * gran->parameters["PBz"];
				POTOK[8] += sign_potok * gran->parameters["PdivB"];
			}

			

			double rho3, u3, v3, w3, bx3, by3, bz3, p3;
			double rho = cell->parameters[now1]["rho"];
			double vx = cell->parameters[now1]["Vx"];
			double vy = cell->parameters[now1]["Vy"];
			double vz = cell->parameters[now1]["Vz"];
			double p = cell->parameters[now1]["p"];
			double bx = cell->parameters[now1]["Bx"];
			double by = cell->parameters[now1]["By"];
			double bz = cell->parameters[now1]["Bz"];
			double dsk = scalarProductFast(vx, vy, vz, bx, by, bz);

			rho3 = rho * Volume / Volume2
				- time * POTOK[0] / Volume2;

			u3 = (rho * vx * Volume / Volume2 - time * (POTOK[1] + (bx / cpi4) * POTOK[8]) / Volume2) / rho3;
			v3 = (rho * vy * Volume / Volume2 - time * (POTOK[2] + (by / cpi4) * POTOK[8]) / Volume2) / rho3;
			w3 = (rho * vz * Volume / Volume2 - time * (POTOK[3] + (bz / cpi4) * POTOK[8]) / Volume2) / rho3;

			bx3 = bx * Volume / Volume2 - time * (POTOK[5] + vx * POTOK[8]) / Volume2;
			by3 = by * Volume / Volume2 - time * (POTOK[6] + vy * POTOK[8]) / Volume2;
			bz3 = bz * Volume / Volume2 - time * (POTOK[7] + vz * POTOK[8]) / Volume2;

			p3 = (((p / this->phys_param->g1 + 0.5 * rho * kvv(vx, vy, vz) + kvv(bx, by, bz) / 25.13274122871834590768) * Volume / Volume2
				- time * (POTOK[4] + (dsk / cpi4) * POTOK[8]) / Volume2) -
				0.5 * rho3 * kvv(u3, v3, w3) - kvv(bx3, by3, bz3) / 25.13274122871834590768) * this->phys_param->g1;

			if (std::isnan(rho3) || std::fpclassify(rho3) == FP_SUBNORMAL)
			{
				whach(Volume);
				whach(Volume2);
				whach(rho);
				whach(time);

				for (short unsigned int ik = 0; ik < 9; ik++)
				{
					cout << "ik: " << POTOK[ik] << endl;
				}
				exit(-1);
			}


			if (p3 < 0.0000000001)
			{
				p3 = 0.000001;
			}

			cell->parameters[now2]["rho"] = rho3;
			cell->parameters[now2]["Vx"] = u3;
			cell->parameters[now2]["Vy"] = v3;
			cell->parameters[now2]["Vz"] = w3;
			cell->parameters[now2]["Bx"] = bx3;
			cell->parameters[now2]["By"] = by3;
			cell->parameters[now2]["Bz"] = bz3;
			cell->parameters[now2]["p"] = p3;
		
		}
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
		// ���������� ���������� ���������
		size_t size = i->parameters[0].size();
		out.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

		// ���������� ������ ���� ����-��������
		for (const auto& pair : i->parameters[0]) {
			// ������� ���������� ����� ����� � ��� ����
			size_t key_size = pair.first.size();
			out.write(reinterpret_cast<const char*>(&key_size), sizeof(size_t));
			out.write(pair.first.c_str(), key_size);

			// ����� ���������� ��������
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
		// ������ ���������� ���������
		size_t size;
		in.read(reinterpret_cast<char*>(&size), sizeof(size_t));

		for (size_t i = 0; i < size; ++i) {
			// ������ ����
			size_t key_size;
			in.read(reinterpret_cast<char*>(&key_size), sizeof(size_t));

			std::vector<char> key_buffer(key_size);
			in.read(key_buffer.data(), key_size);
			std::string key(key_buffer.begin(), key_buffer.end());

			// ������ ��������
			double value;
			in.read(reinterpret_cast<char*>(&value), sizeof(double));

			ii->parameters[0][key] = value;
		}
	}
}
