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

	// Задаём начальные условия на сетке
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
	// Все настройки расчёта считываются из файла Setter.txt
}