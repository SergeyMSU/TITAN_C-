#include "Setka.h"


void Setka::Set_MK_Zone(void)
{
	this->Renumerate();
	Eigen::Vector3d Centr;

	this->Cell_Center->MK_zone_r = 1;
	this->Cell_Center->MK_zone_phi = 0;

	// Задаём зону для каждой ячейки
	for (auto& cell : this->All_Cell)
	{
		Centr[0] = cell->center[0][0];
		Centr[1] = cell->center[0][1];
		Centr[2] = cell->center[0][2];

		if (cell->type == Type_cell::Zone_1)
		{
			cell->MK_zone_r = 1;
			cell->MK_zone_phi = 0;
		}
		else if (cell->type == Type_cell::Zone_2)
		{
			cell->MK_zone_r = 2;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
			}
			else
			{
				cell->MK_zone_phi = 2;
			}
		}
		else if (cell->type == Type_cell::Zone_3)
		{
			cell->MK_zone_r = 3;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
			}
			else
			{
				cell->MK_zone_phi = 2;
			}
		}
		else if (cell->type == Type_cell::Zone_4)
		{
			cell->MK_zone_r = 4;
			if (Centr[0] > 0)
			{
				cell->MK_zone_phi = 1;
			}
			else
			{
				cell->MK_zone_phi = 2;
			}
		}
		else
		{
			cout << "Error 86743207564" << endl;
			exit(-1);
		}
	}

	this->MK_Grans.resize(7);

	// 1 зона
	for (auto& gr : this->Gran_TS)
	{
		this->MK_Grans[0].push_back(gr);
	}

	// 2 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (Centr[0] < -0.000001) continue;

		if (gr->type2 == Type_Gran_surf::HP) this->MK_Grans[1].push_back(gr);


		if (gr->type2 == Type_Gran_surf::TS) this->MK_Grans[1].push_back(gr);
		

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2) this->MK_Grans[1].push_back(gr);
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_2)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[1].push_back(gr);
				}
			}
		}
	}

	// 3 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (Centr[0] > 0.000001) continue;

		if (gr->type2 == Type_Gran_surf::HP) this->MK_Grans[2].push_back(gr);

		if (gr->type2 == Type_Gran_surf::TS) this->MK_Grans[2].push_back(gr);

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2) this->MK_Grans[2].push_back(gr);
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_3) 
			{
				this->MK_Grans[2].push_back(gr);
			}

			
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_2)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[2].push_back(gr);
				}
			}
		}
	}

	// 4 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::HP)
		{
			if(Centr[0] > 0) this->MK_Grans[3].push_back(gr);
		}

		if (gr->type2 == Type_Gran_surf::BS) this->MK_Grans[3].push_back(gr);

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[3].push_back(gr);
				}
			}

		}
	}

	// 5 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::HP)
		{
			if (Centr[0] < 0.0) this->MK_Grans[4].push_back(gr);
		}

		if (gr->type == Type_Gran::Outer_Soft)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3) this->MK_Grans[4].push_back(gr);
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_2 &&
				gr->cells[1]->type == Type_cell::Zone_3) this->MK_Grans[4].push_back(gr);


			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_3)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[4].push_back(gr);
				}
			}

			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if(gr->cells[0]->center[0][0] < 0) this->MK_Grans[4].push_back(gr);
			}

		}
	}

	// 6 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type2 == Type_Gran_surf::BS)
		{
			this->MK_Grans[5].push_back(gr);
		}

		if (gr->type == Type_Gran::Outer_Hard)
		{
			if(Centr[0] > 0) this->MK_Grans[5].push_back(gr);
		}

		if (gr->cells.size() == 2)
		{
			if (gr->cells[0]->type == Type_cell::Zone_4 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[5].push_back(gr);
				}
			}

		}
	}

	// 7 зона
	for (auto& gr : this->All_Gran)
	{
		Centr[0] = gr->center[0][0];
		Centr[1] = gr->center[0][1];
		Centr[2] = gr->center[0][2];

		if (gr->type == Type_Gran::Outer_Hard )
		{
			if (Centr[0] < 0) this->MK_Grans[6].push_back(gr);
		}

		if (gr->type != Type_Gran::Us && Centr[0] < 0)
		{
			if(gr->cells[0]->type == Type_cell::Zone_4) this->MK_Grans[6].push_back(gr);
		}

		if (gr->cells.size() == 2 && gr->type2 == Type_Gran_surf::Us)
		{
			if (gr->cells[0]->type == Type_cell::Zone_3 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] < 0) this->MK_Grans[6].push_back(gr);
			}

			if (gr->cells[0]->type == Type_cell::Zone_4 &&
				gr->cells[1]->type == Type_cell::Zone_4)
			{
				if (gr->cells[0]->center[0][0] * gr->cells[1]->center[0][0] < 0.0)
				{
					this->MK_Grans[6].push_back(gr);
				}
			}

		}


	}
}