#include "Setka.h"


void Setka::Write_file_for_FCMHD(void)
{
	this->Renumerate();

	cout << "Start: Write_file_for_FCMHD" << endl;

	ofstream file("FCMHD_A3.bin", ios::binary);
	if (!file.is_open())
	{
		cout << "ERROR FCMHD_2.bin" << endl;
		exit(-1);
	}

	int n1 = this->All_Cell.size();
	int n2 = this->All_Gran.size();

	file.write(reinterpret_cast<const char*>(&this->phys_param->T_all), sizeof(double));
	file.write(reinterpret_cast<const char*>(&n1), sizeof(int));
	file.write(reinterpret_cast<const char*>(&n2), sizeof(int));

	size_t host_N_cell = this->All_Cell.size();
	size_t host_N_gran = this->All_Gran.size();
	const vector<string> param_order = { "rho", "Vx", "Vy", "Vz", "p", "Bx", "By", "Bz" };

	// host_Cell_par
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (const string& param : param_order)
		{
			double value = this->All_Cell[i]->parameters[0].at(param);
			if (i == 0)
			{
				cout << "S1 - " << value << endl;
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "host_Cell_par success" << endl;

	// host_Cell_center
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (short int k = 0; k < 3; k++)
		{
			double value = this->All_Cell[i]->center[0][k];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "host_Cell_center success" << endl;

	// host_Cell_Volume
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (short int k = 0; k < 1; k++)
		{
			double value = this->All_Cell[i]->volume[k];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "host_Cell_Volume success" << endl;

	// host_Cell_gran
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (short int k = 0; k < 6; k++)
		{
			int value = this->All_Cell[i]->grans[k]->number;
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	cout << "host_Cell_gran success" << endl;

	// host_Gran_normal
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 3; k++)
		{
			double value = this->All_Gran[i]->normal[0][k];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "host_Gran_normal success" << endl;

	// host_Gran_square
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 1; k++)
		{
			double value = this->All_Gran[i]->area[0];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "host_Gran_square success" << endl;


	// host_Gran_center
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 3; k++)
		{
			double value = this->All_Gran[i]->center[0][k];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "host_Gran_center success" << endl;

	// host_Gran_neighbour
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 2; k++)
		{
			Cell* cc = nullptr;
			if (this->All_Gran[i]->cells.size() >= k + 1) cc = this->All_Gran[i]->cells[k];
			int value = 0;
			if (cc != nullptr)
			{
				value = cc->number;
			}

			if (value > host_N_cell)
			{
				cout << "ERROR 38u4rh8fur   " << value << " " << host_N_cell << endl;
				cout << cc->number << endl;
				cout << this->All_Gran[i]->cells.size() << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	cout << "host_Gran_neighbour success" << endl;

	// host_Gran_neighbour_TVD
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 2; k++)
		{
			Cell* cc = nullptr;
			if (this->All_Gran[i]->cells_TVD.size() >= k + 1) cc = this->All_Gran[i]->cells_TVD[k];
			int value = 0;
			if (cc != nullptr)
			{
				value = cc->number;
			}

			if (i == 1)
			{
				cout << "S2 - " << value << endl;
			}

			if (value < 0 || value > host_N_cell)
			{
				cout << "ERRORRcv 8u39874fy8eh" << endl;
				cout << value << endl;
				exit(-1);
			}

			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	cout << "host_Gran_neighbour_TVD success" << endl;

	double vv = 121.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

	// host_Gran_type
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		for (short int k = 0; k < 1; k++)
		{
			auto cc = this->All_Gran[i];
			int value = 0;
			if (cc->type == Type_Gran::Us) value = 1;
			if (cc->type == Type_Gran::Inner_Hard) value = 2;
			if (cc->type == Type_Gran::Outer_Soft) value = 3;
			if (cc->type == Type_Gran::Outer_Hard) value = 4;

			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	cout << "host_Gran_type success" << endl;

	vv = 122.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

	const vector<string> param_order2 = { "Prho", "PVx", "PVy", "PVz", "Pp", "PBx", "PBy", "PBz", "PdivB" };

	// host_Gran_POTOK
	for (size_t i = 0; i < host_N_gran; ++i)
	{
		auto cc = this->All_Gran[i];
		string name;
		double time = Culc_Gran_Potok(cc, 0, 3, name, 1.0);
		double value = 0.0;

		for (const string& param : param_order2)
		{
			value = cc->parameters[param];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "host_Gran_POTOK success" << endl;

	// Проверяющий параметр
	vv = 123.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

	file.close();
}

void Setka::Read_file_for_FCMHD(void)
{
	this->Renumerate();

	std::ifstream file("FCMHD_A3.1_out.bin", std::ios::binary);
	if (!file.is_open()) 
	{
		file.open("CUDA_FORT/FCMHD_A3.1_out.bin", std::ios::binary);
	}

	if (!file.is_open()) 
	{
		throw std::runtime_error("Error opening file: FCMHD______out.bin");
	}

	int n1 = this->All_Cell.size();
	int n2 = this->All_Gran.size();

	file.read(reinterpret_cast<char*>(&this->phys_param->T_all), sizeof(double));

	size_t host_N_cell = this->All_Cell.size();
	size_t host_N_gran = this->All_Gran.size();
	const vector<string> param_order = { "rho", "Vx", "Vy", "Vz", "p", "Bx", "By", "Bz" };

	// host_Cell_par
	for (size_t i = 0; i < host_N_cell; ++i)
	{
		for (const string& param : param_order)
		{
			double value;
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			this->All_Cell[i]->parameters[0][param] = value;
		}
	}

	// Проверяющий параметр
	double vv;
	file.read(reinterpret_cast<char*>(&vv), sizeof(double));

	cout << "Proverka (321) = " << vv << endl;

	file.close();
}

