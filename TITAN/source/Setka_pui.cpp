#include "Setka.h"


void Setka::Init_h0_and_read_from_file(void)
{
	this->h0_pui.resize(this->phys_param->pui_h0_n);
	std::ifstream infile("h0_for_pui.bin", std::ios::binary | std::ios::in);

	if (!infile.is_open()) {
		std::cerr << "Error: Cannot open file " << "h0_for_pui.bin" << " for reading" << std::endl;
		exit(-1);
	}

	try {
		// Считываем размер вектора
		size_t size = 0;
		infile.read(reinterpret_cast<char*>(&size), sizeof(size_t));

		if (infile.gcount() != sizeof(size_t)) {
			std::cerr << "Error: Failed to read vector size" << std::endl;
			exit(-1);
		}

		// Изменяем размер вектора и считываем данные
		h0_pui.resize(size);
		infile.read(reinterpret_cast<char*>(h0_pui.data()), size * sizeof(double));

		if (infile.gcount() != static_cast<std::streamsize>(size * sizeof(double))) {
			std::cerr << "Error: Failed to read vector data" << std::endl;
			exit(-1);
		}

		// Проверяем контрольное число
		int vb_check = 0;
		infile.read(reinterpret_cast<char*>(&vb_check), sizeof(int));

		if (infile.gcount() != sizeof(int)) {
			std::cerr << "Error: Failed to read verification number" << std::endl;
			exit(-1);
		}

		if (vb_check != 123) {
			std::cerr << "Error: Verification failed. File may be corrupted." << std::endl;
			exit(-1);
		}
	}
	catch (const std::exception& e) 
	{
		std::cerr << "Error during file reading: " << e.what() << std::endl;
		exit(-1);
	}
}

void Setka::Culc_h0_for_pui(void)
{
	vector<double> h0_pui(this->phys_param->pui_h0_n);

#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < this->phys_param->pui_h0_n; k++)
	{
		double S = 0.0;
		double UH = (k + 0.5) * this->phys_param->pui_wR / this->phys_param->pui_h0_n;

		for (int i = 0; i < 1000; i++) {
			double w = (i + 0.5) * this->phys_param->pui_wR / 1000.0;

			for (int j = 0; j < 180; j++)
			{
				double the = j * const_pi / 180.0;
				double u = sqrt(pow(w * sin(the), 2) + pow(w * cos(the) - UH, 2));
				S = max(S, u * this->phys_param->sigma(u) /
					((w + this->phys_param->pui_h0_wc) * this->phys_param->sigma(w + this->phys_param->pui_h0_wc)));
			}
		}
		h0_pui[k] = S;
	}


	// Запись результатов в файл
	std::ofstream outfile("h0_for_pui.bin", std::ios::binary | std::ios::out);

	if (!outfile.is_open()) {
		std::cerr << "Error  tybherbvhenbevrrtb564656 " << std::endl;
		return;
	}

	// Записываем размер вектора (количество элементов)
	size_t size = h0_pui.size();
	outfile.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

	// Записываем данные вектора
	outfile.write(reinterpret_cast<const char*>(h0_pui.data()), size * sizeof(double));

	// Для проверку успешности считывания
	int vb = 123;
	outfile.write(reinterpret_cast<const char*>(&vb), sizeof(vb));

	outfile.close();

	std::ofstream outfile2("h0_PUI.txt");
	for (int k = 0; k < this->phys_param->pui_h0_n; k++)
	{
		double UH = (k - 0.5) * this->phys_param->pui_wR / this->phys_param->pui_h0_n;
		outfile2 << UH << " " << h0_pui[k] << std::endl;
	}
	outfile2.close();
}

double Setka::PUI_get_h0(const double& w)
{
	const double Wmax = this->phys_param->pui_wR;
	short int N = 0;
	N = this->h0_pui.size(); 
	if (N == 0) return 0.0;

	double cell_size = Wmax / N;
	int left_index = static_cast<int>(w / cell_size);
	if (left_index == N)
	{
		return this->h0_pui[N - 1];
	}

	// Координаты центров ячеек
	double left_center = (left_index + 0.5) * cell_size;
	double right_center = (left_index + 1.5) * cell_size;

	if (left_index == 0 && w < left_center) {
		// Экстраполяция от первой ячейки
		double next_center = (1.5) * cell_size;
		return this->h0_pui[0] + (this->h0_pui[1] - this->h0_pui[0]) *
			(w - left_center) / (next_center - left_center);
	}
	else if (left_index == N - 1 && w > right_center)
	{
		double prev_center = (N - 1.5) * cell_size;
		return this->h0_pui[N - 2] + (this->h0_pui[N - 1] - this->h0_pui[N - 2]) *
			(w - prev_center) / (right_center - prev_center);
	}

	// Линейная интерполяция между центрами ячеек
	double t = (w - left_center) / (right_center - left_center);
	return this->h0_pui[left_index] * (1.0 - t) + this->h0_pui[left_index + 1] * t;
}

void Setka::Delete_h0(void)
{
	this->h0_pui.clear();
	this->h0_pui.resize(0);
}

bool Setka::Get_pui_SS(vector<double>& pui_Sm, vector<double>& pui_Sp1, vector<double>& pui_Sp2, 
	short int ii, const double& x, const double& y, const double& z,
	Setka& S_MK, Interpol& SI_MK, Cell_handle& prev_cell, Cell_handle& next_cell)
{
	vector<int> num_cell(4);
	vector<double> koeff_cell(4);

	//cout << "Get_pui_Sm  1" << endl;
	bool b = SI_MK.Get_real_cells(x, y, z, num_cell, koeff_cell, prev_cell, next_cell);
	//cout << "Get_pui_Sm  2" << endl;

	if (b == false) return false;

	prev_cell = next_cell;

	Cell* CC;

	for (int jj = 0; jj < this->phys_param->pui_nW; jj++)
	{
		double S = 0.0;
		double SS = 0.0;
		double SSS = 0.0;

		for (int i = 0; i < num_cell.size(); i++)
		{
			int j = num_cell[i];
			double k = koeff_cell[i];

			CC = S_MK.All_Cell[j];
			S += CC->pui_Sm[jj] * k;
			SS += CC->pui_Sp(0, jj) * k;
			if (ii == 2)
			{
				int rows = CC->pui_Sp.rows();
				if (rows >= 2)
				{
					SSS += CC->pui_Sp(1, jj) * k;
				}
				else
				{
					SSS += 0.0;
				}
			}
		}

		pui_Sm[jj] = S;
		pui_Sp1[jj] = SS;
		if (ii == 2) pui_Sp2[jj] = SSS;
	}

	return true;
}

bool Setka::Get_pui_Sm(double& pui_Sm, int n, const double& x, const double& y, const double& z,
	Setka& S_MK, Interpol& SI_MK, Cell_handle& prev_cell, Cell_handle& next_cell)
{
	vector<int> num_cell(4);
	vector<double> koeff_cell(4);

	//cout << "Get_pui_Sm  1" << endl;
	bool b = SI_MK.Get_real_cells(x, y, z, num_cell, koeff_cell, prev_cell, next_cell);
	//cout << "Get_pui_Sm  2" << endl;

	if (b == false) return false;

	prev_cell = next_cell;

	double S = 0.0;
	//cout << "Get_pui_Sm  3" << endl;
	for (int i = 0; i < num_cell.size(); i++)
	{
		int j = num_cell[i];
		double k = koeff_cell[i];

		//cout << "Get_pui_Sm  4" << endl;

		//cout << "j = " << j << "   " << S_MK.All_Cell.size() << endl;
		//cout << S_MK.All_Cell[j]->pui_Sm.size() << "   " << n << endl;
		S += S_MK.All_Cell[j]->pui_Sm[n] * k;


		//cout << "Get_pui_Sm  5" << endl;
	}
	//cout << "Get_pui_Sm  6" << endl;

	pui_Sm = S;
	return true;
}

bool Setka::Get_pui_Sp(double& pui_Sp, short int ii, int n, const double& x, const double& y, const double& z,
	Setka& S_MK, Interpol& SI_MK, Cell_handle& prev_cell, Cell_handle& next_cell)
{
	vector<int> num_cell;
	vector<double> koeff_cell;

	bool b = SI_MK.Get_real_cells(x, y, z, num_cell, koeff_cell, prev_cell, next_cell);

	if (b == false) return false;

	prev_cell = next_cell;

	double S = 0.0;
	for (int i = 0; i < num_cell.size(); i++)
	{
		int j = num_cell[i];
		double k = koeff_cell[i];

		S += S_MK.All_Cell[j]->pui_Sp(ii, n) * k;
	}

	pui_Sp = S;
	return true;
}