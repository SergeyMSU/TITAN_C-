#include "Cell.h"


void Cell::Init_mas_pogl(short int n, short int nH)
{
	this->mas_pogl.resize(nH, n);
	this->mas_pogl.setZero();
}

void Cell::Delete_mas_pogl()
{
	this->mas_pogl.resize(0, 0);
}

void Cell::write_mas_pogl_ToFile(Phys_param* phys_param)
{
	std::string filename = phys_param->pogl_folder + "/poglosh_" + to_string(this->number) + ".bin";


	std::ofstream file(filename, std::ios::binary);
	if (!file.is_open()) 
	{
		std::cerr << "Error ergeryr5y45tertgwergwtw  " << filename << std::endl;
		std::cerr << "Error code: " << strerror(errno) << std::endl; // Добавьте эту строку
		return;
		exit(-1);
	}

	try 
	{

		// Записываем размерности матрицы (rows и cols)
		int rows = this->mas_pogl.rows();
		int cols = this->mas_pogl.cols();
		file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
		file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));

		// Записываем данные mas_pogl
		file.write(reinterpret_cast<const char*>(this->mas_pogl.data()),
			rows * cols * sizeof(double));

		int pvr = 1345;
		file.write(reinterpret_cast<const char*>(&pvr), sizeof(pvr));


		file.close();
	}
	catch (const std::exception& e) {
		std::cerr << "Error  4y546uy46y65ye5y4" << e.what() << std::endl;
		file.close();
		exit(-1);
	}
}

void Cell::read_mas_pogl_FromFile(Phys_param* phys_param)
{
	std::string filename = phys_param->pogl_folder + "/poglosh_" + to_string(this->number) + ".bin";

	if (file_exists(filename) == false)
	{
		//cout << "Not file: " << filename << endl;
		// cout << "Not file: " << this->center[0][0] << " " << this->center[0][1] << " " << this->center[0][2] << endl;
		return;
	}

	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open())
	{
		std::cerr << "Error tgrthrtgrthfth" << filename << std::endl;
		exit(-10);
	}

	try {

		// Читаем размерности матрицы
		int rows, cols;
		file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
		file.read(reinterpret_cast<char*>(&cols), sizeof(cols));

		if (rows != this->mas_pogl.rows() || cols != this->mas_pogl.cols()) {
			// Вариант 1: Выбросить исключение
			cout << "Error  ergetrfgrewfwregwerg " << endl;
			exit(-1);
		}

		// Читаем данные
		file.read(reinterpret_cast<char*>(this->mas_pogl.data()),
			rows * cols * sizeof(double));


		int pvr;
		file.read(reinterpret_cast<char*>(&pvr), sizeof(pvr));
		if (pvr != 1345)
		{
			cout << "Error  efdwef34rfwefwerfewwfwf" << endl;
			exit(-1);
		}

		//cout << "0 Sum = " << this->mas_pogl.sum() << endl;

		// Проверка есть ли ненулевой элемент:
		if (false)
		{
			bool has_non_zero = this->mas_pogl.isZero();

			if (has_non_zero)
			{
				cout << "All null  " << this->center[0][0] << " " << this->center[0][1] << " " << this->center[0][2] <<  endl;
			}
		}

		file.close();
	}
	catch (const std::exception& e) {
		std::cerr << "Error rthrtyhgrt45ty45t " << e.what() << std::endl;
		file.close();
		exit(-19);
	}
}

int Cell::pogl_mas_number(const double& Ve, const double& L, const double& R, short int n)
{
	// Проверка на выход за границы
	if (Ve < L) return 0; // левее левой границы
	if (Ve >= R) return n - 1; // правее правой границы (или последняя ячейка)

	// Вычисление номера ячейки
	double cellWidth = (R - L) / n;
	int index = static_cast<int>((Ve - L) / cellWidth);

	// Обеспечиваем, что индекс в диапазоне [0, n-1]
	return std::min(std::max(index, 0), n - 1);
}


void Cell::Init_pui_integral(short int n, short int zone)
{
	this->F_integr_pui_1.resize(n);
	this->nu_integr_pui_1.resize(n);
	this->Mz_integr_pui_1.resize(n);
	this->E_integr_pui_1.resize(n);

	for (size_t i = 0; i < n; i++)
	{
		this->F_integr_pui_1[i] = 0.0;
		this->nu_integr_pui_1[i] = 0.0;
		this->Mz_integr_pui_1[i] = 0.0;
		this->E_integr_pui_1[i] = 0.0;
	}

	if (zone == 2)
	{
		this->F_integr_pui_2.resize(n);
		this->nu_integr_pui_2.resize(n);
		this->Mz_integr_pui_2.resize(n);
		this->E_integr_pui_2.resize(n);

		for (size_t i = 0; i < n; i++)
		{
			this->F_integr_pui_2[i] = 0.0;
			this->nu_integr_pui_2[i] = 0.0;
			this->Mz_integr_pui_2[i] = 0.0;
			this->E_integr_pui_2[i] = 0.0;
		}
	}
}

void Cell::Delete_pui_integral(void)
{
	this->F_integr_pui_1.resize(0);
	this->nu_integr_pui_1.resize(0);
	this->Mz_integr_pui_1.resize(0);
	this->E_integr_pui_1.resize(0);

	this->F_integr_pui_2.resize(0);
	this->nu_integr_pui_2.resize(0);
	this->Mz_integr_pui_2.resize(0);
	this->E_integr_pui_2.resize(0);
}

void Cell::Init_f_pui(short int n, short int zone)
{
	this->f_pui_1.resize(n);

	for (size_t i = 0; i < n; i++)
	{
		this->f_pui_1[i] = 0.0;
	}

	if (zone == 2)
	{
		this->f_pui_2.resize(n);

		for (size_t i = 0; i < n; i++)
		{
			this->f_pui_2[i] = 0.0;
		}
	}
}

void Cell::Delete_f_pui(void)
{
	this->f_pui_1.resize(0);
	this->f_pui_2.resize(0);
}

void Cell::Init_S(short int k, short int n)
{
	this->pui_Sm.resize(n);
	this->pui_Sp.resize(k, n);

	for (size_t i = 0; i < n; i++)
	{
		this->pui_Sm[i] = 0.0;
	}
	this->pui_Sp.setZero();
}

void Cell::write_S_ToFile(void)
{
	std::string filename = "data_SpSm/func_cells_SpSm_" + to_string(this->number) + ".bin";


	std::ofstream file(filename, std::ios::binary);
	if (!file.is_open()) {
		std::cerr << "Error jdheirgfvbierboifeurfvef " << filename << std::endl;
		exit(-1);
	}

	try {
		// Записываем размеры матриц
		int rows = this->pui_Sp.rows();
		int cols = this->pui_Sp.cols();
		file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
		file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));

		int size = this->pui_Sm.size();
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));


		// Записываем данные pui_Sm
		file.write(reinterpret_cast<const char*>(this->pui_Sm.data()),
			size * sizeof(double));

		// Записываем данные pui_Sp
		file.write(reinterpret_cast<const char*>(this->pui_Sp.data()),
			rows * cols * sizeof(double));

		file.close();
	}
	catch (const std::exception& e) {
		std::cerr << "Error  KJHfuiyregbuyifbeihrfyugwuhuerf" << e.what() << std::endl;
		file.close();
		exit(-1);
	}
}

void Cell::write_pui_integral_ToFile(void)
{
	std::string filename = "data_pui_intergal/func_cells_pui_integral_" + to_string(this->number) + ".bin";


	std::ofstream file(filename, std::ios::binary);
	if (!file.is_open()) 
	{
		std::cerr << "Error egertgerg2345745yerg " << filename << std::endl;
		exit(-1);
	}

	try {
		// Записываем размеры матриц
		int size1 = this->F_integr_pui_1.size();
		int size2 = this->F_integr_pui_2.size();
		file.write(reinterpret_cast<const char*>(&size1), sizeof(size1));
		file.write(reinterpret_cast<const char*>(&size2), sizeof(size2));

		file.write(reinterpret_cast<const char*>(this->F_integr_pui_1.data()),
			size1 * sizeof(double));
		file.write(reinterpret_cast<const char*>(this->nu_integr_pui_1.data()),
			size1 * sizeof(double));
		file.write(reinterpret_cast<const char*>(this->Mz_integr_pui_1.data()),
			size1 * sizeof(double));
		file.write(reinterpret_cast<const char*>(this->E_integr_pui_1.data()),
			size1 * sizeof(double));

		file.write(reinterpret_cast<const char*>(this->F_integr_pui_2.data()),
			size2 * sizeof(double));
		file.write(reinterpret_cast<const char*>(this->nu_integr_pui_2.data()),
			size2 * sizeof(double));
		file.write(reinterpret_cast<const char*>(this->Mz_integr_pui_2.data()),
			size2 * sizeof(double));
		file.write(reinterpret_cast<const char*>(this->E_integr_pui_2.data()),
			size2 * sizeof(double));
		

		int prv = 468;
		file.write(reinterpret_cast<const char*>(&prv), sizeof(prv));

		file.close();
	}
	catch (const std::exception& e) {
		std::cerr << "Error  t4y6545tgdfgwer34t" << e.what() << std::endl;
		file.close();
		exit(-1);
	}
}

void Cell::read_pui_integral_FromFile(Phys_param*& Phys_param)
{
	std::string filename = "data_pui_intergal/func_cells_pui_integral_" + to_string(this->number) + ".bin";

	if (file_exists(filename) == false) return;

	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open())
	{
		std::cout << "Error ji86y4tre4364ye5tg" << filename << std::endl;
		exit(-10);
	}

	try {
		// Читаем размеры матриц
		int size1, size2;
		file.read(reinterpret_cast<char*>(&size1), sizeof(size1));
		file.read(reinterpret_cast<char*>(&size2), sizeof(size2));

		// Проверяем корректность прочитанных размеров
		if (size1 < 0 || size2 < 0) {
			std::cout << "Error: Invalid sizes read from file: " << size1 << ", " << size2 << std::endl;
			file.close();
			exit(-13);
		}

		if (size1 != this->F_integr_pui_1.size() || size2 != this->F_integr_pui_2.size())
		{
			std::cout << "Error reerhyy6576456t45uhrdgrt" << std::endl;
			cout << size1 << " " << size2 << " " << this->F_integr_pui_1.size() << " " << 
				this->F_integr_pui_2.size() << endl;
			file.close();
			exit(-10);
		}


		file.read(reinterpret_cast<char*>(this->F_integr_pui_1.data()),
			size1 * sizeof(double));
		file.read(reinterpret_cast<char*>(this->nu_integr_pui_1.data()),
			size1 * sizeof(double));
		file.read(reinterpret_cast<char*>(this->Mz_integr_pui_1.data()),
			size1 * sizeof(double));
		file.read(reinterpret_cast<char*>(this->E_integr_pui_1.data()),
			size1 * sizeof(double));

		file.read(reinterpret_cast<char*>(this->F_integr_pui_2.data()),
			size2 * sizeof(double));
		file.read(reinterpret_cast<char*>(this->nu_integr_pui_2.data()),
			size2 * sizeof(double));
		file.read(reinterpret_cast<char*>(this->Mz_integr_pui_2.data()),
			size2 * sizeof(double));
		file.read(reinterpret_cast<char*>(this->E_integr_pui_2.data()),
			size2 * sizeof(double));

		int pvr;
		file.read(reinterpret_cast<char*>(&pvr), sizeof(pvr));
		if (pvr != 468)
		{
			cout << "Error  93498tryh38hgibuwe4fghp9iehg" << endl;
			exit(-1);
		}

		file.close();
	}
	catch (const std::exception& e) 
	{
		std::cout << "Error ry4536345te5y34ty34et;9 " << e.what() << std::endl;
		file.close();
		exit(-19);
	}

	if (this->number == 378182)
	{
		this->print_nu_integr_pui(Phys_param, "_init_");
	}

}

void Cell::pui_integral_Culc(Phys_param* phys_param)
{
	// Вычисляем сначала норму  \rho_w (см. документацию "PUI")
	for (int ijk = 0; ijk < 2; ijk++)
	{
		if (ijk == 1 && F_integr_pui_2.size() == 0) continue;

		double SS = 0.0;
		unsigned short int nn1 = 1000;
		double ff, w, S, UH;

		for (int i = 0; i < nn1; i++)
		{
			w = (i + 0.5) * phys_param->pui_wR / nn1;
			ff = pui_get_f(w, ijk, phys_param->pui_wR);
			SS += ff * kv(w) * (w + phys_param->pui_h0_wc) * phys_param->sigma(w + phys_param->pui_h0_wc)
				* (phys_param->pui_wR / nn1);
		}

		double S1 = 0.0;
		double w_val = 0.0;
		for (int i = 0; i < phys_param->pui_F_n; i++)
		{
			//cout << "i = " << i << endl;
			double S2 = (i + 0.5) * 1.0 / phys_param->pui_F_n;

			while (S1 < S2)
			{
				w_val += (phys_param->pui_wR / nn1);
				//cout << "get " << w_val << " " << phys_param->pui_wR << endl;
				ff = pui_get_f(w_val, ijk, phys_param->pui_wR);
				//cout << "ff = " << ff << endl;
				S1 += ff * kv(w_val) * (w_val + phys_param->pui_h0_wc) *
					phys_param->sigma(w_val + phys_param->pui_h0_wc) * (phys_param->pui_wR / nn1) / SS;
			}

			if (ijk == 0) this->F_integr_pui_1[i] = w_val;
			if (ijk == 1) this->F_integr_pui_2[i] = w_val;
		}


		// Далее надо считать частоту перезарядки и источники импульса и энергии
		double the, u;
		nn1 = 500;
		short int nn2 = 90;
		for (int i = 0; i < phys_param->pui_F_n; i++)
		{
			//cout << "i = " << i << endl;
			S = 0.0;
			SS = 0.0;
			S1 = 0.0;
			UH = (i + 0.5) * phys_param->pui_wR / phys_param->pui_F_n;
			for (int j = 0; j < nn1; j++)
			{
				//cout << "j = " << i << endl;
				w = (j + 0.5) * phys_param->pui_wR / nn1;
				ff = pui_get_f(w, ijk, phys_param->pui_wR);
				for (int k = 0; k < nn2; k++)
				{
					//cout << "j = " << i << endl;
					the = const_pi * k / nn2;
					u = sqrt(kv(w * sin(the)) + kv(w * cos(the) - UH));
					u = max(u, 0.000001);
					double ks = const_pi * ff * kv(w) * u * phys_param->sigma(u) * sin(the) * (const_pi / nn2) * (phys_param->pui_wR / nn1);
					S = S + 2.0 * ks;
					SS = SS + 2.0 * w * cos(the) * ks;
					S1 = S1 + kv(w) * ks;

					if (this->number == 2201 && i == 0)
					{
						//cout << "u = " << u << " " << ks << " " << ff << " " << w << " " << the << endl;
					}
				}
			}

			if (ijk == 0)
			{
				if (this->number == 2201)
				{
					//cout << "S = " << S << endl;
				}
				this->nu_integr_pui_1[i] = S;
				this->Mz_integr_pui_1[i] = SS;
				this->E_integr_pui_1[i] = S1;
			}
			else if (ijk == 1)
			{
				this->nu_integr_pui_2[i] = S;
				this->Mz_integr_pui_2[i] = SS;
				this->E_integr_pui_2[i] = S1;
			}
			else
			{
				cout << "Error iehrg9uhet7834yt" << endl;
				exit(-1);
			}
		}
	}

}

void Cell::print_F_integr_pui(string name)
{
	std::ofstream file(name + "_1_print_F_integr_pui_" + to_string(this->number) + ".txt");
	if (!file.is_open()) {
		cout << "Error uehgrifbghvuoyerfowhefwef" << endl;
		exit(-1);
	}

	std::ofstream file3(name + "_1_print_F_integr_interpol_pui_" + to_string(this->number) + ".txt");
	if (!file3.is_open()) {
		cout << "Error uehgrifbghvuoyerfowhefwef" << endl;
		exit(-1);
	}

	double pui_wR = 1.0;

	for (int i = 0; i < F_integr_pui_1.size(); i++)
	{
		double w = (i + 0.5) * pui_wR / F_integr_pui_1.size();
		double w2 = this->PUI_get_F_integer(w, 0);
		file << w << "\t" << F_integr_pui_1[i] << "\t" << w2 << "\n";
	}

	for (double w = 0.0; w <= 1.0; w = w + 0.001)
	{
		double w2 = this->PUI_get_F_integer(w, 0);
		file3 << w << "\t" << w2 << "\n";
	}

	file.close();
	file3.close();

	std::ofstream  file2(name + "_2_print_F_integr_pui_" + to_string(this->number) + ".txt");
	if (!file2.is_open()) {
		cout << "Error dthgretget345" << endl;
		exit(-1);
	}

	std::ofstream  file4(name + "_2_print_F_integr_interpol_pui_" + to_string(this->number) + ".txt");
	if (!file4.is_open()) {
		cout << "Error dthgretget345" << endl;
		exit(-1);
	}

	for (int i = 0; i < F_integr_pui_2.size(); i++)
	{
		double w = (i + 0.5) * pui_wR / F_integr_pui_2.size();
		double w2 = this->PUI_get_F_integer(w, 1);
		file2 << w << "\t" << F_integr_pui_2[i] << "\t" << w2 << "\n";
	}

	for (double w = 0.0; w <= 1.0; w = w + 0.001)
	{
		double w2 = this->PUI_get_F_integer(w, 1);
		file4 << w << "\t" << w2 << "\n";
	}

	file2.close();
	file4.close();
}


void Cell::print_nu_integr_pui(Phys_param* phys_param, string name)
{
	std::ofstream file(name + "_1_print_nu_integr_pui_" + to_string(this->number) + ".txt");
	if (!file.is_open()) {
		cout << "Error uehgrifbghvuoyerfowhefwef" << endl;
		exit(-1);
	}

	for (int i = 0; i < nu_integr_pui_1.size(); i++) 
	{
		double w = (i + 0.5) * phys_param->pui_wR / nu_integr_pui_1.size();
		file << w << "\t" << nu_integr_pui_1[i] << "\n";
	}

	file.close();

	std::ofstream  file2(name + "_2_print_nu_integr_pui_" + to_string(this->number) + ".txt");
	if (!file2.is_open()) {
		cout << "Error dthgretget345" << endl;
		exit(-1);
	}

	for (int i = 0; i < nu_integr_pui_2.size(); i++)
	{
		double w = (i + 0.5) * phys_param->pui_wR / nu_integr_pui_2.size();
		file2 << w << "\t" << nu_integr_pui_2[i] << "\n";
	}

	file2.close();
}

void Cell::write_pui_ToFile(void)
{
	std::string filename = "data_pui/func_cells_pui_" + to_string(this->number) + ".bin";


	std::ofstream file(filename, std::ios::binary);
	if (!file.is_open()) {
		std::cerr << "Error thrthrtghery4523234saedfw423 " << filename << std::endl;
		exit(-1);
	}

	try {
		// Записываем размеры матриц
		int size1 = this->f_pui_1.size();
		int size2 = this->f_pui_2.size();
		file.write(reinterpret_cast<const char*>(&size1), sizeof(size1));
		file.write(reinterpret_cast<const char*>(&size2), sizeof(size2));

		file.write(reinterpret_cast<const char*>(this->f_pui_1.data()),
			size1 * sizeof(double));

		file.write(reinterpret_cast<const char*>(this->f_pui_2.data()),
			size2 * sizeof(double));

		file.close();
	}
	catch (const std::exception& e) {
		std::cerr << "Error  rthr6yu45756874t5ruhjrt435" << e.what() << std::endl;
		file.close();
		exit(-1);
	}
}

void Cell::read_S_FromFile(const double& n_H_lism) 
{
	std::string filename = "data_SpSm/func_cells_SpSm_" + to_string(this->number) + ".bin";

	if (file_exists(filename) == false) return;

	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open()) 
	{
		std::cerr << "Error 98u98y8fwgrihfbuoerfre" << filename << std::endl;
		exit(-10);
	}

	try {
		// Читаем размеры матриц
		int rows, cols, size;
		file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
		file.read(reinterpret_cast<char*>(&cols), sizeof(cols));

		file.read(reinterpret_cast<char*>(&size), sizeof(size));

		// Проверяем соответствие размеров
		if (rows != this->pui_Sp.rows() || cols != this->pui_Sp.cols())
		{
			std::cerr << "Error ejighieurgerg54t4t5"
				<< "fail: " << rows << "  x   " << cols << ", "
				<< "expect: " << this->pui_Sp.rows() << "  x  " << this->pui_Sp.cols() << std::endl;
			file.close();
			exit(-10);
		}

		// Читаем данные pui_Sm
		file.read(reinterpret_cast<char*>(this->pui_Sm.data()),
			size * sizeof(double));

		// Читаем данные pui_Sp
		file.read(reinterpret_cast<char*>(this->pui_Sp.data()),
			rows * cols * sizeof(double));

		file.close();

		// Умножение каждого элемента вектора на n_H_lism
		for (size_t i = 0; i < this->pui_Sm.size(); ++i) 
		{
			this->pui_Sm[i] *= n_H_lism;
		}
		this->pui_Sp = this->pui_Sp * n_H_lism;

	}
	catch (const std::exception& e) {
		std::cerr << "Error 384try78geurgfuegf " << e.what() << std::endl;
		file.close();
		exit(-19);
	}
}

void Cell::read_pui_FromFile(void)
{
	std::string filename = "data_pui/func_cells_pui_" + to_string(this->number) + ".bin";

	if (file_exists(filename) == false) return;

	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open())
	{
		std::cerr << "Error fyjrtyhrtyhert234234asrq3" << filename << std::endl;
		exit(-10);
	}

	try {
		// Читаем размеры матриц
		int size1, size2;
		file.read(reinterpret_cast<char*>(&size1), sizeof(size1));
		file.read(reinterpret_cast<char*>(&size2), sizeof(size2));

		if(size1 != this->f_pui_1.size() || size2 != this->f_pui_2.size())
		{
			std::cerr << "Error ergertgruy567456456" << std::endl;
			cout << size1 << " " << size2 << " " << this->f_pui_1.size() << " " << this->f_pui_2.size() << endl;
			file.close();
			exit(-10);
		}


		file.read(reinterpret_cast<char*>(this->f_pui_1.data()),
			size1 * sizeof(double));

		file.read(reinterpret_cast<char*>(this->f_pui_2.data()),
			size2 * sizeof(double));

		file.close();

	}
	catch (const std::exception& e) {
		std::cerr << "Error ergw4y74567345tfgrgujk89;9 " << e.what() << std::endl;
		file.close();
		exit(-19);
	}
}

double nu_exchenge(const double& u, const double& rho, const double& c, Phys_param*& phys)
{
	// u - модуль разности скоростей 
	if (u / c > 7.0)
	{
		double uz = Velosity_1(u, c);
		return rho * uz * phys->sigma(uz) / phys->par_Kn;
	}
	else
	{
		return (rho * phys->MK_int_1(u, c)) / phys->par_Kn;  // Пробуем вычислять интеграллы численно
	}
}

double Sm_maxwell(const double& w, const double& UHx, const double& UHy, const double& UHz, 
				const double& rhoH, const double& cH, const double& Ux, const double& Uy, const double& Uz, Phys_param*& phys)
{
	// Параметры интегрирования
	const int n_theta = 100;  // Количество шагов по ?
	const int n_phi = 200;    // Количество шагов по ?
	const double r = 1.0;     // Радиус (может быть функцией)

	// Границы интегрирования
	const double theta_min = 0.0;      // от 0
	const double theta_max = const_pi;     // до ?
	const double phi_min = 0.0;        // от 0  
	const double phi_max = 2.0 * const_pi; // до 2?

	// Шаги интегрирования
	const double d_theta = (theta_max - theta_min) / n_theta;
	const double d_phi = (phi_max - phi_min) / n_phi;

	double integral = 0.0;

	// Цикл интегрирования по сферическим углам
	for (int i = 0; i < n_theta; ++i) 
	{
		double theta = theta_min + (i + 0.5) * d_theta;  // середина интервала

		for (int j = 0; j < n_phi; ++j) 
		{
			double phi = phi_min + (j + 0.5) * d_phi;    // середина интервала
			double volume_element = std::sin(theta) * d_theta * d_phi;

			double wx = w * sin(theta) * cos(phi);
			double wy = w * sin(theta) * sin(phi);
			double wz = w * cos(theta);

			double wwx = Ux + wx;
			double wwy = Uy + wy;
			double wwz = Uz + wz;

			// Добавляем вклад в интеграл
			integral += nu_exchenge(norm2(UHx - wwx, UHy - wwy, UHz - wwz), rhoH, cH, phys) * volume_element;
		}
	}

	return integral / (4.0 * const_pi);
}

double Sp_maxwell(const double& w, const double& UHx, const double& UHy, const double& UHz,
	const double& rhoH, const double& cH, const double& Ux, const double& Uy, const double& Uz,
	const double& rho, const double& cp, Phys_param*& phys)
{
	// Параметры интегрирования
	const int n_theta = 100;  // Количество шагов по ?
	const int n_phi = 200;    // Количество шагов по ?
	const double r = 1.0;     // Радиус (может быть функцией)

	// Границы интегрирования
	const double theta_min = 0.0;      // от 0
	const double theta_max = const_pi;     // до ?
	const double phi_min = 0.0;        // от 0  
	const double phi_max = 2.0 * const_pi; // до 2?

	// Шаги интегрирования
	const double d_theta = (theta_max - theta_min) / n_theta;
	const double d_phi = (phi_max - phi_min) / n_phi;

	double integral = 0.0;

	// Цикл интегрирования по сферическим углам
	for (int i = 0; i < n_theta; ++i)
	{
		double theta = theta_min + (i + 0.5) * d_theta;  // середина интервала

		for (int j = 0; j < n_phi; ++j)
		{
			double phi = phi_min + (j + 0.5) * d_phi;    // середина интервала
			double volume_element = std::sin(theta) * d_theta * d_phi;

			double wx = w * sin(theta) * cos(phi);
			double wy = w * sin(theta) * sin(phi);
			double wz = w * cos(theta);

			double wwx = Ux + wx;
			double wwy = Uy + wy;
			double wwz = Uz + wz;

			double f = maxwell(rhoH, cH, UHx, UHy, UHz, wwx, wwy, wwz);
			double nu = nu_exchenge(w, rho, cp, phys);

			// Добавляем вклад в интеграл
			integral += f * nu * volume_element;
		}
	}

	return integral / (4.0 * const_pi);
}


void Cell::print_SmSp(double Wmax, string nam, Phys_param*& phys)
{
	ofstream fout;
	string name_f = "Tecplot_SpSm_" + nam  + "__" + to_string(this->number) + ".txt";
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = u, Sm, Sm_maxwell_H1, Sm_maxwell_H2, Sm_maxwell_H3, Sm_maxwell_H4, Sm_maxwell_all, Sp1, Sp2, Sp_summ, Sp_maxwell_H1, Sp_maxwell_H2, Sp_maxwell_H3, Sp_maxwell_H4, Sp_maxwell" << endl;
	int size = this->pui_Sm.size();
	double dx = Wmax / size;

	double UH1x = this->parameters[0]["MK_Vx_H1"];
	double UH1y = this->parameters[0]["MK_Vy_H1"];
	double UH1z = this->parameters[0]["MK_Vz_H1"];
	double rhoH1 = this->parameters[0]["MK_n_H1"];// *phys->par_n_H_LISM;
	//double cH1 = sqrt(2.0 * this->parameters[0]["MK_p_H1"] / rhoH1);
	double cH1 = sqrt(this->parameters[0]["MK_T_H1"]);

	double UH2x = this->parameters[0]["MK_Vx_H2"];
	double UH2y = this->parameters[0]["MK_Vy_H2"];
	double UH2z = this->parameters[0]["MK_Vz_H2"];
	double rhoH2 = this->parameters[0]["MK_n_H2"];// *phys->par_n_H_LISM;
	//double cH2 = sqrt(2.0 * this->parameters[0]["MK_p_H2"] / rhoH2);
	double cH2 = sqrt(this->parameters[0]["MK_T_H2"]);

	double UH3x = this->parameters[0]["MK_Vx_H3"];
	double UH3y = this->parameters[0]["MK_Vy_H3"];
	double UH3z = this->parameters[0]["MK_Vz_H3"];
	double rhoH3 = this->parameters[0]["MK_n_H3"];// *phys->par_n_H_LISM;
	//double cH3 = sqrt(2.0 * this->parameters[0]["MK_p_H3"] / rhoH3);
	double cH3 = sqrt(this->parameters[0]["MK_T_H3"]);

	double UH4x = this->parameters[0]["MK_Vx_H4"];
	double UH4y = this->parameters[0]["MK_Vy_H4"];
	double UH4z = this->parameters[0]["MK_Vz_H4"];
	double rhoH4 = this->parameters[0]["MK_n_H4"];// *phys->par_n_H_LISM;
	//double cH4 = sqrt(2.0 * this->parameters[0]["MK_p_H4"] / rhoH4);
	double cH4 = sqrt(this->parameters[0]["MK_T_H4"]);

	double Ux = this->parameters[0]["Vx"];
	double Uy = this->parameters[0]["Vy"];
	double Uz = this->parameters[0]["Vz"];
	double rho = this->parameters[0]["rho"];// *phys->par_n_H_LISM;
	double cp = sqrt(2.0 * this->parameters[0]["p"] / rho);


	cout << "rhoH = " << rhoH1 << " " << rhoH2 << " " << rhoH3 << " " << rhoH4 << endl;


	for (int i = 0; i < size; ++i)
	{
		double center = (i + 0.5) * dx;  // центр ячейки - Скорость w

		double Sm1 = Sm_maxwell(center, UH1x, UH1y, UH1z, rhoH1, cH1, Ux, Uy, Uz, phys);
		double Sm2 = Sm_maxwell(center, UH2x, UH2y, UH2z, rhoH2, cH2, Ux, Uy, Uz, phys);
		double Sm3 = Sm_maxwell(center, UH3x, UH3y, UH3z, rhoH3, cH3, Ux, Uy, Uz, phys);
		double Sm4 = Sm_maxwell(center, UH4x, UH4y, UH4z, rhoH4, cH4, Ux, Uy, Uz, phys);
		double Sm = Sm1 + Sm2 + Sm3 + Sm4;

		double Sp1 = Sp_maxwell(center, UH1x, UH1y, UH1z, rhoH1, cH1, Ux, Uy, Uz, rho, cp, phys);
		double Sp2 = Sp_maxwell(center, UH2x, UH2y, UH2z, rhoH2, cH2, Ux, Uy, Uz, rho, cp, phys);
		double Sp3 = Sp_maxwell(center, UH3x, UH3y, UH3z, rhoH3, cH3, Ux, Uy, Uz, rho, cp, phys);
		double Sp4 = Sp_maxwell(center, UH4x, UH4y, UH4z, rhoH4, cH4, Ux, Uy, Uz, rho, cp, phys);
		double Sp = Sp1 + Sp2 + Sp3 + Sp4;


		fout << center << " " << this->pui_Sm[i] << " " << Sm1 << " " << Sm2 << " " << Sm3 << " " << Sm4 << " " << Sm << " "
			<< this->pui_Sp(0, i) << " " << this->pui_Sp(1, i) << " " << this->pui_Sp(0, i) + this->pui_Sp(1, i)  << " " 
			<< Sp1 << " " << Sp2 << " " << Sp3 << " " << Sp4 << " " << Sp << std::endl;
	}

	fout.close();
}

void Cell::print_pui(double Wmax, string nam)
{
	ofstream fout;
	string name_f = "Tecplot_pui_" + nam + "__" + to_string(this->number) + ".txt";
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = u, f1, f2, f" << endl;
	int size = this->f_pui_1.size();
	double dx = Wmax / size;
	for (int i = 0; i < size; ++i)
	{
		double center = (i + 0.5) * dx;  // центр ячейки
		double f2 = 0.0;
		if (this->f_pui_2.size() > i) f2 = this->f_pui_2[i];
		fout << center << " " << this->f_pui_1[i] << " " << f2 << " " <<
			this->f_pui_1[i] + f2 << std::endl;
	}

	fout.close();

	fout.open("interpol_" + name_f);
	fout << "TITLE = HP  VARIABLES = u, f1, f2, f" << endl;
	size = this->f_pui_1.size();
	dx = Wmax / ( 9 * size);
	for (int i = 0; i < 10 * size; ++i)
	{
		double center = (i + 0.5) * dx;  // центр ячейки
		double f2 = 0.0;
		if (this->f_pui_2.size() > 0) f2 = this->pui_get_f(center, 1, Wmax);
		fout << center << " " << this->pui_get_f(center, 0, Wmax) << " " << f2 << " " <<
			this->pui_get_f(center, 0, Wmax) + f2 << std::endl;
	}

	fout.close();
}

void Cell::culc_pui_n_T(const double& pui_wR)
{
	double S = 0.0;
	double S2 = 0.0;

	double rho = this->parameters[0]["rho"];
	double rho_He = this->parameters[0]["rho_He"];

	while (true)
	{
		S = 0.0;
		S2 = 0.0;

		short int pui_nw = this->f_pui_1.size();
		for (short int i = 0; i < this->f_pui_1.size(); i++)
		{
			double w = (i + 0.5) * pui_wR / this->f_pui_1.size();
			S = S + this->f_pui_1[i] * 4.0 * const_pi * kv(w) * (pui_wR / pui_nw);
			S2 = S2 + this->f_pui_1[i] * 4.0 * const_pi * pow4(w) * (pui_wR / pui_nw);
		}

		if (S > 0.0000001)
		{
			S2 = S2 / (S * 3.0);
		}
		else
		{
			S2 = 1.0;
		}
		this->parameters[0]["MK_rho_Pui_1"] = S;
		this->parameters[0]["MK_T_Pui_1"] = S2;


		pui_nw = this->f_pui_2.size();
		S = S2 = 0.0;

		for (short int i = 0; i < this->f_pui_2.size(); i++)
		{
			double w = (i + 0.5) * pui_wR / this->f_pui_2.size();
			S = S + this->f_pui_2[i] * 4.0 * const_pi * kv(w) * (pui_wR / pui_nw);
			S2 = S2 + this->f_pui_2[i] * 4.0 * const_pi * pow4(w) * (pui_wR / pui_nw);
		}

		if (S > 0.0000001)
		{
			S2 = S2 / (S * 3.0);
		}
		else
		{
			S2 = 1.0;
		}
		this->parameters[0]["MK_rho_Pui_2"] = S;
		this->parameters[0]["MK_T_Pui_2"] = S2;

		if (this->parameters[0]["MK_rho_Pui_1"] + this->parameters[0]["MK_rho_Pui_2"] > (rho - rho_He) * 0.995)
		{
			double kkl = ((rho - rho_He) * 0.99499) / 
				(this->parameters[0]["MK_rho_Pui_1"] + this->parameters[0]["MK_rho_Pui_2"]);

			if (kkl > 1.0)
			{
				cout << "Error dgiuehrfguiheiprhf934rrr   " << kkl <<  endl;
				exit(-1);
			}
			for (short int i = 0; i < this->f_pui_1.size(); i++)
			{
				this->f_pui_1[i] *= kkl;
			}

			for (short int i = 0; i < this->f_pui_2.size(); i++)
			{
				this->f_pui_2[i] *= kkl;
			}
		}
		else
		{
			break;
		}
	}




}


double Cell::pui_get_nu(const double& w, short int ii, const double& Wmax)
{
	// Выбор нужного вектора
	const auto& vec = (ii == 0) ? this->nu_integr_pui_1 : (ii == 1) ? this->nu_integr_pui_2 : [&]() -> const vector<double>&
		{
			cout << "ERROR ey4556uy564556ty4y453dfsde" << endl;
			exit(-1);
			return f_pui_1; // заглушка, никогда не выполнится
		}();

	const short int N = vec.size();
	if (N == 0) return 0.0;

	// Проверка граничных значений
	if (w <= 0.0) return vec[0];
	if (w >= Wmax) return vec[N - 1];

	const double cell_size = Wmax / N;
	int left_index = static_cast<int>(w / cell_size);

	// Корректировка индексов
	if (w < (left_index + 0.5) * cell_size)
	{
		left_index--;
	}
	int right_index = left_index + 1;

	// Проверка скорректированных индексов
	if (left_index < 0) return vec[0];
	if (right_index >= N) return vec[N - 1];

	// Линейная интерполяция
	const double left_center = (left_index + 0.5) * cell_size;
	const double right_center = (right_index + 0.5) * cell_size;
	const double t = (w - left_center) / (right_center - left_center);

	return max(vec[left_index] * (1.0 - t) + vec[right_index] * t, 0.0);
}

double Cell::pui_get_f(const double& w, short int ii, const double& Wmax)
{
	// Выбор нужного вектора
	const auto& vec = (ii == 0) ? this->f_pui_1 : (ii == 1) ? this->f_pui_2 : [&]() -> const vector<double>&
		{
		cout << "ERROR ey4556uy564556ty4y453dfsde" << endl;
		exit(-1);
		return f_pui_1; // заглушка, никогда не выполнится
		}();

	const short int N = vec.size();
	if (N == 0) return 0.0;

	// Проверка граничных значений
	if (w <= 0.0) return vec[0];
	if (w >= Wmax) return vec[N - 1];

	const double cell_size = Wmax / N;
	int left_index = static_cast<int>(w / cell_size);

	// Корректировка индексов
	if (w < (left_index + 0.5) * cell_size) 
	{
		left_index--;
	}
	int right_index = left_index + 1;

	// Проверка скорректированных индексов
	if (left_index < 0) return vec[0];
	if (right_index >= N) return vec[N - 1];

	// Линейная интерполяция
	const double left_center = (left_index + 0.5) * cell_size;
	const double right_center = (right_index + 0.5) * cell_size;
	const double t = (w - left_center) / (right_center - left_center);

	return max(vec[left_index] * (1.0 - t) + vec[right_index] * t, 0.0);
}

double Cell::PUI_get_F_integer(const double& w, short int ii)
{
	const double Wmax = 1.0;
	// Выбор нужного вектора
	const auto& vec = (ii == 0) ? this->F_integr_pui_1 : (ii == 1) ? this->F_integr_pui_2 : [&]() -> const vector<double>&
		{
			cout << "ERROR ey4556uy564556ty4y453dfsde" << endl;
			exit(-1);
			return this->F_integr_pui_1; // заглушка, никогда не выполнится
		}();

	const short int N = vec.size();
	if (N == 0) return 0.0;

	// Проверка граничных значений
	if (w <= 0.0) return vec[0];
	if (w >= Wmax) return vec[N - 1];

	const double cell_size = Wmax / N;
	int left_index = static_cast<int>(w / cell_size);

	// Корректировка индексов
	if (w < (left_index + 0.5) * cell_size)
	{
		left_index--;
	}
	int right_index = left_index + 1;

	// Проверка скорректированных индексов
	if (left_index < 0) return vec[0];
	if (right_index >= N) return vec[N - 1];

	// Линейная интерполяция
	const double left_center = (left_index + 0.5) * cell_size;
	const double right_center = (right_index + 0.5) * cell_size;
	const double t = (w - left_center) / (right_center - left_center);

	return max(vec[left_index] * (1.0 - t) + vec[right_index] * t, 0.0);
}

void Cell::MK_pui_charge_exchange_velocity(Sensor* sens, Setka* SS, Phys_param* Phys,
	const double& Upx, const double& Upy, 
	const double& Upz, const double& UHx, const double& UHy, const double& UHz, 
	double& VHx, double& VHy, double& VHz, short int num_pui)
{
	// !? Функция перезарядки - вычисляет новые скорости атома после перезарядки

	double ksi1, ksi2, ksi3;
	double UH, w, the, u, h0, phi;
	double vx, vy, vz;

	Eigen::Vector3d ex, ey, ez;

	ez[0] = UHx - Upx;
	ez[1] = UHy - Upy;
	ez[2] = UHz - Upz;

	UH = sqrt(kv(Upx - UHx) + kv(Upy - UHy) + kv(Upz - UHz));
	ez = ez / UH;
	h0 = SS->PUI_get_h0(UH);

	get_bazis(ez, ex, ey);

	while (true)
	{
		ksi1 = sens->MakeRandom();
		ksi2 = sens->MakeRandom();
		ksi3 = sens->MakeRandom();

		w = this->PUI_get_F_integer(ksi1, num_pui);
		the = acos(1.0 - 2.0 * ksi2);
		u = sqrt(kv(w) * kv(sin(the)) + kv(w * cos(the) - UH));
		if (u * Phys->sigma(u) / ((w + Phys->pui_h0_wc) * Phys->sigma(w + Phys->pui_h0_wc) * h0) >= ksi3) break;
	}

	ksi1 = sens->MakeRandom();
	phi = ksi1 * 2.0 * const_pi;

	vx = w * sin(the) * cos(phi);
	vy = w * sin(the) * sin(phi);
	vz = w * cos(the);

	VHx = Upx + vx * ex(0) + vy * ey(0) + vz * ez(0);
	VHy = Upy + vx * ex(1) + vy * ey(1) + vz * ez(1);
	VHz = Upz + vx * ex(2) + vy * ey(2) + vz * ez(2);
}

void Cell::Get_RBF_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par)
{
	Eigen::Vector3d point;
	Eigen::Vector3d point2;

	point << x, y, z;

	if (this->interpolate_alpha.find("rho") != this->interpolate_alpha.end())
	{
		double result = 0.0;
		for (size_t i = 0; i < this->yzels.size(); ++i) 
		{
			point2 << this->yzels[i]->coord[0][0], this->yzels[i]->coord[0][1],
				this->yzels[i]->coord[0][2];
			result += this->interpolate_alpha["rho"][i] * rbfKernel((point - point2).norm());
		}
		par["rho"] = result;
	}
}

void Cell::Get_IDW_interpolation(const double& x, const double& y, const double& z, 
	unordered_map<string, double>& par, Phys_param* phys_param)
{
	Eigen::Vector3d point;
	Eigen::Vector3d point2;

	point << x, y, z;

	unordered_map<string, double> sum_weights;
	unordered_map<string, double> sum_weighted_values;

	for (auto& nam : phys_param->param_names)
	{
		sum_weights[nam] = 0.0;
		sum_weighted_values[nam] = 0.0;
	}

	//double sum_weights = 0.0;
	//double sum_weighted_values = 0.0;

	for (size_t i = 0; i < this->yzels.size(); ++i)
	{
			
		point2 << this->yzels[i]->coord[0][0], 
			this->yzels[i]->coord[0][1], this->yzels[i]->coord[0][2];
		double dist = 0.0001 + (point - point2).norm();
		double weight = 1.0 / std::pow(dist, 1.0);

		for (auto& nam : phys_param->param_names)
		{
			if (this->yzels[i]->parameters.find(nam) != this->yzels[i]->parameters.end())
			{

				sum_weights[nam] += weight;
				sum_weighted_values[nam] += weight * this->yzels[i]->parameters[nam];
				//cout << this->yzels[i]->parameters["rho"] << endl;
			}
		}
	}


	// Добавим центр ячейки
	point2 << this->center[0][0],
		this->center[0][1], this->center[0][2];
	double dist = 0.0001 + (point - point2).norm();
	double weight = 1.0 / std::pow(dist, 1.0);

	for (auto& nam : phys_param->param_names)
	{
		if (this->parameters[0].find(nam) != this->parameters[0].end())
		{

			sum_weights[nam] += weight;
			sum_weighted_values[nam] += weight * this->parameters[0][nam];
			//cout << this->yzels[i]->parameters["rho"] << endl;
		}
	}

	for (auto& nam : phys_param->param_names)
	{
		if (sum_weights[nam] > 0.00001)
		{
			par[nam] = sum_weighted_values[nam] / sum_weights[nam];
		}
		else
		{
			par[nam] = 0.0;
		}
	}
}

Cell* Cell::Get_Sosed(Gran* gr)
{
	if (gr == nullptr)
	{
		cout << "Error 5434323121" << endl;
		exit(-1);
	}

	if(gr->cells.size() == 1) return nullptr;

	if (gr->cells[0]->number == this->number)
	{
		return gr->cells[1];
	}
	else
	{
		return gr->cells[0];
	}
	return nullptr;
}

double determinant4x4(const Eigen::Vector4d& a, const Eigen::Vector4d& b, 
	const Eigen::Vector4d& c, const Eigen::Vector4d& d) 
{
	Eigen::Matrix4d mat;
	mat << a, b, c, d;
	return mat.determinant();
}

bool isPointInsideTetrahedron(const Eigen::Vector3d& P, const Eigen::Vector3d& A, const Eigen::Vector3d& B, 
	const Eigen::Vector3d& C, const Eigen::Vector3d& D) 
{
	// Добавляем 4-ю координату (1) для работы с 4?4 матрицами
	Eigen::Vector4d A4(A(0), A(1), A(2), 1.0);
	Eigen::Vector4d B4(B(0), B(1), B(2), 1.0);
	Eigen::Vector4d C4(C(0), C(1), C(2), 1.0);
	Eigen::Vector4d D4(D(0), D(1), D(2), 1.0);
	Eigen::Vector4d P4(P(0), P(1), P(2), 1.0);

	double detABCD = determinant4x4(A4, B4, C4, D4);
	if (std::abs(detABCD) < 1e-10) {
		std::cerr << "Degenerate tetrahedron (zero volume)!" << std::endl;
		cout << A4(0) << " " << A4(1) << " " << A4(2) << endl;
		cout << B4(0) << " " << B4(1) << " " << B4(2) << endl;
		cout << C4(0) << " " << C4(1) << " " << C4(2) << endl;
		cout << D4(0) << " " << D4(1) << " " << D4(2) << endl;
		cout << "__________________________________________" << endl;
		return false;
	}

	double u = determinant4x4(P4, B4, C4, D4) / detABCD;
	if (u < -1e-10) return false;
	double v = determinant4x4(A4, P4, C4, D4) / detABCD;
	if (v < -1e-10) return false;
	double w = determinant4x4(A4, B4, P4, D4) / detABCD;
	if (w < -1e-10) return false;
	double t = determinant4x4(A4, B4, C4, P4) / detABCD;
	if (t < -1e-10) return false;

	return (std::fabs(u + v + w + t - 1.0) < 1e-6);
}


bool Cell::Belong_point(const double& x, const double& y, const double& z, short int now, bool fast, Cell*& Next)
{
	// fast = false - правильный алгоритм
	// fast = true - быстрый алгоритм (но имеет погрешность)

	Eigen::Vector3d P(x, y, z); // Точка
	Next = nullptr;  // Next - следующая ячейка определяется только если fast == true
	// алгоритм для fast == false не позволяет найти следующую ячейку

	//double xc = this->center[now][0];
	//double yc = this->center[now][1];
	//double zc = this->center[now][2];

	// Честное определение принадлежности
	if (fast == false)
	{
		Eigen::Vector3d A(this->center[now][0], this->center[now][1], this->center[now][2]);
		Eigen::Vector3d B, C, D;
		bool is_in = false;
		for (auto& gr : this->grans)
		{
			//double x1 = gr->center[now][0];
			//double y1 = gr->center[now][1];
			//double z1 = gr->center[now][2];
			if (false) // старый вариант, где грань = 4 треугольника
			{
				B << gr->center[now][0], gr->center[now][1], gr->center[now][2];
				int n_yz = gr->yzels.size();
				for (short int i = 0; i < n_yz; i++)
				{
					int j = i + 1;
					if (j >= n_yz) j = 0;
					C << gr->yzels[i]->coord[now][0], gr->yzels[i]->coord[now][1], gr->yzels[i]->coord[now][2];
					D << gr->yzels[j]->coord[now][0], gr->yzels[j]->coord[now][1], gr->yzels[j]->coord[now][2];
					is_in = isPointInsideTetrahedron(P, A, B, C, D);
					if (is_in == true) return true;
				}
			}
			else // грань = 2 треугольника
			{
				if (gr->yzels.size() != 4)
				{
					cout << "Error 9765486573" << endl;
					exit(-1);
				}

				B << gr->yzels[0]->coord[now][0], gr->yzels[0]->coord[now][1], gr->yzels[0]->coord[now][2];
				C << gr->yzels[1]->coord[now][0], gr->yzels[1]->coord[now][1], gr->yzels[1]->coord[now][2];
				D << gr->yzels[2]->coord[now][0], gr->yzels[2]->coord[now][1], gr->yzels[2]->coord[now][2];
				is_in = isPointInsideTetrahedron(P, A, B, C, D);
				if (is_in == true) return true;

				
				C << gr->yzels[3]->coord[now][0], gr->yzels[3]->coord[now][1], gr->yzels[3]->coord[now][2];
				is_in = isPointInsideTetrahedron(P, A, B, C, D);
				if (is_in == true) return true;

			}
		}

		return false;
	}
	else
	{
		Eigen::Vector3d B, C, D;

		for (auto& gr : this->grans)
		{
			if (gr->cells.size() != 2) continue;

			B << gr->center[now][0], gr->center[now][1], gr->center[now][2];
			C << gr->normal[now][0], gr->normal[now][1], gr->normal[now][2];

			if (gr->cells[0]->number == this->number) C = -C;

			D = P - B;

			if (C.dot(D) < -1e-10) 
			{
				Next = this->Get_Sosed(gr);
				return false;
			}
		}
		return true;
	}

	return false;
}

void Cell::Set_Cell_Geo_for_MK(void)
{
	Eigen::Vector3d A;
	Eigen::Vector3d B;
	Eigen::Vector3d C;


	// Считаем l_size   характерный размер ячейки
	if (true)
	{
		double S = 0.0;

		A[0] = this->yzels[0]->coord[0][0];
		A[1] = this->yzels[0]->coord[0][1];
		A[2] = this->yzels[0]->coord[0][2];

		B[0] = this->yzels[0]->coord[1][0];
		B[1] = this->yzels[0]->coord[1][1];
		B[2] = this->yzels[0]->coord[1][2];
		C = B - A;

		S = C.norm();

		B[0] = this->yzels[0]->coord[4][0];
		B[1] = this->yzels[0]->coord[4][1];
		B[2] = this->yzels[0]->coord[4][2];
		C = B - A;
		S = min(S, C.norm());


		B[0] = this->yzels[0]->coord[2][0];
		B[1] = this->yzels[0]->coord[2][1];
		B[2] = this->yzels[0]->coord[2][2];
		C = B - A;
		S = min(S, C.norm());

		this->geo_parameters["l_size"] = S;
	}
}


void Cell::Culc_center(unsigned short int st_time)
{
	double xc, yc, zc;
	int Ny = this->yzels.size();
	xc = 0.0;
	yc = 0.0;
	zc = 0.0;

	// Вычисляем центр грани
	for (auto& i : this->yzels)
	{
		xc += i->coord[st_time][0];
		yc += i->coord[st_time][1];
		zc += i->coord[st_time][2];
	}
	xc /= Ny;
	yc /= Ny;
	zc /= Ny;
	this->center[st_time][0] = xc;
	this->center[st_time][1] = yc;
	this->center[st_time][2] = zc;
}

void Cell::Culc_volume(unsigned short int st_time, unsigned short int method)
{
	// 0 - быстрый вариант (нужны посчанные площади граней и их центры)
	// 1 - медленный вариант (нужны посчитанные только центры граней)
	// для ровных граней оба методы работают одинаково, но чем более кривая грань, 
	// тем большее расхождения получается.
	// медленный вариант более правильный (но это не точно). 
	// Думаю всегда можно выбирать вариант "0"
	if (method == 0)
	{
		double xg, yg, zg, h;
		double V = 0.0;
		for (auto& i : this->grans)
		{
			xg = i->center[st_time][0];
			yg = i->center[st_time][1];
			zg = i->center[st_time][2];
			h = fabs(scalarProductFast(i->normal[st_time][0], i->normal[st_time][1],
				i->normal[st_time][2],
				this->center[st_time][0] - xg, this->center[st_time][1] - yg,
				this->center[st_time][2] - zg));
			V += h * i->area[st_time];
		}
		this->volume[st_time] = V / 3.0; 
		return;
	}
	if (method == 1)
	{
		int j1;
		double V = 0.0;
		for (auto& i : this->grans)
		{
			if (false) // старый вариант где грань = 4 треугольника
			{
				for (int j = 0; j < i->yzels.size(); j++)
				{
					j1 = j + 1;
					if (j1 >= i->yzels.size()) j1 = 0;
					V += tetrahedronVolume(this->center[st_time][0], this->center[st_time][1],
						this->center[st_time][2],
						i->center[st_time][0], i->center[st_time][1], i->center[st_time][2],
						i->yzels[j]->coord[st_time][0], i->yzels[j]->coord[st_time][1],
						i->yzels[j]->coord[st_time][2],
						i->yzels[j1]->coord[st_time][0], i->yzels[j1]->coord[st_time][1],
						i->yzels[j1]->coord[st_time][2]);
				}
			}
			else
			{
				V += tetrahedronVolume(this->center[st_time][0], this->center[st_time][1],
					this->center[st_time][2],
					i->yzels[0]->coord[st_time][0], i->yzels[0]->coord[st_time][1],
					i->yzels[0]->coord[st_time][2],
					i->yzels[1]->coord[st_time][0], i->yzels[1]->coord[st_time][1],
					i->yzels[1]->coord[st_time][2],
					i->yzels[2]->coord[st_time][0], i->yzels[2]->coord[st_time][1],
					i->yzels[2]->coord[st_time][2]);
				V += tetrahedronVolume(this->center[st_time][0], this->center[st_time][1],
					this->center[st_time][2],
					i->yzels[0]->coord[st_time][0], i->yzels[0]->coord[st_time][1],
					i->yzels[0]->coord[st_time][2],
					i->yzels[3]->coord[st_time][0], i->yzels[3]->coord[st_time][1],
					i->yzels[3]->coord[st_time][2],
					i->yzels[2]->coord[st_time][0], i->yzels[2]->coord[st_time][1],
					i->yzels[2]->coord[st_time][2]);
			}
		}
		this->volume[st_time] = V;
		return;
	}
	else
	{
		cout << "Error  0967452122" << endl;
		return;
	}
}

double Cell ::func_R(unsigned short int i_time)
{
	return norm2(this->center[i_time][0], this->center[i_time][1], this->center[i_time][2]);
}

void Cell::Tecplot_print_cell(void)
{
	// name - это имя лучей
	ofstream fout;
	string name_f = "Tecplot_cell_in_3D_.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	fout << "ZONE T=HP, N = " << this->yzels.size() << ", E = " << this->yzels.size() - 1 << ", F=FEPOINT, ET=LINESEG" << endl;

	for (auto& j : this->yzels)
	{
		fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
	}

	for (int m = 0; m < this->yzels.size() - 1; m++)
	{
		fout << m + 1 << " " << m + 2 << endl;
	}

	fout << endl;


	fout.close();
}

void Cell::MK_Add_moment(MK_particle& P, const double& cp, const double& u, const double& mu_ex,
	const double& u1, const double& u2, const double& u3, const double& skalar, Phys_param* phys_param)
{
	string name_H;

	// Определяем сорт водорода для записи моментов
	switch (P.sort)
	{
	case 1:
		name_H = "H1";
		break;
	case 2:
		name_H = "H2";
		break;
	case 3:
		name_H = "H3";
		break;
	case 4:
		name_H = "H4";
		break;
	case 5:
		name_H = "H5";
		break;
	case 6:
		name_H = "H6";
		break;
	case 7:
		name_H = "H7";
		break;
	case 8:
		name_H = "H8";
		break;
	case 9:
		name_H = "H9";
		break;
	case 10:
		name_H = "H10";
		break;
	default:
		cout << "ERROR NO SUCH MOMENT j9egrhg9u34980tuf9hwe9prggfewr" << endl;
		exit(-1);
		break;
	}


	if (u / cp > 7.0)
	{
		double uz = Velosity_1(u, cp);
		double uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * const_pi * sqrt_pi);
		double uz_E = Velosity_3(u, cp);

		this->mut.lock();
		this->parameters[0]["MK_IVx_H"] -= mu_ex * uz_M * u1 / u;
		this->parameters[0]["MK_IVy_H"] -= mu_ex * uz_M * u2 / u;
		this->parameters[0]["MK_IVz_H"] -= mu_ex * uz_M * u3 / u;
		this->parameters[0]["MK_IT_H"] += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - uz_M * skalar / u);

		this->parameters[0]["MK_IVx_" + name_H] -= mu_ex * uz_M * u1 / u;
		this->parameters[0]["MK_IVy_" + name_H] -= mu_ex * uz_M * u2 / u;
		this->parameters[0]["MK_IVz_" + name_H] -= mu_ex * uz_M * u3 / u;
		this->parameters[0]["MK_IT_" + name_H] += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - uz_M * skalar / u);

		this->mut.unlock();
	}
	else
	{
		double k1 = phys_param->MK_int_1(u, cp);
		double k2 = phys_param->MK_int_2(u, cp);
		double k3 = phys_param->MK_int_3(u, cp);

		this->mut.lock();
		this->parameters[0]["MK_IVx_H"] += mu_ex * (k2 / k1) * u1 / u;
		this->parameters[0]["MK_IVy_H"] += mu_ex * (k2 / k1) * u2 / u;
		this->parameters[0]["MK_IVz_H"] += mu_ex * (k2 / k1) * u3 / u;
		this->parameters[0]["MK_IT_H"] += mu_ex * (-0.5 * k3 / k1 + k2 / k1 * skalar / u);

		this->parameters[0]["MK_IVx_" + name_H] += mu_ex * (k2 / k1) * u1 / u;
		this->parameters[0]["MK_IVy_" + name_H] += mu_ex * (k2 / k1) * u2 / u;
		this->parameters[0]["MK_IVz_" + name_H] += mu_ex * (k2 / k1) * u3 / u;
		this->parameters[0]["MK_IT_" + name_H] += mu_ex * (-0.5 * k3 / k1 + k2 / k1 * skalar / u);
		this->mut.unlock();
	}

}

void Cell::MK_Add_particle(MK_particle& P, const double& time, Phys_param* phys_param)
{
	// Поглощение
	if (this->mas_pogl.size() > 0)
	{
		double R = sqrt(kv(this->center[0][0]) +
			kv(this->center[0][1]) +
			kv(this->center[0][2]));
		double Ve = (P.Vel[0] * this->center[0][0] +
			P.Vel[1] * this->center[0][1] +
			P.Vel[2] * this->center[0][2]) / R;
		int ii = pogl_mas_number(Ve, phys_param->pogl_L, phys_param->pogl_R, phys_param->pogl_n);
		this->mut.lock();
		this->mas_pogl(P.sort - 1, ii) += time * P.mu;
		this->mut.unlock();
	}


	string name_H;

	// Определяем сорт водорода для записи моментов
	switch (P.sort)
	{
	case 1:
		name_H = "H1";
		break;
	case 2:
		name_H = "H2";
		break;
	case 3:
		name_H = "H3";
		break;
	case 4:
		name_H = "H4";
		break;
	case 5:
		name_H = "H5";
		break;
	case 6:
		name_H = "H6";
		break;
	case 7:
		name_H = "H7";
		break;
	case 8:
		name_H = "H8";
		break;
	case 9:
		name_H = "H9";
		break;
	case 10:
		name_H = "H10";
		break;
	default:
		cout << "ERROR NO SUCH MOMENT j9egrhg9u34980tuf9hwe9prggfewr" << endl;
		exit(-1);
		break;
	}

	this->mut.lock();

	this->parameters[0]["MK_n_H"] += time * P.mu;

	this->parameters[0]["MK_n_" + name_H] += time * P.mu;
	this->parameters[0]["MK_Vx_" + name_H] += time * P.mu * P.Vel[0];
	this->parameters[0]["MK_Vy_" + name_H] += time * P.mu * P.Vel[1];
	this->parameters[0]["MK_Vz_" + name_H] += time * P.mu * P.Vel[2];
	this->parameters[0]["MK_T_" + name_H] += time * P.mu * kvv(P.Vel[0], P.Vel[1], P.Vel[2]);

	this->mut.unlock();
}

void Cell::MK_Add_pui_source(MK_particle& P, const double& wr, const double& nu_ex, const double& mu,
	const double& time, Phys_param* phys_param, short int zone, short int parent)
{
	// zone = 1, 2, 3, 4
	// parent - 0, 1, 2 от каго рождён атом? тепловой протон, pui1, pui2

	int index = static_cast<int>(wr / phys_param->pui_wR * phys_param->pui_nW);
	if (index < 0)  index = 0;
	if (index >= phys_param->pui_nW)  index = phys_param->pui_nW - 1;

	this->mut.lock();
	this->pui_Sm[index] += mu * time;
	this->mut.unlock();

	// Надо понять, в какой S записываем 
	short int k = 0;
	if (zone == 1)
	{
		k = phys_param->proton_arise_1(P.sort - 1, parent);
	}
	else if (zone == 2)
	{
		k = phys_param->proton_arise_2(P.sort - 1, parent);
	}
	else if (zone == 3)
	{
		k = phys_param->proton_arise_3(P.sort - 1, parent);
	}
	else if (zone == 4)
	{
		k = phys_param->proton_arise_4(P.sort - 1, parent);
	}
	else
	{
		cout << "ERROR j9egrhg9u34980tuf9hwe9prggfewr" << endl;
		exit(-1);
	}

	//cout << "Proverka:   " <<  P.sort << "   " << parent << "   " << zone << "   " << k << endl;


	if (k == 0) return;
	k = k - 1;
	
	this->mut.lock();
	this->pui_Sp(k, index) +=  mu * time * nu_ex;
	this->mut.unlock();
}

void Cell::MK_calc_Sm(Phys_param* phys_param)
{
	vector<double> pui_Sm2(phys_param->pui_nW);
	for (auto& i : pui_Sm2) i = 0.0;
	double dthe = const_pi / 40.0;
	double Vh, ff, the, d;

	for (size_t ij = 0; ij < phys_param->pui_nW; ij++)
	{
		double w = ((ij + 0.5) * phys_param->pui_wR / phys_param->pui_nW);
		for (size_t j = 0; j < phys_param->pui_nW; j++)
		{
			Vh = ((j + 0.5) * phys_param->pui_wR / phys_param->pui_nW);
			ff = this->pui_Sm[j];
			if (ff <= 0.0) continue;
			for (size_t k = 1; k <= 40; k++)
			{
				the = dthe * k;
				d = sqrt(kv(Vh) + kv(w) - 2.0 * w * Vh * cos(the));
				if (d > 0.000000001)
				{
					pui_Sm2[ij] += ff * d * phys_param->sigma(d) * sin(the) * dthe * 2.0 * const_pi;
				}
			}
		}
		pui_Sm2[ij] = pui_Sm2[ij] / (4.0 * const_pi);
	}

	for (size_t ij = 0; ij < phys_param->pui_nW; ij++)
	{
		this->pui_Sm[ij] = pui_Sm2[ij] / phys_param->par_Kn;// *phys_param->par_n_H_LISM;
		// Я решил не домножать источники на концентрацию атомов. Нужно на неё домжить непосредственно при использовании
	}

}

void Cell::MK_normir_Moments(Phys_param* phys_param)
{

	if (this->mas_pogl.size() > 0)
	{
		this->mas_pogl /= this->volume[0];
	}

	// Общие моменты
	if (this->parameters[0].find("MK_n_H") != this->parameters[0].end())
	{
		this->parameters[0]["MK_n_H"] /= this->volume[0];
	}

	if (this->parameters[0].find("MK_IVx_H") != this->parameters[0].end())
	{
		this->parameters[0]["MK_IVx_H"] *= (phys_param->par_n_H_LISM/this->volume[0]);
		this->parameters[0]["MK_IVy_H"] *= (phys_param->par_n_H_LISM/this->volume[0]);
		this->parameters[0]["MK_IVz_H"] *= (phys_param->par_n_H_LISM/this->volume[0]);
		this->parameters[0]["MK_IT_H"] *= (phys_param->par_n_H_LISM/this->volume[0]);
	}

	// Моменты по сортам водорода
	vector<string> names_H = { "H1","H2","H3","H4","H5","H6","H7","H8","H9","H10" };
	for (const auto& name : names_H)
	{
		if (this->parameters[0].find("MK_n_" + name) != this->parameters[0].end())
		{
			if (this->parameters[0]["MK_n_" + name] > 0.0000000001)
			{
				this->parameters[0]["MK_Vx_" + name] /= (this->parameters[0]["MK_n_" + name]);
				this->parameters[0]["MK_Vy_" + name] /= (this->parameters[0]["MK_n_" + name]);
				this->parameters[0]["MK_Vz_" + name] /= (this->parameters[0]["MK_n_" + name]);
				this->parameters[0]["MK_T_" + name] = (2.0 / 3.0) * (this->parameters[0]["MK_T_" + name] / this->parameters[0]["MK_n_" + name] -
					kvv(this->parameters[0]["MK_Vx_" + name], this->parameters[0]["MK_Vy_" + name], this->parameters[0]["MK_Vz_" + name]));
			}
			else
			{
				this->parameters[0]["MK_Vx_" + name] = 0.0;
				this->parameters[0]["MK_Vy_" + name] = 0.0;
				this->parameters[0]["MK_Vz_" + name] = 0.0;
				this->parameters[0]["MK_T_" + name] = 0.0;
			}

			this->parameters[0]["MK_n_" + name] /= this->volume[0];
		}
	}


	if (phys_param->MK_source_S == true)
	{
		for (size_t ij = 0; ij < phys_param->pui_nW; ij++)
		{
			double pui_w1 = ij * phys_param->pui_wR / phys_param->pui_nW;
			double pui_w2 = (ij + 1) * phys_param->pui_wR / phys_param->pui_nW;
			this->pui_Sm[ij] /= (this->volume[0]);
			for (size_t ki = 0; ki < this->pui_Sp.rows(); ki++)
			{
				if (ki > 1)
				{
					cout << "Error lokal 9803u45th9erfgre" << endl;
				}
				this->pui_Sp(ki, ij) /= (this->volume[0] * 4.0 * const_pi * (1.0 / 3.0) * (pow3(pui_w2) - pow3(pui_w1)));
			}
		}
	}
}
