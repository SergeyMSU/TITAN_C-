#include "Setka.h"
#include <algorithm>
#include <filesystem> // Для работы с файловой системой

namespace fs = std::filesystem; // Создаем псевдоним для удобства
using namespace Eigen;

/**
 * @brief Решает уравнение Лапласа для магнитного потенциала методом фиктивных источников (MFS)
 *
 * @param fict_points Координаты фиктивных источников (N x 3)
 * @param bnd_points Координаты граничных точек (M x 3)
 * @param bnd_normals Векторы нормалей в граничных точках (M x 3)
 * @param bnd_Bn Значения нормальной компоненты Bn в граничных точках (M)
 * @param use_svd Использовать SVD (true) или QR (false) для решения. SVD устойчивее к плохой обусловленности.
 * @return VectorXd Вектор амплитуд источников q (N)
 */
VectorXd solveMFS(const MatrixXd& fict_points,    // N x 3
	const MatrixXd& bnd_points,     // M x 3
	const MatrixXd& bnd_normals,    // M x 3
	const VectorXd& bnd_Bn,         // M
	bool use_svd = true) {         // Переключатель метода решения

	const int N = fict_points.rows();  // Количество фиктивных источников
	const int M = bnd_points.rows();   // Количество граничных точек

	// Проверка размеров
	assert(fict_points.cols() == 3 && bnd_points.cols() == 3 && bnd_normals.cols() == 3);
	assert(bnd_normals.rows() == M && bnd_Bn.size() == M);

	// Матрица системы: M x N (плюс одна строка для условия однозначности)
	MatrixXd A = MatrixXd::Zero(M + 1, N);
	VectorXd b = VectorXd::Zero(M + 1);

	cout << "Zapolnyaem " << endl;
	// Заполняем основную часть матрицы A и вектора b
	const double pi4 = 4.0 * M_PI;
	for (int i = 0; i < M; i++) 
	{
		// Нормаль в i-й граничной точке
		Vector3d ni = bnd_normals.row(i);

		for (int j = 0; j < N; j++) {
			// Вектор от фиктивного источника к граничной точке
			Vector3d r_vec = bnd_points.row(i).transpose() - fict_points.row(j).transpose();
			double R = r_vec.norm();

			if (R < 1e-12) {
				// Теоретически не должно происходить, если источники вне области
				A(i, j) = 0.0;
			}
			else {
				// Нормальная производная фундаментального решения: ?G/?n = -(n·r)/(4?R?)
				A(i, j) = -ni.dot(r_vec) / (pi4 * R * R * R);
			}
		}

		// Правая часть: заданное значение нормальной производной (Bn)
		b(i) = bnd_Bn(i);
	}

	cout << "END Zapolnyaem " << endl;

	// Условие для устранения неоднозначности (сумма амплитуд = 0)
	A.row(M).setOnes();  // Последняя строка: все единицы
	b(M) = 0.0;          // Сумма q_j = 0

	// Решение системы A * q = b
	VectorXd q;

	cout << "solve" << endl;
	if (use_svd) 
	{
		// SVD - наиболее устойчивый метод для плохо обусловленных матриц
		JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);

		// Устанавливаем порог для сингулярных значений (можно регулировать)
		double threshold = 1e-8 * svd.singularValues()(0);
		q = svd.solve(b);
	}
	else {
		// QR-разложение - быстрее, но может быть менее устойчивым
		HouseholderQR<MatrixXd> qr(A);
		q = qr.solve(b);
	}
	cout << "end solve" << endl;

	return q;
}

/**
 * @brief Вычисляет потенциал и поле в произвольной точке по найденным амплитудам
 */
void computeField(const VectorXd& q,
	const MatrixXd& fict_points,
	const Vector3d& point,
	double& psi,      // Потенциал в точке
	Vector3d& B) {    // Магнитное поле в точке

	psi = 0.0;
	B.setZero();
	const double pi4 = 4.0 * M_PI;

	for (int j = 0; j < q.size(); j++) 
	{
		Vector3d r_vec = point - fict_points.row(j).transpose();
		double R = r_vec.norm();

		if (R > 1e-12) {
			// Потенциал: G = 1/(4?R)
			psi += q(j) / (pi4 * R);

			// Поле: ?G = -r/(4?R?)
			B += -q(j) * r_vec / (pi4 * R * R * R);
		}
	}
}



void Setka::Algoritm(short int alg, Setka* Smain)
{
	// 1  - Плазма МГД
	// 2  - Монте-Карло (для основной сетки) - старый алгорим, теперь используется № 10
	// 3  - Вычисление f_pui по посчитанным S+ S- 
	// 4  - Вычисление n_pui  и  T_pui  по рассчитанным f_pui
	// 5  - Добавить в ячейки основной сетки значение моментов водорода из Монте-Карло (которые посчитаны для сетки MK)
	// 6  - Вычисление функции h0 для розыгрыша пикапов (она считается один раз для каждого сечения перезарядки)  (СТАРАЯ реализация - надо адаптировать)
	// 7  - Вычисление всех интеграллов в ячейках для розыгрыша пикапов (частота и т.д.) 
	// 8  - Вычисление поглощения вдоль заданных лучей (новая реализация через вспомогательную сетку)
	// 9  - (не работает) Перемасштабирование функций распредления водорода (речь про число ячеек AMR), без потери значений (СТАРАЯ реализация - надо адаптировать)
	// 10 - Монте-Карло (новая реализация через вспомогательную сетку)
	// 101 - Монте-Карло (новая реализация через вспомогательную сетку) - полностью имитационный метод без AMR накапливания
	// 11 - расчёт поверхностных токов на разрывах
	// 12 - расчёт объёмных токов
	// 13 - просмотр источников S+/S- и сравнение их с флюидными источниками
	// 14 - расчёт гипотетических токов в сверхзвуковом ветре от HCS
	// 15 - расчёт геометрии HCS 
	// 16 - расчёт потенциального поля во внутреннем слое и различных энергий
	// 17 - расчёт гипотетических токов в гелиошизе от HCS
	// 18 - расчёт геометрии HCS в гелиошизе (не работает, ничего не видно)
	// 19 - расчёт потенциального поля во внутреннем слое методом контрольных объёмов
	// 20 - расчёт потенциального поля в сверхзвуке методом контрольных объёмов
	// 21 - расчёт потенциального поля во внешнем ударном слое методом контрольных объёмов
	// 22 - расчёт потенциального поля в сверхзвуке методом контрольных объёмов - второй порядок
	// 23 - печатаем мини-интерполяционную сетку и источники Sp Sm для Игоря
	// 24 - печатаем карты в Линии H-alpha
	// 25 - Подготовка расчёта инфракрасных спектров (вычисляет плотность и температуру пыли простым методом)
	// 26 - Расчёт инфракрасных спектров
	// 27 - Расчёт температуры пыли методом Монте-Карло

	cout << "Start Algoritm: " << alg << endl;

	this->Test_geometr();
	this->Calculating_measure(0);
	this->Calculating_measure(1);

	if (alg == 1)
	{
		this->Find_Yzel_Sosed_for_BS();

		this->Smooth_angle_HP();
		this->Smooth_head_HP3();
		this->Smooth_head_TS3();

		//this->Go(true, 1000, 1);

		for (int i = 1; i <= 8 * 1; i++) // 6 * 2   12 * 5
		{
			auto start = std::chrono::high_resolution_clock::now();
			cout << "IIIII = " << i << endl;

			//S1.Go(true, 600, 1); // 400   1
			cout << "All time = " << this->phys_param->ALL_Time << endl;
			cout << "All time (in days) = " << this->phys_param->ALL_Time / 0.00142358 << endl;
			cout << "All time (in years) = " << this->phys_param->ALL_Time / 0.519607 << endl;
			this->Go(false, 400, 1); // 400   1
			if (i % 3000000000 == 0)
			{
				//this->Go(true, 100, 1); // 400   1 
			}
			else
			{
				//this->Go(true, 100, 1); // 400   1 
			}
			this->Smooth_head_HP3();
			this->Smooth_head_TS3();

			//S1.Print_parameters_in_some_point();

			this->Tecplot_print_cell_plane_parameters();
			this->Tecplot_print_all_lush_in_2D();

			this->Tecplot_print_all_gran_in_surface("TS");
			this->Tecplot_print_all_gran_in_surface("HP");
			this->Tecplot_print_all_gran_in_surface("BS");

			// Печать результатов
			if (false)
			{
				this->Save_for_interpolate("For_intertpolate_0059-.bin", false);
				Interpol SS = Interpol("For_intertpolate_0059-.bin");

				this->Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
					Eigen::Vector3d(1.0, 0.0, 0.0), "_(1, 0, 0)_" + to_string(this->phys_param->ALL_Time) + "_", 500.0);

				this->Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
					Eigen::Vector3d(cos(const_pi / 18.0), sin(const_pi / 18.0), 0.0), "_(10 deg, 0)_" + to_string(this->phys_param->ALL_Time) + "_", 500.0);

				this->Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
					Eigen::Vector3d(-1.0, 0.0, 0.0), "_(-1, 0, 0)_" + to_string(this->phys_param->ALL_Time) + "_", 500.0);

				this->Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
					Eigen::Vector3d(0.0, 1.0, 0.0), "_(0, 1, 0)_" + to_string(this->phys_param->ALL_Time) + "_", 500.0);

				this->Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_2d_(0, 0, 1, 0)_" + to_string(this->phys_param->ALL_Time) + "_");
			}

			//this->Go(true, 100, 1);
			//this->Tecplot_print_cell_plane_parameters();

			//this->Init_physics();

			if (i % 9 == 0)
			{
				string namn = "parameters_promeg_11" + to_string(i) + ".bin";
				this->Save_cell_parameters(namn);
			}

			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

			std::cout << "Execution time: " << duration.count() / 1000.0 / 60.0 << " minutes" << std::endl;
		}
	}
	else if (alg == 2)
	{
		// Определим зоны для МК
		this->Set_MK_Zone();

		//Проверим зоны
		if (false)
		{
			this->Tecplot_print_gran_with_condition(0);
			this->Tecplot_print_gran_with_condition(1);
			this->Tecplot_print_gran_with_condition(2);
			this->Tecplot_print_gran_with_condition(3);
			this->Tecplot_print_gran_with_condition(4);
			this->Tecplot_print_gran_with_condition(5);
			this->Tecplot_print_gran_with_condition(6);
		}

		// Удаляем какие-то функции распределения
		if (false)
		{
			for (auto& gr : this->All_Gran)
			{
				for (int ii = 0; ii <= 1; ii++)
				{
					string name_f = this->phys_param->AMR_folder + "/" + "func_grans_AMR_" + to_string(ii) + "_H" +
						to_string(1) + "_" + to_string(gr->number) + ".bin";
					if (std::filesystem::exists(name_f))
					{
						std::filesystem::remove(name_f);
					}
				}
			}
		}



		vector<short int> zones_number;
		vector<double> zones_n_koeff;        // Можно для каждой зоны настроить своё количество частиц

		cout << "Start zones_number push_back" << endl;

		zones_number.push_back(6); zones_n_koeff.push_back(1.0);
		zones_number.push_back(4); zones_n_koeff.push_back(1.0);
		zones_number.push_back(2); zones_n_koeff.push_back(1.0);
		zones_number.push_back(1); zones_n_koeff.push_back(1.0);
		zones_number.push_back(3); zones_n_koeff.push_back(1.0);
		zones_number.push_back(5); zones_n_koeff.push_back(1.0);
		zones_number.push_back(7); zones_n_koeff.push_back(1.0);


		short int ijij = 0;
		for (const auto& zone_play : zones_number)
		{
			cout << "Start zone = " << zone_play << endl;
			this->MK_prepare(zone_play);
			this->MK_go(zone_play, int(this->phys_param->N_per_gran * zones_n_koeff[ijij]), nullptr, Smain);
			this->MK_delete(zone_play);
			ijij++;
		}
	}
	else if (alg == 3)
	{
		bool interpol_SS = false;  // Во время вычисления f_pui надо ли интерполировать S+ S- в каждой точке? Или брать среднее в ячейке (это быстрее)

		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		cout << "Create Setka Smc" << endl;
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		Setka Smc = Setka("SDK_40_2D_Setka.bin", "SDK_40_krug_setka.bin", 40);

		cout << "Create SI_main" << endl;
		// Из основной сетки создаём интерполяционную сетку
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		cout << "Move Setka Smc" << endl;
		// Двигаем поверхности вспомогательной сетки к поверхностям основной
		Smc.Move_to_surf(&SI_main);
		// Точно задаём положение внутренней границы сетки
		Smc.geo->R0 = Smc.phys_param->R_0;

		// Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
		Smc.auto_set_luch_geo_parameter(0, true);
		// Настраиваем новую сетку (также как и основную)   [обязательно]
		if (true)
		{
			// Считаем объёмы, площади и другие геометрические характеристики
			Smc.Calculating_measure(0);
			Smc.Calculating_measure(1);

			// Задаём граничные грани
			Smc.Init_boundary_grans();
		}

		// В сетке для MK очистим ненужные имена переменных 
		if (true)
		{
			Smc.phys_param->param_names.assign(Smc.phys_param->MK_param.begin(), Smc.phys_param->MK_param.end());
		}

		// Заполним сетку МК значениями плазмы из основной сетки (чтобы вместо интерполяции в МК использовать значения в центрах ячеек - так быстрее)
		// переинтерполяция
		if (true)
		{
			Smc.PereInterpolate(&SI_main, false);
		}

		Smc.Test_geometr();


		// Загружаем S+ S- для всей сетки
		for (auto& A : Smc.All_Cell)
		{
			A->Init_S(2, Smc.phys_param->pui_nW);
			A->read_S_FromFile(Smc.phys_param->par_n_H_LISM);
		}

		Smc.Save_for_interpolate("For_intertpolate_work_MK.bin", false);
		Interpol SI_MK = Interpol("For_intertpolate_work_MK.bin");


		// Интерполируем S+ S- с малой сетки на большую
		if (interpol_SS == false)
		{
			vector<double> mas_Sm_(this->phys_param->pui_nW);
			vector<double> mas_Sp1_(this->phys_param->pui_nW);
			vector<double> mas_Sp2_(this->phys_param->pui_nW);

			Cell_handle prev_cell_ = Cell_handle();
			Cell_handle next_cell_ = nullptr;

			for (auto& A : this->All_Cell)
			{
				std::fill(mas_Sm_.begin(), mas_Sm_.end(), 0.0);
				std::fill(mas_Sp1_.begin(), mas_Sp1_.end(), 0.0);
				std::fill(mas_Sp2_.begin(), mas_Sp2_.end(), 0.0);

				short int zone = this->determ_zone(A, 0);
				short int kk = 1;;
				A->Init_S(2, this->phys_param->pui_nW);
				if (zone == 2) kk = 2;

				this->Get_pui_SS(mas_Sm_, mas_Sp1_, mas_Sp2_, kk,
					A->center[0][0], A->center[0][1], A->center[0][2],
					Smc, SI_MK, prev_cell_, next_cell_);

				for (size_t i = 0; i < this->phys_param->pui_nW; ++i)
				{
					A->pui_Sm[i] = mas_Sm_[i];
					A->pui_Sp(0, i) = mas_Sp1_[i];
					A->pui_Sp(1, i) = mas_Sp2_[i];
				}
			}


			// Удаляем все S+ S- (чистим память) на малой сетке
			for (auto& A : Smc.All_Cell)
			{
				A->pui_Sm.resize(0);
				A->pui_Sp.resize(0, 0);
			}
		}



		// Считаем функции распределения
		unsigned int st = 0;
		cout << "Start: Culc PUI" << endl;


		// Удаляем некоторые файлы pui (если что-то не так посчиталось
		if (false)
		{
			for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
			{
				auto A = this->All_Cell[idx];
				short int zone = determ_zone(A, 0);
				if (zone == 2)
				{
					std::string filename = "data_pui/func_cells_pui_" + to_string(A->number) + ".bin";
					fs::remove(filename);
				}
			}
		}


		// Важно - источники S+ S- хранятся на сетке МК
		// Но PUI считаются на основной сетке!

#pragma omp parallel for schedule(dynamic)
		for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
		{
			auto A = this->All_Cell[idx];
#pragma omp critical (gergergerg4) 
			{
				st++;
				if (st % 500 == 0)
				{
					cout << "st = " << st << "   from " << this->All_Cell.size() << endl;
				}
			}

			//std::string filename = "data_pui/func_cells_pui_" + to_string(A->number) + ".bin";
			//if (file_exists(filename) == true) continue;

			short int zone = this->determ_zone(A, 0);
			//cout << "A" << endl;
			A->Init_f_pui(this->phys_param->pui_nW, zone);
			//cout << "B" << endl;
			this->Culc_f_pui_in_cell(A, Smc, SI_main, SI_MK, interpol_SS);
			//cout << "C" << endl;
			A->write_pui_ToFile();
			//cout << "D" << endl;
			A->Delete_f_pui();
			//cout << "F" << endl;
		}
		cout << "End: Culc PUI" << endl;


		// Удаляем все S+ S- (чистим память)
		for (auto& A : Smc.All_Cell)
		{
			A->pui_Sm.resize(0);
			A->pui_Sp.resize(0, 0);
		}

		for (auto& A : this->All_Cell)
		{
			A->pui_Sm.resize(0);
			A->pui_Sp.resize(0, 0);
		}


		this->Print_pui(17.0, 0.0, 0.0);
		this->Print_pui(20.0, 0.0, 0.0);
		this->Print_pui(25.0, 0.0, 0.0);
		this->Print_pui(1.0, 0.0, 0.0);
		this->Print_pui(5.0, 0.0, 0.0);
		this->Print_pui(10.0, 0.0, 0.0);
		this->Print_pui(15.0, 0.0, 0.0);
		this->Print_pui(28.0, 0.0, 0.0);
		this->Print_pui(50.0, 0.0, 0.0);
		this->Print_pui(100.0, 0.0, 0.0);
		this->Print_pui(200.0, 0.0, 0.0);
	}
	else if (alg == 4)
	{
		unsigned int st = 0;
#pragma omp parallel for schedule(dynamic)
		for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
		{
			auto A = this->All_Cell[idx];
#pragma omp critical (gergergerg4) 
			{
				st++;
				if (st % 10000 == 0)
				{
					cout << "st = " << st << "   from " << this->All_Cell.size() << endl;
				}
			}

			std::string filename = "data_pui/func_cells_pui_" + to_string(A->number) + ".bin";
			if (file_exists(filename) != true)
			{
				cout << "Error  Net pui! fergebhr6yveybhe5vte " << endl;
				exit(-1);
			}

			short int zone = determ_zone(A, 0);
			A->Init_f_pui(this->phys_param->pui_nW, zone);
			A->read_pui_FromFile();
			// Здесь вычисляем нужные моменты
			A->culc_pui_n_T(this->phys_param->pui_wR);
			A->Delete_f_pui();
		}

		this->phys_param->param_names.push_back("MK_rho_Pui_1");
		this->phys_param->param_names.push_back("MK_T_Pui_1");
		this->phys_param->param_names.push_back("MK_rho_Pui_2");
		this->phys_param->param_names.push_back("MK_T_Pui_2");

	}
	else if (alg == 5)
	{
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		cout << "Create Setka Smc" << endl;
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		Setka Smc = Setka("SDK_40_2D_Setka.bin", "SDK_40_krug_setka.bin", 40);

		cout << "Create SI_main" << endl;
		// Из основной сетки создаём интерполяционную сетку
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		cout << "Move Setka Smc" << endl;
		// Двигаем поверхности вспомогательной сетки к поверхностям основной
		Smc.Move_to_surf(&SI_main);
		// Точно задаём положение внутренней границы сетки
		Smc.geo->R0 = Smc.phys_param->R_0;

		// Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
		Smc.auto_set_luch_geo_parameter(0, true);

		// Настраиваем новую сетку (также как и основную)   [обязательно]
		if (true)
		{
			// Считаем объёмы, площади и другие геометрические характеристики
			Smc.Calculating_measure(0);
			Smc.Calculating_measure(1);

			// Задаём граничные грани
			Smc.Init_boundary_grans();
		}

		// В сетке для MK очистим ненужные имена переменных 
		if (true)
		{
			Smc.phys_param->param_names.assign(Smc.phys_param->MK_param.begin(), Smc.phys_param->MK_param.end());
		}

		Smc.Test_geometr();

		// Так как в файле хранятся параметры во всех ячейках (даже в тех, которые были за пределом рассчитанной зоны)
		// Нужно сначала скачать все кроме текущей зоны (так как они только что посчитаны), а потом записать все

		if (file_exists(Smc.phys_param->MK_file))
		{
			Smc.Download_cell_MK_parameters(Smc.phys_param->MK_file, -1);
		}
		else
		{
			cout << "Error 94ut9yegfh9perfg8yvowjrgf9348" << endl;
		}

		cout << "Create SI_MK" << endl;
		// Из MK сетки создаём интерполяционную сетку
		Smc.Save_for_interpolate("For_intertpolate_work_MK.bin", false);

		// Переинтерполируем параметры Монте-Карло из вспомогательной сетки в основную
		this->PereInterpolate("For_intertpolate_work_MK.bin", false, true);
	}
	else if (alg == 6)
	{
		this->Culc_h0_for_pui(); // Считаеи h0 и сразу записывает в файл
	}
	else if (alg == 7)
	{
		unsigned int st = 0;
		#pragma omp parallel for schedule(dynamic)
		for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
		{
			#pragma omp critical (gergergerg4) 
			{
				st++;
				if (st % 5000 == 0)
				{
					cout << "st = " << st << "   from " << this->All_Cell.size() << endl;
				}
			}

		auto A = this->All_Cell[idx];
		short int zone = determ_zone(A, 0);

		/*if (A->number != 378182)
		{
			continue;
		}
		else
		{
			cout << "Culc  378182 " << endl;
		}*/

		std::string filename = "data_pui_intergal/func_cells_pui_integral_" + to_string(A->number) + ".bin";
		//if (file_exists(filename) == true) continue;

		A->Init_f_pui(this->phys_param->pui_nW, zone);
		A->read_pui_FromFile();

		A->culc_pui_n_T(this->phys_param->pui_wR);
		A->Init_pui_integral(this->phys_param->pui_F_n, zone);
		A->pui_integral_Culc(this->phys_param);
		A->write_pui_integral_ToFile();
		A->Delete_pui_integral();
		A->Delete_f_pui();
		}

		this->phys_param->param_names.push_back("MK_rho_Pui_1");
		this->phys_param->param_names.push_back("MK_T_Pui_1");
		this->phys_param->param_names.push_back("MK_rho_Pui_2");
		this->phys_param->param_names.push_back("MK_T_Pui_2");
	}
	else if (alg == 8)
	{
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		//Setka* Smc;
		//this->Create_mini_Setka_for_MK(Smc);
		cout << "Create Setka Smc" << endl;
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		Setka Smc = Setka("SDK_40_2D_Setka.bin", "SDK_40_krug_setka.bin", 40);

		cout << "Create SI_main" << endl;
		// Из основной сетки создаём интерполяционную сетку
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		cout << "Move Setka Smc" << endl;
		// Двигаем поверхности вспомогательной сетки к поверхностям основной
		Smc.Move_to_surf(&SI_main);
		// Точно задаём положение внутренней границы сетки
		Smc.geo->R0 = Smc.phys_param->R_0;

		// Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
		Smc.auto_set_luch_geo_parameter(0, true);
		// Настраиваем новую сетку (также как и основную)   [обязательно]
		if (true)
		{
			// Считаем объёмы, площади и другие геометрические характеристики
			Smc.Calculating_measure(0);
			Smc.Calculating_measure(1);

			// Задаём граничные грани
			Smc.Init_boundary_grans();
		}


		// В сетке для MK очистим ненужные имена переменных 
		if (false)
		{
			Smc.phys_param->param_names.assign(Smc.phys_param->MK_param.begin(), Smc.phys_param->MK_param.end());
		}

		// Заполним сетку МК значениями плазмы из основной сетки (чтобы вместо интерполяции в МК использовать значения в центрах ячеек - так быстрее)
		// переинтерполяция
		if (true)
		{
			Smc.PereInterpolate(&SI_main, false, false);
		}

		Smc.Test_geometr();

		cout << "Reading arrays from files" << Smc.phys_param->pogl_folder << "  " << Smc.phys_param->pogl_n << endl;
		cout << static_cast<int>(Smc.phys_param->num_H) << endl;
		cout << Smc.phys_param->R_MK_Max << endl;

		for (size_t idx = 0; idx < Smc.All_Cell.size(); ++idx)
		{
			//cout << "A" << endl;
			auto A = Smc.All_Cell[idx];
			A->Init_mas_pogl(Smc.phys_param->pogl_n, Smc.phys_param->num_H);
			//cout << "B" << endl;
			A->read_mas_pogl_FromFile(Smc.phys_param);
			//cout << "C" << endl;
		}

		Cell* prev = nullptr;
		Cell* A = Smc.Find_cell_point(40.0, 0.0, 0.0, 0, prev);
		cout << "1 Sum = " << A->mas_pogl.sum() << endl;

		// Считываем моменты
		cout << "Reading moments from " << Smc.phys_param->MK_file << endl;
		if (true)
		{
			if (file_exists(Smc.phys_param->MK_file))
			{
				Smc.Download_cell_MK_parameters(Smc.phys_param->MK_file, 1000);
			}
		}

		Smc.Set_MK_Zone();

		cout << "2 Sum = " << A->mas_pogl.sum() << endl;

		cout << "Arrays read successfully" << endl;

		Smc.mas_pogl_Culc(1.0, 0.0, 0.0, "upwind");
		cout << "start mas_pogl_Culc_fluid: " << endl;
		this->mas_pogl_Culc_fluid(1.0, 0.0, 0.0, "upwind");
		cout << "end mas_pogl_Culc_fluid: " << endl;


		Smc.mas_pogl_Culc(0.985132, -0.169234, 0.0295701, "36Oph");
		cout << "start mas_pogl_Culc_fluid: " << endl;
		this->mas_pogl_Culc_fluid(0.985132, -0.169234, 0.0295701, "36Oph");
		cout << "end mas_pogl_Culc_fluid: " << endl;

		Smc.mas_pogl_Culc(1.0, 0.1, 0.0, "sim_upwind");
		this->mas_pogl_Culc_fluid(1.0, 0.1, 0.0, "sim_upwind");
		Smc.mas_pogl_Culc(0.0, 1.0, 0.0, "crosswind1");
		this->mas_pogl_Culc_fluid(0.0, 1.0, 0.0, "crosswind1");
		Smc.mas_pogl_Culc(0.0, 1.0, 1.0, "crosswind2");
		this->mas_pogl_Culc_fluid(0.0, 1.0, 1.0, "crosswind2");
		Smc.mas_pogl_Culc(0.0, 0.0, 1.0, "crosswind3");
		this->mas_pogl_Culc_fluid(0.0, 0.0, 1.0, "crosswind3");
		Smc.mas_pogl_Culc(-1.0, 0.0, 0.0, "downwind");
		this->mas_pogl_Culc_fluid(-1.0, 0.0, 0.0, "downwind");
		Smc.mas_pogl_Culc(-1.0, 1.0, 0.0, "tail1");
		this->mas_pogl_Culc_fluid(-1.0, 1.0, 0.0, "tail1");
		Smc.mas_pogl_Culc(-1.0, 0.70710678, 0.70710678, "tail2");
		this->mas_pogl_Culc_fluid(-1.0, 0.70710678, 0.70710678, "tail2");
		Smc.mas_pogl_Culc(-1.0, 0.0, 1.0, "tail3");
		this->mas_pogl_Culc_fluid(-1.0, 0.0, 1.0, "tail3");

		cout << "Removing arrays" << endl;

		for (size_t idx = 0; idx < Smc.All_Cell.size(); ++idx)
		{
			auto A = Smc.All_Cell[idx];
			A->Delete_mas_pogl();
		}

		Smc.Print_f_proect_in_gran(1);
		Smc.Print_f_proect_in_gran(2);
		Smc.Print_f_proect_in_gran(3);

		Smc.Print_f_proect_in_cell(13.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(20.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(35.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(40.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(45.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(50.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(55.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(60.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(70.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(80.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(90.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(-40.0, 0.0, 0.0);
		Smc.Print_f_proect_in_cell(-80.0, 0.0, 0.0);

		Smc.Print_fH(2, Type_Gran_surf::TS, 1.0, 0.0, 0.0, 7.0 * const_pi / 180.0);
		Smc.Print_fH(2, Type_Gran_surf::HP, 1.0, 0.0, 0.0, 7.0 * const_pi / 180.0);
	}
	else if (alg == 9)
	{
		short int sortH = 5; // Какой сорт водорода будет менять?  1-4
		double Diapazon = 100.0; // Какой новый диапазон функции

		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		cout << "Create Setka Smc" << endl;
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		Setka Smc = Setka("SDK_40_2D_Setka.bin", "SDK_40_krug_setka.bin", 40);

		cout << "Create SI_main" << endl;
		// Из основной сетки создаём интерполяционную сетку
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		cout << "Move Setka Smc" << endl;
		// Двигаем поверхности вспомогательной сетки к поверхностям основной
		Smc.Move_to_surf(&SI_main);
		// Точно задаём положение внутренней границы сетки
		Smc.geo->R0 = Smc.phys_param->R_0;

		// Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
		Smc.auto_set_luch_geo_parameter(0, true);
		// Настраиваем новую сетку (также как и основную)   [обязательно]
		if (true)
		{
			// Считаем объёмы, площади и другие геометрические характеристики
			Smc.Calculating_measure(0);
			Smc.Calculating_measure(1);

			// Задаём граничные грани
			Smc.Init_boundary_grans();

			// Проверки
			if (this->phys_param->is_PUI != Smc.phys_param->is_PUI)
			{
				cout << "Error eijrgfouiehg384tfg7gf" << endl;
				exit(-1);
			}
		}


		Smc.Test_geometr();

		if (true)
		{
			unsigned int in = 0;

			#pragma omp parallel for schedule(dynamic)
			for (size_t idx = 0; idx < Smc.All_Gran.size(); ++idx)
			{
				auto gr = Smc.All_Gran[idx];
				#pragma omp critical (gergergerg4) 
				{
					in++;
					if (in % 10000 == 0)
					{
						cout << "Gran: " << in << "  /  " << Smc.All_Gran.size() << endl;
					}
				}
				for (int ii = 0; ii <= 1; ii++)
				{
					string name_f = this->phys_param->AMR_folder + "/" + "func_grans_AMR_" + to_string(ii) + "_H" +
						to_string(sortH) + "_" + to_string(gr->number) + ".bin";
					//cout << "A1" << endl;
					if (std::filesystem::exists(name_f))
					{
						// Выделяем место под AMR, сколько сортов водорода, столько и места
						if (gr->AMR.size() < this->phys_param->num_H)
						{
							//cout << "A2" << endl;
							gr->AMR.resize(this->phys_param->num_H);
							for (size_t i = 0; i < this->phys_param->num_H; i++)
							{
								gr->AMR[i][0] = nullptr;
								gr->AMR[i][1] = nullptr;
							}
						}
						//cout << "A3" << endl;

						// Считываем AMR
						gr->Read_AMR(ii, sortH, this->phys_param, false);
						auto func = gr->AMR[sortH - 1][ii];
						//cout << "A4" << endl;
						std::vector<AMR_cell*> cells_amr;
						std::vector<std::pair<Point, size_t>> points; // точки и их номера для построения триангуляции
						std::vector <Int_point*> ALL_Cells;     // Точки в которых хранятся параметры
						std::array<double, 3> center;
						Delaunay* Delone;
						double Vx, Vy, Vz;
						unsigned int i = 0;
						//cout << "A5" << endl;
						func->Get_all_cells(cells_amr);
						for (const auto& cell : cells_amr)
						{
							cell->Get_Center(func, center);
							//func->Get_real_koordinate(center[0], center[1], center[2], Vx, Vy, Vz);
							auto A = new Int_point(center[0], center[1], center[2]);
							A->parameters["f"] = cell->getF();
							points.push_back({ {center[0], center[1], center[2]}, i });
							ALL_Cells.push_back(A);
							i++;
						}
						//cout << "A6" << endl;
						Delone = new Delaunay(points.begin(), points.end());

						auto new_func = new AMR_f();
						new_func->AMR_self = new_func;
						//cout << "A7" << endl;
						new_func->AMR_resize(0.0, Diapazon, -Diapazon, Diapazon,                      // ЗДЕСЬ НАПИСАН ДИАПОЗОН ИЗМЕНЕНИЯ
							-Diapazon, Diapazon, 3, 6, 6);
						//cout << "A71" << endl;
						new_func->Copy_and_Refine(ALL_Cells, Delone);
						//cout << "A8" << endl;
						//cout << "Copy_and_Refine:  " << func->Size() << "   " << new_func->Size() << endl;

						func->Delete();
						gr->AMR[sortH - 1][ii] = new_func;
						//cout << "A9" << endl;
						delete Delone;
						for (auto& i : ALL_Cells)
						{
							delete i;
						}
						//cout << "A10" << endl;
						ALL_Cells.clear();
						std::filesystem::remove(name_f);
						new_func->Save(name_f);

						gr->AMR.clear();
					}
				}
			}
		}
	}
	else if (alg == 10)
	{
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		cout << "Create Setka Smc" << endl;
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		Setka Smc = Setka("SDK_40_2D_Setka.bin", "SDK_40_krug_setka.bin", 40);

		cout << "Create SI_main" << endl;
		// Из основной сетки создаём интерполяционную сетку
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		cout << "Move Setka Smc" << endl;
		// Двигаем поверхности вспомогательной сетки к поверхностям основной
		Smc.Move_to_surf(&SI_main);
		// Точно задаём положение внутренней границы сетки
		Smc.geo->R0 = Smc.phys_param->R_0;

		// Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
		Smc.auto_set_luch_geo_parameter(0, true);
		// Настраиваем новую сетку (также как и основную)   [обязательно]
		if (true)
		{
			// Считаем объёмы, площади и другие геометрические характеристики
			Smc.Calculating_measure(0);
			Smc.Calculating_measure(1);

			// Задаём граничные грани
			Smc.Init_boundary_grans();

			// Проверки
			if (this->phys_param->is_PUI != Smc.phys_param->is_PUI)
			{
				cout << "Error eijrgfouiehg384tfg7gf" << endl;
				exit(-1);
			}
		}

		// Визуализация новой сетки для проверки   [опционально]
		if (true)
		{
			Smc.Tecplot_print_all_lush_in_2D();
			Smc.Tecplot_print_2D_setka(0.0, 0.0, 1.0, -0.00001, "Smc_setka_2d_(0, 0, 1, 0)_");
			Smc.Tecplot_print_2D_setka(0.0, 1.0, 0.0, -0.00001, "Smc_setka_2d_(0, 1, 0, 0)_");
			Smc.Tecplot_print_2D_setka(0.0, 1.0, 1.0, -0.00001, "Smc_setka_2d_(0, 1, 1, 0)_");
			Smc.Tecplot_print_all_gran_in_surface("TS");
			Smc.Tecplot_print_all_gran_in_surface("HP");
			Smc.Tecplot_print_all_gran_in_surface("BS");
		}

		// В сетке для MK очистим ненужные имена переменных 
		if (true)
		{
			Smc.phys_param->param_names.assign(Smc.phys_param->MK_param.begin(), Smc.phys_param->MK_param.end());
		}

		// Заполним сетку МК значениями плазмы из основной сетки (чтобы вместо интерполяции в МК использовать значения в центрах ячеек - так быстрее)
		// переинтерполяция
		if (true)
		{
			Smc.PereInterpolate(&SI_main, false);
		}

		Smc.Test_geometr();

		// Настройка всех массивов для расчёта пикапов
		if (Smc.phys_param->is_PUI == true)
		{
			cout << "Download PUI" << endl;
			// Загружаем h0
			Smc.Init_h0_and_read_from_file();
			this->Init_h0_and_read_from_file();

			// Загружаем все интеграллы пикапов
			unsigned int st = 0;
			#pragma omp parallel for schedule(dynamic)
			for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
			{
				#pragma omp critical (gergergerg4) 
				{
					st++;
					if (st % 50000 == 0)
					{
						cout << "step = " << st << "   from " << this->All_Cell.size() << endl;
					}
				}

				auto A = this->All_Cell[idx];
				short int zone = determ_zone(A, 0);
				A->Init_pui_integral(this->phys_param->pui_F_n, zone);
				A->read_pui_integral_FromFile(this->phys_param);
				A->Init_f_pui(this->phys_param->pui_nW, zone);
				A->read_pui_FromFile();
				if (false)//(idx == 2200)
				{
					//A->pui_integral_Culc(this->phys_param);
					A->print_pui(this->phys_param->pui_wR, "2200_pui");

					cout << "FF = " << A->pui_get_f(40.0, 0, phys_param->pui_wR) << " " << 
						A->pui_get_f(10.0, 0, phys_param->pui_wR) << " " << 
						A->pui_get_f(0.0, 0, phys_param->pui_wR) << " " << 
						A->pui_get_f(-10.0, 0, phys_param->pui_wR) << " " << 
						A->pui_get_f(199.8, 0, phys_param->pui_wR) << " " << 
						A->pui_get_f(160.0, 0, phys_param->pui_wR) << " " << endl;
				}
				A->culc_pui_n_T(this->phys_param->pui_wR);
				A->Delete_f_pui();
			}
		}

		// Проверим, загрузились ли массивы
		if (true)
		{
			cout << "Proverka chastot pui" << endl;
			Cell* CC;
			Cell* prev = nullptr;
			CC = this->Find_cell_point(20.0, 0.0, 0.0, 0, prev);
			//CC = this->All_Cell[2200];

			CC->print_nu_integr_pui(this->phys_param);
			CC->print_F_integr_pui();

			double nu = CC->pui_get_nu(5.0, 0, this->phys_param->pui_wR);
			if (nu <= 0.0)
			{
				cout << "Warning iudrhguseroigfsegsr" << endl;
				cout << CC->center[0][0] << " " << CC->center[0][1] << " " << CC->center[0][2] << endl;
				cout << int(CC->type) << endl;
			}
			//cout << "nu1 = " << nu << endl;
			nu = CC->pui_get_nu(5.0, 1, this->phys_param->pui_wR);
			if (nu <= 0.0)
			{
				cout << "Warning ghietgy87tg98e9g98e" << endl;
				cout << CC->center[0][0] << " " << CC->center[0][1] << " " << CC->center[0][2] << endl;
				cout << int(CC->type) << endl;
			}
			//cout << "nu2 = " << nu << endl;

			/*cout << CC->pui_get_nu(50, 0, this->phys_param->pui_wR) << " " <<
				CC->pui_get_nu(50, 1, this->phys_param->pui_wR) << " " <<
				CC->pui_get_nu(0.0, 0, this->phys_param->pui_wR) << " " <<
				CC->pui_get_nu(-10.0, 0, this->phys_param->pui_wR) << " " <<
				CC->pui_get_nu(250, 1, this->phys_param->pui_wR) << " " <<
				CC->pui_get_nu(-50, 1, this->phys_param->pui_wR) << endl;*/
		}

		//return;

		// Удаляем какие-то функции распределения
		if (false)
		{
			for (auto& gr : Smc.All_Gran)
			{
				for (int ii = 0; ii <= 1; ii++)
				{
					string name_f = Smc.phys_param->AMR_folder + "/" + "func_grans_AMR_" + to_string(ii) + "_H" +
						to_string(8) + "_" + to_string(gr->number) + ".bin";

					if (std::filesystem::exists(name_f))
					{
						std::filesystem::remove(name_f);
					}
				}
			}
		}

		cout << "Set MK zone" << endl;
		// Определим зоны для МК
		Smc.Set_MK_Zone();

		//Проверим зоны   [опционально]
		if (true)
		{
			Smc.Tecplot_print_gran_with_condition(0);
			Smc.Tecplot_print_gran_with_condition(1);
			Smc.Tecplot_print_gran_with_condition(2);
			Smc.Tecplot_print_gran_with_condition(3);
			Smc.Tecplot_print_gran_with_condition(4);
			Smc.Tecplot_print_gran_with_condition(5);
			Smc.Tecplot_print_gran_with_condition(6);
		}


		vector<short int> zones_number;
		vector<double> zones_n_koeff;        // Можно для каждой зоны настроить своё количество частиц

		cout << "Start zones_number push_back" << endl;

		//zones_number.push_back(1); zones_n_koeff.push_back(1.0);
		//zones_number.push_back(2); zones_n_koeff.push_back(1.0);

		zones_number.push_back(6); zones_n_koeff.push_back(1.0);
		zones_number.push_back(6); zones_n_koeff.push_back(1.0);
		zones_number.push_back(4); zones_n_koeff.push_back(1.0);
		zones_number.push_back(4); zones_n_koeff.push_back(1.0);
		zones_number.push_back(2); zones_n_koeff.push_back(1.0);
		//zones_number.push_back(2); zones_n_koeff.push_back(1.0);
		//zones_number.push_back(6); zones_n_koeff.push_back(1.0);
		

		short int ijij = 0;
		for (const auto& zone_play : zones_number)
		{
			cout << "Start zone = " << zone_play << endl;
			Smc.MK_prepare(zone_play);
			//Smc.MK_go(zone_play, int(this->phys_param->N_per_gran * zones_n_koeff[ijij]), &SI_main);
			Smc.MK_go(zone_play, int(this->phys_param->N_per_gran * zones_n_koeff[ijij]), nullptr, Smain);
			Smc.MK_delete(zone_play);
			ijij++;
		}

		cout << "Create SI_MK" << endl;
		// Из основной сетки создаём интерполяционную сетку
		Smc.Save_for_interpolate("For_intertpolate_work_MK.bin", false);

		// Переинтерполируем параметры Монте-Карло из вспомогательной сетки в основную
		this->PereInterpolate("For_intertpolate_work_MK.bin", false, true);


		// Очистка
		if (Smc.phys_param->is_PUI == true)
		{
			// Загружаем h0
			Smc.Delete_h0();
			this->Delete_h0();

			// Загружаем все интеграллы пикапов
			for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
			{
				auto A = this->All_Cell[idx];
				short int zone = determ_zone(A, 0);
				A->Delete_pui_integral();
			}
		}

	}
	else if (alg == 101)
	{
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		cout << "Create Setka Smc" << endl;
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		Setka Smc = Setka("SDK_40_2D_Setka.bin", "SDK_40_krug_setka.bin", 40);

		cout << "Create SI_main" << endl;
		// Из основной сетки создаём интерполяционную сетку
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		cout << "Move Setka Smc" << endl;
		// Двигаем поверхности вспомогательной сетки к поверхностям основной
		Smc.Move_to_surf(&SI_main);
		// Точно задаём положение внутренней границы сетки
		Smc.geo->R0 = Smc.phys_param->R_0;

		// Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
		Smc.auto_set_luch_geo_parameter(0, true);
		// Настраиваем новую сетку (также как и основную)   [обязательно]
		if (true)
		{
			// Считаем объёмы, площади и другие геометрические характеристики
			Smc.Calculating_measure(0);
			Smc.Calculating_measure(1);

			// Задаём граничные грани
			Smc.Init_boundary_grans();

			// Проверки
			if (this->phys_param->is_PUI != Smc.phys_param->is_PUI)
			{
				cout << "Error eijrgfouiehg384tfg7gf" << endl;
				exit(-1);
			}
		}

		// Визуализация новой сетки для проверки   [опционально]
		if (true)
		{
			Smc.Tecplot_print_all_lush_in_2D();
			Smc.Tecplot_print_2D_setka(0.0, 0.0, 1.0, -0.00001, "Smc_setka_2d_(0, 0, 1, 0)_");
			Smc.Tecplot_print_2D_setka(0.0, 1.0, 0.0, -0.00001, "Smc_setka_2d_(0, 1, 0, 0)_");
			Smc.Tecplot_print_2D_setka(0.0, 1.0, 1.0, -0.00001, "Smc_setka_2d_(0, 1, 1, 0)_");
			Smc.Tecplot_print_all_gran_in_surface("TS");
			Smc.Tecplot_print_all_gran_in_surface("HP");
			Smc.Tecplot_print_all_gran_in_surface("BS");
		}

		// В сетке для MK очистим ненужные имена переменных 
		if (true)
		{
			Smc.phys_param->param_names.assign(Smc.phys_param->MK_param.begin(), Smc.phys_param->MK_param.end());
		}

		// Заполним сетку МК значениями плазмы из основной сетки (чтобы вместо интерполяции в МК использовать значения в центрах ячеек - так быстрее)
		// переинтерполяция
		if (true)
		{
			Smc.PereInterpolate(&SI_main, false);
		}

		Smc.Test_geometr();

		// Настройка всех массивов для расчёта пикапов
		if (Smc.phys_param->is_PUI == true)
		{
			cout << "Download PUI" << endl;
			// Загружаем h0
			Smc.Init_h0_and_read_from_file();
			this->Init_h0_and_read_from_file();

			// Загружаем все интеграллы пикапов
			unsigned int st = 0;
			#pragma omp parallel for schedule(dynamic)
			for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
			{
				#pragma omp critical (gergergerg4) 
				{
					st++;
					if (st % 50000 == 0)
					{
						cout << "step = " << st << "   from " << this->All_Cell.size() << endl;
					}
				}

			auto A = this->All_Cell[idx];
			short int zone = determ_zone(A, 0);
			A->Init_pui_integral(this->phys_param->pui_F_n, zone);
			A->read_pui_integral_FromFile(this->phys_param);
			A->Init_f_pui(this->phys_param->pui_nW, zone);
			A->read_pui_FromFile();
			if (false)//(idx == 2200)
			{
				//A->pui_integral_Culc(this->phys_param);
				A->print_pui(this->phys_param->pui_wR, "2200_pui");

				cout << "FF = " << A->pui_get_f(40.0, 0, phys_param->pui_wR) << " " <<
					A->pui_get_f(10.0, 0, phys_param->pui_wR) << " " <<
					A->pui_get_f(0.0, 0, phys_param->pui_wR) << " " <<
					A->pui_get_f(-10.0, 0, phys_param->pui_wR) << " " <<
					A->pui_get_f(199.8, 0, phys_param->pui_wR) << " " <<
					A->pui_get_f(160.0, 0, phys_param->pui_wR) << " " << endl;
			}
			A->culc_pui_n_T(this->phys_param->pui_wR);
			A->Delete_f_pui();
			}
		}

		// Проверим, загрузились ли массивы
		if (Smc.phys_param->is_PUI == true)
		{
			cout << "Proverka chastot pui" << endl;
			Cell* CC;
			Cell* prev = nullptr;
			CC = this->Find_cell_point(20.0, 0.0, 0.0, 0, prev);
			//CC = this->All_Cell[2200];

			CC->print_nu_integr_pui(this->phys_param);
			CC->print_F_integr_pui();

			double nu = CC->pui_get_nu(5.0, 0, this->phys_param->pui_wR);
			if (nu <= 0.0)
			{
				cout << "Warning iudrhguseroigfsegsr" << endl;
				cout << CC->center[0][0] << " " << CC->center[0][1] << " " << CC->center[0][2] << endl;
				cout << int(CC->type) << endl;
			}
			//cout << "nu1 = " << nu << endl;
			nu = CC->pui_get_nu(5.0, 1, this->phys_param->pui_wR);
			if (nu <= 0.0)
			{
				cout << "Warning ghietgy87tg98e9g98e" << endl;
				cout << CC->center[0][0] << " " << CC->center[0][1] << " " << CC->center[0][2] << endl;
				cout << int(CC->type) << endl;
			}
		}

		//return;

		// Удаляем какие-то функции распределения
		if (false)
		{
			for (auto& gr : Smc.All_Gran)
			{
				for (int ii = 0; ii <= 1; ii++)
				{
					string name_f = Smc.phys_param->AMR_folder + "/" + "func_grans_AMR_" + to_string(ii) + "_H" +
						to_string(8) + "_" + to_string(gr->number) + ".bin";

					if (std::filesystem::exists(name_f))
					{
						std::filesystem::remove(name_f);
					}
				}
			}
		}

		cout << "Set MK zone" << endl;
		// Определим зоны для МК
		Smc.Set_MK_Zone();

		//Проверим зоны   [опционально]
		if (true)
		{
			Smc.Tecplot_print_gran_with_condition(0);
			Smc.Tecplot_print_gran_with_condition(1);
			Smc.Tecplot_print_gran_with_condition(2);
			Smc.Tecplot_print_gran_with_condition(3);
			Smc.Tecplot_print_gran_with_condition(4);
			Smc.Tecplot_print_gran_with_condition(5);
			Smc.Tecplot_print_gran_with_condition(6);
		}


		vector<short int> zones_number;
		vector<double> zones_n_koeff;        // Можно для каждой зоны настроить своё количество частиц

		cout << "Start zones_number push_back" << endl;


		zones_number.push_back(6); zones_n_koeff.push_back(1.0);

		short int ijij = 0;
		for (const auto& zone_play : zones_number)
		{
			cout << "Start zone = " << zone_play << endl;
			Smc.MK_prepare(zone_play, false);
			//Smc.MK_go(zone_play, int(this->phys_param->N_per_gran * zones_n_koeff[ijij]), &SI_main);
			Smc.MK_go_Imit(zone_play, int(this->phys_param->N_per_gran * zones_n_koeff[ijij]), nullptr, Smain);
			Smc.MK_delete(zone_play, false);
			ijij++;
		}

		// Очистка
		if (Smc.phys_param->is_PUI == true)
		{
			// Загружаем h0
			Smc.Delete_h0();
			this->Delete_h0();

			// Загружаем все интеграллы пикапов
			for (size_t idx = 0; idx < this->All_Cell.size(); ++idx)
			{
				auto A = this->All_Cell[idx];
				short int zone = determ_zone(A, 0);
				A->Delete_pui_integral();
			}
		}

		}
	else if (alg == 11)
	{
		ofstream fout;
		string name_f;

		//HP
		if (true)
		{
			name_f = "HP_J_dissipation_6_year.txt";
			fout.open(name_f);
			fout << "TITLE = HP  VARIABLES = x, y, z, r, phi, the, Jx, Jy, Jz, |J|, J2x, J2y, J2z, |J2|, Bx_L, By_L, Bz_L, Bx_R, By_R, Bz_R" << endl;
			fout << "ZONE T=HP, N = " << this->Gran_HP.size() * 4 << ", E = " << this->Gran_HP.size() << ", F=FEPOINT, ET=quadrilateral" << endl;

			for (const auto& i : this->Gran_HP)
			{
				auto A = i->cells[0];
				auto B = i->cells[1];
				Eigen::Vector3d n, B1, B2, cc, J2, J, BB1, BB2;

				n[0] = i->normal[0][0];
				n[1] = i->normal[0][1];
				n[2] = i->normal[0][2];

				if (A->type == Type_cell::Zone_3)
				{
					A = i->cells[1];
					B = i->cells[0];

					n[0] = -i->normal[0][0];
					n[1] = -i->normal[0][1];
					n[2] = -i->normal[0][2];
				}


				B1[0] = A->parameters[0]["Bx"] / 6.0;
				B1[1] = A->parameters[0]["By"] / 6.0;
				B1[2] = A->parameters[0]["Bz"] / 6.0;

				B2[0] = B->parameters[0]["Bx"];
				B2[1] = B->parameters[0]["By"];
				B2[2] = B->parameters[0]["Bz"];


				BB1 = B2 - B1;
				BB2 = B2 + B1;

				//cout << "do = " << B1[0] << endl;
				J = n.cross(BB1);
				J2 = n.cross(BB2);
				//cout << "posle = " << B1[0] << endl;

				J = J * 3.73834;
				J2 = J2 * 3.73834;

				for (auto& j : i->yzels)
				{
					cc[0] = j->coord[0][0];
					cc[1] = j->coord[0][1];
					cc[2] = j->coord[0][2];

					fout << cc[0] << " " << cc[1] << " " << cc[2] << " " << norm2(cc[0], cc[1], cc[2]) << " " <<
						polar_angle(cc[1], cc[2]) << " " << polar_angle(cc[0], norm2(0.0, cc[1], cc[2])) << " " <<
						J[0] << " " << J[1] << " " << J[2] << " " << J.norm() << " " <<
						J2[0] << " " << J2[1] << " " << J2[2] << " " << J2.norm() << " " <<
						B1[0] << " " << B1[1] << " " << B1[2] << " " <<
						B2[0] << " " << B2[1] << " " << B2[2] << " " << endl;
				}
			}


			for (int k = 0; k < this->Gran_HP.size(); k++)
			{
				fout << 4 * k + 1 << " " << 4 * k + 2 << " " << 4 * k + 3 << " " << 4 * k + 4 << endl;
			}


			fout.close();
		}

		// TS
		if (false)
		{
			name_f = "TS_J_with_polatiry.txt";
			fout.open(name_f);
			fout << "TITLE = HP  VARIABLES = x, y, z, phi, the, Jx, Jy, Jz, |J|" << endl;
			fout << "ZONE T=HP, N = " << this->Gran_TS.size() * 4 << ", E = " << this->Gran_TS.size() << ", F=FEPOINT, ET=quadrilateral" << endl;

			for (const auto& i : this->Gran_TS)
			{
				auto A = i->cells[0];
				auto B = i->cells[1];
				Eigen::Vector3d n, B1, B2, cc;

				n[0] = i->normal[0][0];
				n[1] = i->normal[0][1];
				n[2] = i->normal[0][2];

				B1[0] = A->parameters[0]["Bx"];
				B1[1] = A->parameters[0]["By"];
				B1[2] = A->parameters[0]["Bz"];

				B2[0] = B->parameters[0]["Bx"];
				B2[1] = B->parameters[0]["By"];
				B2[2] = B->parameters[0]["Bz"];

				Eigen::Vector3d J = n.cross(B2 - B1);

				J = J * 4.3614;

				bool llk = true;
				if (i->center[0][0] * 0.0891029508867553 + i->center[0][1] * 0.7044237408557898 + i->center[0][2] * (-0.7041646522383865) < 0.0)
				{
					J = -J;
					llk = false;
				}

				for (auto& j : i->yzels)
				{
					cc[0] = j->coord[0][0];
					cc[1] = j->coord[0][1];
					cc[2] = j->coord[0][2];

					Eigen::Vector3d JJ;
					Eigen::Vector3d PP;

					PP << 0.0891029508867553, 0.7044237408557898, -0.7041646522383865;
					PP = PP * 25.0;
					if (llk == false)
					{
						PP *= -1.0;
					}
					
					JJ = PP - cc;


					if (JJ.norm() > 7.0)
					{
						JJ = J;
					}

					fout << cc[0] << " " << cc[1] << " " << cc[2] << " " <<
						polar_angle(cc[1], cc[2]) << " " << polar_angle(cc[0], norm2(0.0, cc[1], cc[2])) << " " <<
						JJ[0] << " " << JJ[1] << " " << JJ[2] << " " << J.norm() << endl;
				}
			}


			for (int k = 0; k < this->Gran_TS.size(); k++)
			{
				fout << 4 * k + 1 << " " << 4 * k + 2 << " " << 4 * k + 3 << " " << 4 * k + 4 << endl;
			}


			fout.close();
		}

		// BS
		if (false)
		{
			name_f = "BS_J.txt";
			fout.open(name_f);
			fout << "TITLE = HP  VARIABLES = x, y, z, r, phi, the, Jx, Jy, Jz, |J|" << endl;
			fout << "ZONE T=HP, N = " << this->Gran_BS.size() * 4 << ", E = " << this->Gran_BS.size() << ", F=FEPOINT, ET=quadrilateral" << endl;

			for (const auto& i : this->Gran_BS)
			{
				auto A = i->cells[0];
				auto B = i->cells[1];
				Eigen::Vector3d n, B1, B2, cc;

				n[0] = i->normal[0][0];
				n[1] = i->normal[0][1];
				n[2] = i->normal[0][2];

				B1[0] = A->parameters[0]["Bx"];
				B1[1] = A->parameters[0]["By"];
				B1[2] = A->parameters[0]["Bz"];

				B2[0] = B->parameters[0]["Bx"];
				B2[1] = B->parameters[0]["By"];
				B2[2] = B->parameters[0]["Bz"];

				Eigen::Vector3d J = n.cross(B2 - B1);

				J = J * 4.3614;

				for (auto& j : i->yzels)
				{
					cc[0] = j->coord[0][0];
					cc[1] = j->coord[0][1];
					cc[2] = j->coord[0][2];

					if (n[0] > 0.2)
					{
						fout << cc[0] << " " << cc[1] << " " << cc[2] << " " << norm2(cc[0], cc[1], cc[2]) << " " <<
							polar_angle(cc[1], cc[2]) << " " << polar_angle(cc[0], norm2(0.0, cc[1], cc[2])) << " " <<
							J[0] << " " << J[1] << " " << J[2] << " " << J.norm() << endl;
					}
					else
					{
						fout << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 
							polar_angle(cc[1], cc[2]) << " " << polar_angle(cc[0], norm2(0.0, cc[1], cc[2])) << " " <<
							J[0] << " " << J[1] << " " << J[2] << " " << J.norm() << endl;
					}
				}
			}


			for (int k = 0; k < this->Gran_BS.size(); k++)
			{
				fout << 4 * k + 1 << " " << 4 * k + 2 << " " << 4 * k + 3 << " " << 4 * k + 4 << endl;
			}


			fout.close();
		}
	}
	else if (alg == 12)
	{
		// Вычислим ротор в центре каждой ячейки
		this->Edges_create();
		//this->Culc_usual_rotors_in_cell();
		this->Culc_usual_rotors_in_cell_2();
		//this->Culc_usual_rotors_in_cell_from_interpol();


		// Надо улучшить ротеры вблизи разрывов
		for (auto& gr : this->Gran_TS)
		{
			auto C1 = gr->cells[0];
			auto C3 = gr->cells[1];
			auto C2 = gr->cells_TVD[0];
			auto C4 = gr->cells_TVD[1];

			C1->parameters[0]["rotB_x"] = C2->parameters[0]["rotB_x"];
			C1->parameters[0]["rotB_y"] = C2->parameters[0]["rotB_y"];
			C1->parameters[0]["rotB_z"] = C2->parameters[0]["rotB_z"];

			C3->parameters[0]["rotB_x"] = C4->parameters[0]["rotB_x"];
			C3->parameters[0]["rotB_y"] = C4->parameters[0]["rotB_y"];
			C3->parameters[0]["rotB_z"] = C4->parameters[0]["rotB_z"];


			if (false)
			{
				C3->parameters[0]["rotB_x"] = C1->parameters[0]["rotB_x"];
				C3->parameters[0]["rotB_y"] = C1->parameters[0]["rotB_y"];
				C3->parameters[0]["rotB_z"] = C1->parameters[0]["rotB_z"];
			}
			else
			{
				C1->parameters[0]["rotB_x"] = C3->parameters[0]["rotB_x"];
				C1->parameters[0]["rotB_y"] = C3->parameters[0]["rotB_y"];
				C1->parameters[0]["rotB_z"] = C3->parameters[0]["rotB_z"];
			}
		}

		for (auto& gr : this->Gran_HP)
		{
			auto C1 = gr->cells[0];
			auto C2 = gr->cells_TVD[0];

			auto C3 = gr->cells[1];
			auto C4 = gr->cells_TVD[1];

			C1->parameters[0]["rotB_x"] = C2->parameters[0]["rotB_x"];
			C1->parameters[0]["rotB_y"] = C2->parameters[0]["rotB_y"];
			C1->parameters[0]["rotB_z"] = C2->parameters[0]["rotB_z"];

			C3->parameters[0]["rotB_x"] = C4->parameters[0]["rotB_x"];
			C3->parameters[0]["rotB_y"] = C4->parameters[0]["rotB_y"];
			C3->parameters[0]["rotB_z"] = C4->parameters[0]["rotB_z"];

			// Для внутреннего слоя

			if (norm2(C4->center[0][0], C4->center[0][1], C4->center[0][2]) >
				norm2(C2->center[0][0], C2->center[0][1], C2->center[0][2]))
			{
				if (false) // внутри
				{
					C3->parameters[0]["rotB_x"] = C1->parameters[0]["rotB_x"];
					C3->parameters[0]["rotB_y"] = C1->parameters[0]["rotB_y"];
					C3->parameters[0]["rotB_z"] = C1->parameters[0]["rotB_z"];
				}
				else
				{
					C1->parameters[0]["rotB_x"] = C3->parameters[0]["rotB_x"];
					C1->parameters[0]["rotB_y"] = C3->parameters[0]["rotB_y"];
					C1->parameters[0]["rotB_z"] = C3->parameters[0]["rotB_z"];
				}
			}
		}


		this->Save_for_interpolate("For_intertpolate_0059-.bin", false);
		Interpol SS = Interpol("For_intertpolate_0059-.bin");

		cout << "AAA" << endl;

		this->Tecplot_print_2D(&SS, 0.0177656909751554, 0.7057402284561816, 0.7082479157489927, -0.00001, "_IHG_meridional_", false,
			Eigen::Vector3d(-0.9958639688067080, 0.0756169508599243, -0.0503689624193315),
			Eigen::Vector3d(0.0891029508867553, 0.7044237408557898, -0.7041646522383865),
			Eigen::Vector3d(0.0, 0.0, 0.0));

		this->Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_XY_plane_", false,
			Eigen::Vector3d(1.0, 0.0, 0.0),
			Eigen::Vector3d(0.0, 1.0, 0.0),
			Eigen::Vector3d(0.0, 0.0, 0.0));

		this->Tecplot_print_2D(&SS, 0.0, 1.0, 0.0, -0.00001, "_XZ_plane_", false,
			Eigen::Vector3d(1.0, 0.0, 0.0),
			Eigen::Vector3d(0.0, 0.0, 1.0),
			Eigen::Vector3d(0.0, 0.0, 0.0));


		// Рисует тетраэдры в текплот
		if (false)
		{
			this->Save_for_interpolate_one_zone_only("For_intertpolate_work.bin", Type_cell::Zone_3);
			Interpol SS = Interpol("For_intertpolate_work.bin");

			ofstream fout;
			string name_f = "3D_setka_J.txt";
			fout.open(name_f);
			fout << "TITLE = HP  VARIABLES = x, y, z, Jx, Jy, Jz, |J|" << endl;

			fout << "ZONE T=\"Tetrahedra\", N=" << SS.points_1.size()
				<< ", E=" << SS.Delone_1->number_of_finite_cells()
				<< ", DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON\n";

			for (size_t i = 0; i < SS.points_1.size(); ++i)
			{
				const Point& p = SS.points_1[i].first;
				Int_point* params = SS.Cells_1[i];

				fout << p.x() << " " << p.y() << " " << p.z() << " " <<
					4.15368 * params->parameters["rotB_x"] << " " << 4.15368 * params->parameters["rotB_y"] << " " 
					<< 4.15368 * params->parameters["rotB_z"] << " " <<
					4.15368 * norm2(params->parameters["rotB_x"], params->parameters["rotB_y"], params->parameters["rotB_z"]) << endl;
			}

			// Выводим коннективность тетраэдров
			for (Delaunay::Cell_iterator cit = SS.Delone_1->cells_begin();
				cit != SS.Delone_1->cells_end(); ++cit)
			{
				// Пропускаем бесконечные ячейки
				if (SS.Delone_1->is_infinite(cit))
				{
					continue;
				}

				// Для каждой вершины тетраэдра находим ее глобальный индекс
				std::vector<size_t> vertex_indices;
				for (int i = 0; i < 4; ++i)
				{
					size_t vertex_point = cit->vertex(i)->info() + 1;
					vertex_indices.push_back(vertex_point);
				}

				// Записываем индексы вершин тетраэдра
				fout << vertex_indices[0] << " " << vertex_indices[1] << " "
					<< vertex_indices[2] << " " << vertex_indices[3] << "\n";
			}

			fout.close();
		}

		// Трассируем линии тока (вокруг TS раньше было, с HP так не работает
		if (false)
		{
			std::ofstream file("I_IHS.txt");
			file << "VARIABLES = X, Y, Z, I" << std::endl;
			int line_count = 0;

			this->Save_for_interpolate("For_intertpolate_work.bin", false);
			Interpol SS = Interpol("For_intertpolate_work.bin");

			unsigned int NN = 0;

			cout << "Start trasser" << endl;

			#pragma omp parallel for schedule(dynamic)
			for (auto& gr : this->Gran_HP)
			{
				bool bb;
				Cell* CC, * prev;
				prev = nullptr;
				std::unordered_map<string, double> parameters;
				int my_N;

				#pragma omp critical (dsds1) 
				{
					NN++;
					my_N = NN;
					if (NN % 100 == 0)
					{
						cout << "Trasser:  " << NN << "  from " << this->Gran_TS.size() << endl;
					}
				}

				if (my_N % 10 != 0) continue;

				if (gr->center[0][0] < 0) continue;

				double x, y, z, v;
				x = gr->center[0][0] + 0.1 * gr->normal[0][0];
				y = gr->center[0][1] + 0.1 * gr->normal[0][1];
				z = gr->center[0][2] + 0.1 * gr->normal[0][2];

				bb = SS.Get_param(x, y, z, parameters);

				if (bb == false)
				{
					cout << "Error bb  wergfwe4fr3gvbt4tevrgtefetw" << endl;
					continue;
				}

				if (parameters["rotB_x"] * x + parameters["rotB_y"] * y + parameters["rotB_z"] * z > 0.0)
				{
					v = 1.0;
				}
				else
				{
					v = -1.0;
				}

				std::vector<std::vector<double>> line;

				unsigned int kk = 0;
				while (true)
				{
					kk++;
					if (kk > 1000000) break;
					double nnn = 4.15368 * norm2(parameters["rotB_x"], parameters["rotB_y"], parameters["rotB_z"]);
					line.push_back({ x, y, z, nnn });

					x = x + 0.02 * v * parameters["rotB_x"] / nnn;
					y = y + 0.02 * v * parameters["rotB_y"] / nnn;
					z = z + 0.02 * v * parameters["rotB_z"] / nnn;

					CC = nullptr;
					CC = Find_cell_point(x, y, z, 0, prev);
					if (CC == nullptr) break;
					if (CC->type != Type_cell::Zone_3) break;

					bb = SS.Get_param(x, y, z, parameters);
					if (bb == false) break;
				}

#pragma omp critical (dsds2) 
				{
					int size_l = 0;
					for (const auto& point : line)
					{
						if (point[0] < -52.24) continue;
						size_l++;
					}
					size_l = line.size();

					file << "ZONE T=\"Line" << line_count++ << "\" I=" << size_l << " F=POINT" << std::endl;
					for (const auto& point : line)
					{
						//if (point[0] < -52.24) continue;
						file << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << std::endl;
					}
				}

			}

			file.close();
		}

		// Трассируем линию тока
		if (false)
		{
			std::ofstream file("Line_1.txt");
			file << "VARIABLES = X, Y, Z, I" << std::endl;
			int line_count = 0;

			this->Save_for_interpolate("For_intertpolate_work.bin", false);
			Interpol SS = Interpol("For_intertpolate_work.bin");

			unsigned int NN = 0;

			cout << "Start trasser" << endl;

			bool bb;
			Cell* CC, * prev;
			prev = nullptr;
			std::unordered_map<string, double> parameters;
			int my_N;



			double x, y, z, v;
			x = 50.0;
			y = 0.0;
			z = 0.0;

			bb = SS.Get_param(x, y, z, parameters);

			if (bb == false)
			{
				cout << "Error bb  wergfwe4fr3gvbt4tevrgtefetw" << endl;
				return;
			}

			v = 1.0;

			std::vector<std::vector<double>> line;

			unsigned int kk = 0;
			while (true)
			{
				kk++;
				if (kk > 1000000) break;
				double nnn = 4.15368 * norm2(parameters["rotB_x"], parameters["rotB_y"], parameters["rotB_z"]);
				line.push_back({ x, y, z, nnn });

				x = x + 0.02 * v * parameters["rotB_x"] / nnn;
				y = y + 0.02 * v * parameters["rotB_y"] / nnn;
				z = z + 0.02 * v * parameters["rotB_z"] / nnn;

				CC = nullptr;
				CC = Find_cell_point(x, y, z, 0, prev);
				if (CC == nullptr) break;
				if (CC->type != Type_cell::Zone_3) break;

				bb = SS.Get_param(x, y, z, parameters);
				if (bb == false) break;
			}

			int size_l = 0;
			size_l = line.size();

			file << "ZONE T=\"Line" << line_count++ << "\" I=" << size_l << " F=POINT" << std::endl;
			for (const auto& point : line)
			{
				file << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << std::endl;
			}

			

			file.close();
		}


	}
	else if (alg == 13)
	{
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		cout << "Create Setka Smc" << endl;
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		Setka Smc = Setka("SDK_40_2D_Setka.bin", "SDK_40_krug_setka.bin", 40);
		Smc.name = "Mini_for_MK";

		cout << "Create SI_main" << endl;
		// Из основной сетки создаём интерполяционную сетку
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		cout << "Move Setka Smc" << endl;
		// Двигаем поверхности вспомогательной сетки к поверхностям основной
		Smc.Move_to_surf(&SI_main);
		// Точно задаём положение внутренней границы сетки
		Smc.geo->R0 = Smc.phys_param->R_0;

		// Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
		Smc.auto_set_luch_geo_parameter(0, true);
		// Настраиваем новую сетку (также как и основную)   [обязательно]
		if (true)
		{
			// Считаем объёмы, площади и другие геометрические характеристики
			Smc.Calculating_measure(0);
			Smc.Calculating_measure(1);

			// Задаём граничные грани
			Smc.Init_boundary_grans();
		}

		// Заполним сетку МК значениями плазмы из основной сетки (чтобы вместо интерполяции в МК использовать значения в центрах ячеек - так быстрее)
		// переинтерполяция
		if (true)
		{
			Smc.PereInterpolate(&SI_main, false);
		}

		Smc.Test_geometr();


		Smc.Download_cell_MK_parameters(Smc.phys_param->MK_file, -10);

		Smc.Print_SpSm(17.0, 0.0, 0.0);
		Smc.Print_SpSm(10.0, 0.0, 0.0);
		Smc.Print_SpSm(15.0, 0.0, 0.0);
		Smc.Print_SpSm(17.0, 0.0, 0.0);
		Smc.Print_SpSm(20.0, 0.0, 0.0);
		Smc.Print_SpSm(22.0, 0.0, 0.0);
		Smc.Print_SpSm(25.0, 0.0, 0.0);
		Smc.Print_SpSm(27.0, 0.0, 0.0);
		Smc.Print_SpSm(30.0, 0.0, 0.0);
		Smc.Print_SpSm(40.0, 0.0, 0.0);
	}
	else if (alg == 14)
	{
		this->Save_for_interpolate("For_intertpolate_0059-.bin", false);
		Interpol SS = Interpol("For_intertpolate_0059-.bin");

		this->Tecplot_print_2D_for_HCS_potencial_1_zone(&SS, 0.0177656909751554, 0.7057402284561816, 0.7082479157489927, -0.00001, "_IHG_meridional_HCS_", false,
			Eigen::Vector3d(-0.9958639688067080, 0.0756169508599243, -0.0503689624193315),
			Eigen::Vector3d(0.0891029508867553, 0.7044237408557898, -0.7041646522383865),
			Eigen::Vector3d(0.0, 0.0, 0.0));
	}
	else if (alg == 15)
	{
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		unsigned int N_p = 90; // Количество точек на экваторе
		unsigned int N_l = 360; // Число слоёв
		unsigned int N_step = 3; // Число шагов по времени до создания новых точек

		// Вектор для хранения всех слоев
		std::vector<std::vector<Eigen::Vector3d>> all_layers;
		all_layers.reserve(N_l);

		double ddt = 0.0004 / N_step;
		// Глобальный цикл
		for (int step = 0; step < N_l; step++)
		{
			cout << "step = " << step << "   from: " << N_l << endl;
			std::vector<Eigen::Vector3d> points;
			points.reserve(N_p);

			// Создаём точки
			for (int i = 0; i < N_p; ++i)
			{
				// Угол в экваториальной плоскости
				double phi = 2.0 * const_pi * i / N_p;

				// Создаем точку на экваторе (z=0)
				Eigen::Vector3d point(this->phys_param->R_0 * cos(phi), this->phys_param->R_0 * sin(phi), 0.0);

				// Первое вращение: на угол alpha вокруг оси X
				Eigen::AngleAxisd rotation1(const_pi / 18.0, Eigen::Vector3d::UnitX());
				point = rotation1 * point;

				// Второе вращение: на угол beta вокруг оси Z
				Eigen::AngleAxisd rotation2(step * ddt * N_step * 2.0 * const_pi / 0.03843666, Eigen::Vector3d::UnitZ());
				point = rotation2 * point;

				Eigen::Vector3d point2 = this->phys_param->Matr * point;

				points.push_back(point2);
			}
			all_layers.push_back(points);


			// Теперь передвигаем точки
			for (int j = 0; j < N_step; ++j)
			{

				#pragma omp parallel for
				for (int i = 0; i < (int)all_layers.size(); ++i)
				{
					std::unordered_map<string, double> parameters;
					std::array<Cell_handle, 6> prev_cell;
					std::array<Cell_handle, 6> next_cell;
					for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();


					for (int j = 0; j < (int)all_layers[i].size(); ++j)
					{
						Eigen::Vector3d& point = all_layers[i][j];

						double r = 1.0;
						if (point.norm() < 3.0 * this->phys_param->R_0) r = 3.0 * this->phys_param->R_0 / point.norm();
						bool fine_int = SI_main.Get_param(r * point(0), r * point(1), r * point(2), parameters, prev_cell, next_cell);

						if (fine_int == false)
						{
							fine_int = SI_main.Get_param(r * point(0) * 0.997, r * point(1) * 0.999, r * point(2) * 0.999, parameters, prev_cell, next_cell);
							if (fine_int == false)
							{
								cout << "erorr euiegh87eg8ferg" << endl;
								exit(-1);
							}
						}

						for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i];

						point(0) += parameters["Vx"] * ddt;
						point(1) += parameters["Vy"] * ddt;
						point(2) += parameters["Vz"] * ddt;
					}
				}
			}
		}

		// Запись в Техплот
		if (true)
		{
			std::unordered_map<string, double> parameters;

			ofstream tecfile("HCS_3d.txt");
			if (!tecfile.is_open()) 
			{
				cerr << "Error jgiofueh9fgh3e489fergл " << endl;
				return;
			}

			tecfile << "TITLE = \"Heliospheric Current Sheet Surface\"" << endl;
			tecfile << "VARIABLES = \"X\", \"Y\", \"Z\", \"|J|\", \"Jx\", \"Jy\", \"Jz\"" << endl;

			unsigned int N_quads;
			N_quads = (N_l - 1) * N_p;  // Замкнутая поверхность

			// Записываем зону с четырехугольниками
			tecfile << "ZONE T=\"Surface\", N=" << N_l * N_p
				<< ", E=" << N_quads
				<< ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;

			// Записываем все точки с номером слоя как переменную
			for (unsigned int layer = 0; layer < N_l; ++layer) 
			{
				for (unsigned int point = 0; point < N_p; ++point)
				{
					const Eigen::Vector3d& p = all_layers[layer][point];

					double r = 1.0;
					if (p.norm() < 3.0 * this->phys_param->R_0) r = 3.0 * this->phys_param->R_0 / p.norm();
					bool fine_int = SI_main.Get_param(r * p(0), r * p(1), r * p(2), parameters);
					if (fine_int == false)
					{
						fine_int = SI_main.Get_param(r * p(0) * 0.997, r * p(1) * 0.999, r * p(2) * 0.999, parameters);
						if (fine_int == false)
						{
							cout << "erorr euiegh87eg8ferg" << endl;
							exit(-1);
						}
					}

					vector<Eigen::Vector3d> neighbors;

					int pp, pm;
					pp = point + 1;
					pm = point - 1;
					if (pp >= N_p) pp = 0;
					if (pm < 0) pm = N_p - 1;

					neighbors.push_back(all_layers[layer][pp]);
					neighbors.push_back(all_layers[layer][pm]);

					if (layer < N_l - 1)
					{
						neighbors.push_back(all_layers[layer + 1][pp]);
						neighbors.push_back(all_layers[layer + 1][point]);
						neighbors.push_back(all_layers[layer + 1][pm]);
					}

					if (layer > 0)
					{
						neighbors.push_back(all_layers[layer - 1][pp]);
						neighbors.push_back(all_layers[layer - 1][point]);
						neighbors.push_back(all_layers[layer - 1][pm]);
					}

					/*cout << "Start -----------------------" << endl;
					cout << "p = " << p[0] << " " << p[1] << " " << p[2] << endl;
					for (auto& aa : neighbors)
					{
						cout << "neighbor = " << aa[0] << " " << aa[1] << " " << aa[2] << endl;
					}*/
					Eigen::Vector3d nn = computeSurfaceNormal(p, neighbors);
					//cout << "End ------------------------- " << endl;

					if (nn[0] * 3.8759783635738505 + nn[1] * 30.64243272722684 + nn[2] * -30.631162372369808 < 0) nn = nn * -1.0;

					const double dim_j = 1.743;

					Eigen::Vector3d BB(parameters["Bx"], parameters["By"], parameters["Bz"]);
					Eigen::Vector3d jj = 2.0 * nn.cross(BB);

					tecfile << p.x() << " " << p.y() << " " << p.z()
						<< " " << 2.0 * norm2(parameters["Bx"], parameters["By"], parameters["Bz"]) * dim_j << " " 
						<< jj[0] << " " << jj[1] << " " << jj[2] << endl;
				}
			}

			// Записываем коннективити (соединения)
			// Tecplot использует 1-индексацию
			for (unsigned int layer = 0; layer < N_l - 1; ++layer) 
			{
				for (unsigned int point = 0; point < N_p - 1; ++point) 
				{
					unsigned int idx1 = layer * N_p + point + 1;       // Текущий слой, текущая точка
					unsigned int idx2 = layer * N_p + point + 1 + 1;   // Текущий слой, следующая точка
					unsigned int idx3 = (layer + 1) * N_p + point + 1 + 1; // Следующий слой, следующая точка
					unsigned int idx4 = (layer + 1) * N_p + point + 1;     // Следующий слой, текущая точка

					tecfile << idx1 << " " << idx2 << " " << idx3 << " " << idx4 << endl;
				}

				// Если поверхность замкнута, добавляем четырехугольник между последней и первой точкой
				if (true) 
				{
					unsigned int point = N_p - 1;
					unsigned int idx1 = layer * N_p + point + 1;       // Текущий слой, последняя точка
					unsigned int idx2 = layer * N_p + 0 + 1;           // Текущий слой, первая точка
					unsigned int idx3 = (layer + 1) * N_p + 0 + 1;     // Следующий слой, первая точка
					unsigned int idx4 = (layer + 1) * N_p + point + 1; // Следующий слой, последняя точка

					tecfile << idx1 << " " << idx2 << " " << idx3 << " " << idx4 << endl;
				}
			}

			tecfile.close();

		}

	}
	else if (alg == 16)
	{
		this->Set_MK_Zone();


		// this->MK_Grans[zone_MK - 1].size();
		// cell->MK_zone = 2;
		double Volume = 0.0;
		double E_B = 0.0;
		double E_B_pot = 0.0;
		double E_int = 0.0;
		double E_kin = 0.0;

		const int N = this->MK_Grans[2 - 1].size();   // Число фиктивных источников
		const int M = N;                              // Число граничных точек

		MatrixXd fict_points(N, 3);
		MatrixXd bnd_points(M, 3);
		MatrixXd bnd_normals(M, 3);
		VectorXd bnd_Bn(M);

		//    @param fict_points Координаты фиктивных источников(N x 3)
		//	* @param bnd_points Координаты граничных точек(M x 3)
		//	* @param bnd_normals Векторы нормалей в граничных точках(M x 3)
		//	* @param bnd_Bn Значения нормальной компоненты Bn в граничных точках(M)
		//	* @param use_svd Использовать SVD(true) или QR(false) для решения.SVD устойчивее к плохой обусловленности.
		//	* @return VectorXd Вектор амплитуд источников q(N)

		int i = 0;
		for (auto& gr : this->MK_Grans[2 - 1])
		{
			Cell* A, * B;
			double normal = 1.0;
			if (gr->cells[0]->MK_zone != 2)
			{
				A = gr->cells[0];
				B = gr->cells[1];
				normal = -1.0;
			}
			else
			{
				A = gr->cells[1];
				B = gr->cells[0];
			}
			// A - снаружи
			// B - внутри


			/*fict_points(i, 0) = A->center[0][0];
			fict_points(i, 1) = A->center[0][1];
			fict_points(i, 2) = A->center[0][2];*/

			fict_points(i, 0) = gr->center[0][0];
			fict_points(i, 1) = gr->center[0][1];
			fict_points(i, 2) = gr->center[0][2];

			bnd_normals(i, 0) = gr->normal[0][0] * normal;
			bnd_normals(i, 1) = gr->normal[0][1] * normal;
			bnd_normals(i, 2) = gr->normal[0][2] * normal;

			if (gr->type2 == Type_Gran_surf::HP)
			{
				/*bnd_points(i, 0) = gr->center[0][0];
				bnd_points(i, 1) = gr->center[0][1];
				bnd_points(i, 2) = gr->center[0][2];*/

				bnd_points(i, 0) = B->center[0][0];
				bnd_points(i, 1) = B->center[0][1];
				bnd_points(i, 2) = B->center[0][2];

				bnd_Bn(i) = 0.0;
			}
			else
			{
				bnd_points(i, 0) = B->center[0][0];
				bnd_points(i, 1) = B->center[0][1];
				bnd_points(i, 2) = B->center[0][2];

				bnd_Bn(i) = (B->parameters[0]["Bx"] * gr->normal[0][0] +
					B->parameters[0]["By"] * gr->normal[0][1] +
					B->parameters[0]["Bz"] * gr->normal[0][2]) * normal;
			}

			i++;
		}

		// Решение системы
		cout << "Solve " << endl;
		VectorXd q = solveMFS(fict_points, bnd_points, bnd_normals, bnd_Bn, false);
		cout << "End Solve " << endl;



		for (auto& cc : this->All_Cell)
		{
			if (cc->MK_zone != 2) continue;

			Vector3d test_point(cc->center[0][0], cc->center[0][1], cc->center[0][2]);
			double psi;
			Vector3d B_pot;
			computeField(q, fict_points, test_point, psi, B_pot);

			Volume += cc->volume[0];
			E_B += kvv(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) / (8.0 * const_pi) * cc->volume[0];
			E_B_pot += kv(B_pot.norm()) / (8.0 * const_pi) * cc->volume[0];
			E_int += (cc->parameters[0]["p"]) / (this->phys_param->gamma - 1.0) * cc->volume[0];
			E_kin += cc->parameters[0]["rho"] * kvv(cc->parameters[0]["Vx"], cc->parameters[0]["Vy"], cc->parameters[0]["Vz"]) / (2.0) * cc->volume[0];
		}

		cout << "E = " << endl;
		cout << "E_B = " << E_B/ Volume << endl;
		cout << "E_B_pot = " << E_B_pot/ Volume << endl;
		cout << "E_B - E_B_pot = " << (E_B - E_B_pot)/ Volume << endl;
		cout << "E_int = " << E_int/ Volume << endl;
		cout << "E_kin = " << E_kin/ Volume << endl;
		cout << "Volume = " << Volume << endl;
	}
	else if (alg == 17)
	{
		this->Save_for_interpolate("For_intertpolate_0059-.bin", false);
		Interpol SS = Interpol("For_intertpolate_0059-.bin");

		this->Tecplot_print_2D_for_HCS_potencial_1_zone(&SS, 0.0, 0.0, 1.0, -0.00001, "_IHG_meridional_potencial_HCS_", false,
			Eigen::Vector3d(-0.9958639688067077, 0.07561695085992419, -0.05036896241933166),
			Eigen::Vector3d(0.08910295088675518, 0.7044237408557894, -0.7041646522383864),
			Eigen::Vector3d(0.0, 0.0, 0.0));

		this->Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_IHG_meridional_", false,
			Eigen::Vector3d(-0.9958639688067077, 0.07561695085992419, -0.05036896241933166),
			Eigen::Vector3d(0.08910295088675518, 0.7044237408557894, -0.7041646522383864),
			Eigen::Vector3d(0.0, 0.0, 0.0));
	}
	else if (alg == 18)
	{
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		unsigned int N_p = 90; // Количество точек на экваторе
		unsigned int N_l = 360; // Число слоёв
		unsigned int N_step = 3; // Число шагов по времени до создания новых точек

		// Вектор для хранения всех слоев
		std::vector<std::vector<Eigen::Vector3d>> all_layers;
		all_layers.reserve(N_l);

		double ddt = 0.0004 / N_step;
		// Глобальный цикл
		for (int step = 0; step < N_l; step++)
		{
			cout << "step = " << step << "   from: " << N_l << endl;
			std::vector<Eigen::Vector3d> points;
			points.reserve(N_p);

			// Создаём точки
			for (int i = 0; i < N_p; ++i)
			{
				// Угол в экваториальной плоскости
				double phi = 2.0 * const_pi * i / N_p;

				// Создаем точку на экваторе (z=0)
				Eigen::Vector3d point(this->phys_param->R_0 * cos(phi), this->phys_param->R_0 * sin(phi), 0.0);

				// Первое вращение: на угол alpha вокруг оси X
				Eigen::AngleAxisd rotation1(const_pi / 18.0, Eigen::Vector3d::UnitX());
				point = rotation1 * point;

				// Второе вращение: на угол beta вокруг оси Z
				Eigen::AngleAxisd rotation2(step * ddt * N_step * 2.0 * const_pi / 0.03843666, Eigen::Vector3d::UnitZ());
				point = rotation2 * point;

				Eigen::Vector3d point2 = this->phys_param->Matr * point;

				points.push_back(point2);
			}
			all_layers.push_back(points);


			// Теперь передвигаем точки
			for (int j = 0; j < N_step; ++j)
			{
				#pragma omp parallel for
				for (int i = 0; i < (int)all_layers.size(); ++i)
				{
					std::unordered_map<string, double> parameters;
					std::array<Cell_handle, 6> prev_cell;
					std::array<Cell_handle, 6> next_cell;
					for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();


					for (int j = 0; j < (int)all_layers[i].size(); ++j)
					{
						Eigen::Vector3d& point = all_layers[i][j];

						if (polar_angle(point(0), norm2(0.0, point(1), point(2))) > const_pi / 9.0) continue;

						double r = 1.0;
						if (point.norm() < 3.0 * this->phys_param->R_0) r = 3.0 * this->phys_param->R_0 / point.norm();
						bool fine_int = SI_main.Get_param(r * point(0), r * point(1), r * point(2), parameters, prev_cell, next_cell);

						if (fine_int == false)
						{
							fine_int = SI_main.Get_param(r * point(0) * 0.997, r * point(1) * 0.999, r * point(2) * 0.999, parameters, prev_cell, next_cell);
							if (fine_int == false)
							{
								cout << "erorr euiegh87eg8ferg" << endl;
								exit(-1);
							}
						}

						for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i];

						point(0) += parameters["Vx"] * ddt;
						point(1) += parameters["Vy"] * ddt;
						point(2) += parameters["Vz"] * ddt;
					}
				}
			}
		}

		// Теперь двигаем точки
		cout << "Move only" << endl;

		bool bj = false;
		while (bj == false)
		{
			bj = true;

			#pragma omp parallel for
			for (int i = 0; i < (int)all_layers.size(); ++i)
			{
				std::unordered_map<string, double> parameters;
				std::array<Cell_handle, 6> prev_cell;
				std::array<Cell_handle, 6> next_cell;
				for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();
				Cell* prevC = nullptr;


				for (int j = 0; j < (int)all_layers[i].size(); ++j)
				{
					Eigen::Vector3d& point = all_layers[i][j];

					if (polar_angle(point(0), norm2(0.0, point(1), point(2))) > const_pi / 9.0) continue;

					double r = 1.0;
					bool fine_int = SI_main.Get_param(r * point(0), r * point(1), r * point(2), parameters, prev_cell, next_cell);

					if (fine_int == false)
					{
						fine_int = SI_main.Get_param(r * point(0) * 0.997, r * point(1) * 0.999, r * point(2) * 0.999, parameters, prev_cell, next_cell);
						if (fine_int == false)
						{
							cout << "erorr euiegh87eg8ferg" << endl;
							exit(-1);
						}
					}

					for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i];

					point(0) += parameters["Vx"] * ddt;
					point(1) += parameters["Vy"] * ddt;
					point(2) += parameters["Vz"] * ddt;

					if (i == (int)all_layers[i].size() - 1)
					{
						auto CCC = this->Find_cell_point(point(0), point(1), point(2), 0, prevC);
						if (CCC->type == Type_cell::Zone_1)
						{
							bj = false;
						}
					}
				}
			}
		}




		// Запись в Техплот
		if (true)
		{
			std::unordered_map<string, double> parameters;

			ofstream tecfile("HCS_3d.txt");
			if (!tecfile.is_open())
			{
				cerr << "Error jgiofueh9fgh3e489fergл " << endl;
				return;
			}

			tecfile << "TITLE = \"Heliospheric Current Sheet Surface\"" << endl;
			tecfile << "VARIABLES = \"X\", \"Y\", \"Z\", \"|J|\", \"Jx\", \"Jy\", \"Jz\"" << endl;

			unsigned int N_quads;
			N_quads = (N_l - 1) * N_p;  // Замкнутая поверхность

			// Записываем зону с четырехугольниками
			tecfile << "ZONE T=\"Surface\", N=" << N_l * N_p
				<< ", E=" << N_quads
				<< ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;

			// Записываем все точки с номером слоя как переменную
			for (unsigned int layer = 0; layer < N_l; ++layer)
			{
				for (unsigned int point = 0; point < N_p; ++point)
				{
					const Eigen::Vector3d& p = all_layers[layer][point];

					if (polar_angle(p(0), norm2(0.0, p(1), p(2))) > const_pi / 9.0)
					{
						tecfile << 0.0 << " " << 0.0 << " " << 0.0
							<< " " << 0.0 << " "
							<< 0.0 << " " << 0.0 << " " << 0.0 << endl;
						continue;
					}

					double r = 1.0;
					if (p.norm() < 3.0 * this->phys_param->R_0) r = 3.0 * this->phys_param->R_0 / p.norm();
					bool fine_int = SI_main.Get_param(r * p(0), r * p(1), r * p(2), parameters);
					if (fine_int == false)
					{
						fine_int = SI_main.Get_param(r * p(0) * 0.997, r * p(1) * 0.999, r * p(2) * 0.999, parameters);
						if (fine_int == false)
						{
							cout << "erorr euiegh87eg8ferg" << endl;
							exit(-1);
						}
					}

					vector<Eigen::Vector3d> neighbors;

					int pp, pm;
					pp = point + 1;
					pm = point - 1;
					if (pp >= N_p) pp = 0;
					if (pm < 0) pm = N_p - 1;

					neighbors.push_back(all_layers[layer][pp]);
					neighbors.push_back(all_layers[layer][pm]);

					if (layer < N_l - 1)
					{
						neighbors.push_back(all_layers[layer + 1][pp]);
						neighbors.push_back(all_layers[layer + 1][point]);
						neighbors.push_back(all_layers[layer + 1][pm]);
					}

					if (layer > 0)
					{
						neighbors.push_back(all_layers[layer - 1][pp]);
						neighbors.push_back(all_layers[layer - 1][point]);
						neighbors.push_back(all_layers[layer - 1][pm]);
					}

					/*cout << "Start -----------------------" << endl;
					cout << "p = " << p[0] << " " << p[1] << " " << p[2] << endl;
					for (auto& aa : neighbors)
					{
						cout << "neighbor = " << aa[0] << " " << aa[1] << " " << aa[2] << endl;
					}*/
					Eigen::Vector3d nn = computeSurfaceNormal(p, neighbors);
					//cout << "End ------------------------- " << endl;

					if (nn[0] * 3.8759783635738505 + nn[1] * 30.64243272722684 + nn[2] * -30.631162372369808 < 0) nn = nn * -1.0;

					const double dim_j = 1.74456;

					Eigen::Vector3d BB(parameters["Bx"], parameters["By"], parameters["Bz"]);
					Eigen::Vector3d jj = 2.0 * nn.cross(BB);

					tecfile << p.x() << " " << p.y() << " " << p.z()
						<< " " << 2.0 * norm2(parameters["Bx"], parameters["By"], parameters["Bz"]) * dim_j << " "
						<< jj[0] << " " << jj[1] << " " << jj[2] << endl;
				}
			}

			// Записываем коннективити (соединения)
			// Tecplot использует 1-индексацию
			for (unsigned int layer = 0; layer < N_l - 1; ++layer)
			{
				for (unsigned int point = 0; point < N_p - 1; ++point)
				{
					unsigned int idx1 = layer * N_p + point + 1;       // Текущий слой, текущая точка
					unsigned int idx2 = layer * N_p + point + 1 + 1;   // Текущий слой, следующая точка
					unsigned int idx3 = (layer + 1) * N_p + point + 1 + 1; // Следующий слой, следующая точка
					unsigned int idx4 = (layer + 1) * N_p + point + 1;     // Следующий слой, текущая точка

					tecfile << idx1 << " " << idx2 << " " << idx3 << " " << idx4 << endl;
				}

				// Если поверхность замкнута, добавляем четырехугольник между последней и первой точкой
				if (true)
				{
					unsigned int point = N_p - 1;
					unsigned int idx1 = layer * N_p + point + 1;       // Текущий слой, последняя точка
					unsigned int idx2 = layer * N_p + 0 + 1;           // Текущий слой, первая точка
					unsigned int idx3 = (layer + 1) * N_p + 0 + 1;     // Следующий слой, первая точка
					unsigned int idx4 = (layer + 1) * N_p + point + 1; // Следующий слой, последняя точка

					tecfile << idx1 << " " << idx2 << " " << idx3 << " " << idx4 << endl;
				}
			}

			tecfile.close();

		}

		}
	else if (alg == 19)
	{
		this->Set_MK_Zone();

		this->phys_param->param_names.push_back("gr_x");
		this->phys_param->param_names.push_back("gr_y");
		this->phys_param->param_names.push_back("gr_z");

		for (auto& cc : this->All_Cell)
		{
			cc->parameters[0]["phi_1"] = 0.0;
			cc->parameters[0]["phi_2"] = 0.0;

			cc->parameters[0]["gr_x"] = 0.0;
			cc->parameters[0]["gr_y"] = 0.0;
			cc->parameters[0]["gr_z"] = 0.0;
		}

		std::ifstream in_file("cell_data.bin", std::ios::binary);

		if (in_file.is_open()) 
		{
			for (auto& cc : this->All_Cell) 
			{
				double gr_x, gr_y, gr_z;

				in_file.read(reinterpret_cast<char*>(&gr_x), sizeof(double));
				in_file.read(reinterpret_cast<char*>(&gr_y), sizeof(double));
				in_file.read(reinterpret_cast<char*>(&gr_z), sizeof(double));

				cc->parameters[0]["gr_x"] = gr_x;
				cc->parameters[0]["gr_y"] = gr_y;
				cc->parameters[0]["gr_z"] = gr_z;
			}
			in_file.close();
		}


		for (auto& cc : this->All_Cell)
		{
			if (cc->MK_zone == 2) cc->MK_zone = 20;

			if (cc->MK_zone == 3)
			{
				if(cc->center[0][0] > this->geo->L6) cc->MK_zone = 20;
			}
		}



		int ZONE_now = 20;  // В какой зоне сейчас считаем


		int step = 0;
		double dphi_ = 100.0;
		double x_max = 0.0;
		double y_max = 0.0;
		double z_max = 0.0;

		while (true)
		{
			step++;
			if (step % 100 == 0)
			{
				cout << "step = " << step << "  nevyazka = " << dphi_ << endl;
				cout << "x_max = " << x_max << "  y_max = " << y_max << "  z_max = " << z_max << endl;
			}

			// Пробегаемся по граням считаем нужные потоки
#pragma omp parallel for schedule(dynamic)
			for (auto& gr : this->All_Gran)
			{
				if (gr->cells.size() != 2) continue;
				if (gr->cells[0]->MK_zone != ZONE_now && gr->cells[1]->MK_zone != ZONE_now) continue;

				Cell* A, * B;
				A = gr->cells[0];
				B = gr->cells[1];

				if (gr->type2 == Type_Gran_surf::HP)
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					gr->parameters["dPhi"] = 0.0;
					gr->parameters["dSS"] = gr->area[0] / dl;
					continue;
				}


				if (gr->cells[0]->MK_zone == ZONE_now && gr->cells[1]->MK_zone == ZONE_now)
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					double dl1 = norm2(A->center[0][0] - gr->center[0][0], A->center[0][1] - gr->center[0][1], A->center[0][2] - gr->center[0][2]);
					double dl2 = norm2(B->center[0][0] - gr->center[0][0], B->center[0][1] - gr->center[0][1], B->center[0][2] - gr->center[0][2]);

					if (step < 1000)
					{
						gr->parameters["dPhi"] = (B->parameters[0]["phi_1"] - A->parameters[0]["phi_1"]) / dl * gr->area[0];
					}
					else
					{
						gr->parameters["dPhi"] = ((A->parameters[0]["gr_x"] * dl2 + B->parameters[0]["gr_x"] * dl1) / (dl1 + dl2) * gr->normal[0][0] +
							(A->parameters[0]["gr_y"] * dl2 + B->parameters[0]["gr_y"] * dl1) / (dl1 + dl2) * gr->normal[0][1] +
							(A->parameters[0]["gr_z"] * dl2 + B->parameters[0]["gr_z"] * dl1) / (dl1 + dl2) * gr->normal[0][2]) * gr->area[0];
					}
					
					//gr->parameters["dPhi"] = (B->parameters[0]["phi_1"] - A->parameters[0]["phi_1"]) / dl * gr->area[0];
					gr->parameters["dSS"] = gr->area[0] / dl;
				}
				else
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					if (gr->cells[0]->MK_zone == ZONE_now)
					{
						gr->parameters["dPhi"] = -(A->parameters[0]["Bx"] * gr->normal[0][0] + A->parameters[0]["By"] * gr->normal[0][1] +
							A->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
					}
					else
					{
						gr->parameters["dPhi"] = -(B->parameters[0]["Bx"] * gr->normal[0][0] + B->parameters[0]["By"] * gr->normal[0][1] +
							B->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
					}
					gr->parameters["dSS"] = gr->area[0] / dl;
				}
			}

			// Пробегаемся по ячейкам
#pragma omp parallel for schedule(dynamic)
			for (auto& cc : this->All_Cell)
			{
				if (cc->MK_zone != ZONE_now) continue;
				double dPhi = 0.0;
				double dSS = 0.0;

				for (auto& gr : cc->grans)
				{
					dSS += gr->parameters["dSS"];
					if (gr->cells[0]->number == cc->number)
					{
						dPhi += gr->parameters["dPhi"];
					}
					else
					{
						dPhi -= gr->parameters["dPhi"];
					}
				}

				cc->parameters[0]["phi_2"] = cc->parameters[0]["phi_1"] + 1.0 / dSS * dPhi;
				
			}

			dphi_ = 0.0;


#pragma omp parallel for schedule(dynamic)
			for (auto& cc : this->All_Cell)
			{
				if (cc->MK_zone != ZONE_now) continue;

				double ddd = fabs(cc->parameters[0]["phi_1"] - cc->parameters[0]["phi_2"]);
				cc->parameters[0]["phi_1"] = cc->parameters[0]["phi_2"];
				if (ddd > dphi_)
				{
					#pragma omp critical (sdfs)
					{
						if (ddd > dphi_)
						{
							dphi_ = ddd;
							x_max = cc->center[0][0];
							y_max = cc->center[0][1];
							z_max = cc->center[0][2];
						}
					}
				}
			}

			/*if (dphi_ < 0.1)
			{
				cout << "Yspex" << endl;
				break;
			}*/
			if (step > 10000)
			{
				cout << "Yspex" << endl;
				break;
			}

			if (true)
			{
				double Volume = 0.0;
				double E_B = 0.0;
				double E_B_pot = 0.0;
				double dE_B = 0.0;
				double E_int = 0.0;
				double E_kin = 0.0;

				double Volume2 = 0.0;
				double E_B2 = 0.0;
				double E_B_pot2 = 0.0;
				double dE_B2 = 0.0;
				double E_int2 = 0.0;
				double E_kin2 = 0.0;

#pragma omp parallel for schedule(dynamic)
				for (auto& cc : this->All_Cell)
				{
					if (cc->MK_zone != ZONE_now) continue;

					double gr_x = 0.0;
					double gr_y = 0.0;
					double gr_z = 0.0;
					double BB = 0.0;
					bool bnb = false;

					for (auto& gr : cc->grans)
					{
						if (gr->cells[0]->MK_zone != ZONE_now) bnb = true;
						if (gr->cells[1]->MK_zone != ZONE_now) bnb = true;

						double d1 = norm2(gr->cells[0]->center[0][0] - gr->center[0][0], gr->cells[0]->center[0][1] - gr->center[0][1],
							gr->cells[0]->center[0][2] - gr->center[0][2]);
						double d2 = norm2(gr->cells[1]->center[0][0] - gr->center[0][0], gr->cells[1]->center[0][1] - gr->center[0][1],
							gr->cells[1]->center[0][2] - gr->center[0][2]);

						BB = (gr->cells[0]->parameters[0]["phi_1"] * d2 + gr->cells[1]->parameters[0]["phi_1"] * d1) / (d1 + d2);
						if (gr->cells[0]->number == cc->number)
						{
							gr_x += BB * gr->normal[0][0] * gr->area[0];
							gr_y += BB * gr->normal[0][1] * gr->area[0];
							gr_z += BB * gr->normal[0][2] * gr->area[0];
						}
						else
						{
							gr_x -= BB * gr->normal[0][0] * gr->area[0];
							gr_y -= BB * gr->normal[0][1] * gr->area[0];
							gr_z -= BB * gr->normal[0][2] * gr->area[0];
						}
					}

					gr_x = gr_x / cc->volume[0];
					gr_y = gr_y / cc->volume[0];
					gr_z = gr_z / cc->volume[0];

					if (bnb == true)
					{
						//gr_x = -1.0; //-cc->parameters[0]["Bx"];
						//gr_y = 0.0;  //-cc->parameters[0]["By"];
						//gr_z = 0.0;  // -cc->parameters[0]["Bz"];

						gr_x = -cc->parameters[0]["Bx"];
						gr_y = -cc->parameters[0]["By"];
						gr_z = -cc->parameters[0]["Bz"];
					}

					cc->parameters[0]["gr_x"] = gr_x;
					cc->parameters[0]["gr_y"] = gr_y;
					cc->parameters[0]["gr_z"] = gr_z;

					if (step % 1000 == 0)
					{
						if (cc->center[0][0] >= 0.0)
						{
							#pragma omp critical (dwedewfwef)
							{
								Volume += cc->volume[0];
								E_B += kvv(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) / (8.0 * const_pi) * cc->volume[0];
								E_B_pot += kvv(gr_x, gr_y, gr_z) / (8.0 * const_pi) * cc->volume[0];
								E_int += (cc->parameters[0]["p"]) / (this->phys_param->gamma - 1.0) * cc->volume[0];
								E_kin += cc->parameters[0]["rho"] * kvv(cc->parameters[0]["Vx"], cc->parameters[0]["Vy"], cc->parameters[0]["Vz"]) / (2.0) * cc->volume[0];
								dE_B += kv(norm2(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) -
									norm2(gr_x, gr_y, gr_z)) / (8.0 * const_pi);
							}
						}
						else
						{
							#pragma omp critical (dwedewfwef2)
							{
								Volume2 += cc->volume[0];
								E_B2 += kvv(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) / (8.0 * const_pi) * cc->volume[0];
								E_B_pot2 += kvv(gr_x, gr_y, gr_z) / (8.0 * const_pi) * cc->volume[0];
								E_int2 += (cc->parameters[0]["p"]) / (this->phys_param->gamma - 1.0) * cc->volume[0];
								E_kin2 += cc->parameters[0]["rho"] * kvv(cc->parameters[0]["Vx"], cc->parameters[0]["Vy"], cc->parameters[0]["Vz"]) / (2.0) * cc->volume[0];
								dE_B2 += kv(norm2(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) -
									norm2(gr_x, gr_y, gr_z)) / (8.0 * const_pi);
							}
						}
					}
				}

				if (step % 1000 == 0)
				{
					cout << "----------------------------------------" << endl;
					cout << "rho E = " << (E_B + E_int + E_kin) / Volume << endl;
					cout << "E_B % = " << E_B / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "dE_B % = " << dE_B / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_B_pot % = " << E_B_pot / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_free % = " << (E_B - E_B_pot) / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_int % = " << E_int / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_kin % = " << E_kin / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "Volume = " << Volume << endl;
					cout << "----------------------------------------" << endl;
					cout << "rho E = " << (E_B2 + E_int2 + E_kin2) / Volume2 << endl;
					cout << "E_B % = " << E_B2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "dE_B % = " << dE_B2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "E_B_pot % = " << E_B_pot2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "E_free % = " << (E_B2 - E_B_pot2) / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "E_int % = " << E_int2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "E_kin % = " << E_kin2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "Volume = " << Volume2 << endl;
					cout << "----------------------------------------" << endl;
				}
			}


		}
		
		// Запись значений в бинарный файл
		if (false)
		{
			std::ofstream out_file("cell_data.bin", std::ios::binary);

			if (out_file.is_open())
			{
				for (auto& cc : this->All_Cell)
				{
					// Записываем в файл
					double gr_x = cc->parameters[0]["gr_x"];
					double gr_y = cc->parameters[0]["gr_y"];
					double gr_z = cc->parameters[0]["gr_z"];

					out_file.write(reinterpret_cast<const char*>(&gr_x), sizeof(double));
					out_file.write(reinterpret_cast<const char*>(&gr_y), sizeof(double));
					out_file.write(reinterpret_cast<const char*>(&gr_z), sizeof(double));
				}
				out_file.close();
			}
		}

		this->Save_for_interpolate("For_intertpolate_0059-.bin", false);
		Interpol SS = Interpol("For_intertpolate_0059-.bin");

		this->Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_2d_(0, 0, 1, 0)_");

		this->Tecplot_print_2D(&SS, 0.0891029508867553, 0.7044237408557898, -0.7041646522383865, -0.00001, "_IHG_polar_", false,
			Eigen::Vector3d(-0.9958639688067080, 0.0756169508599243, -0.0503689624193315),
			Eigen::Vector3d(0.0177656909751554, 0.7057402284561816, 0.7082479157489927),
			Eigen::Vector3d(0.0, 0.0, 0.0));
		
	}
	else if (alg == 20)
	{
		this->Set_MK_Zone();

		for (auto& cc : this->All_Cell)
		{
			cc->parameters[0]["phi_1"] = 0.0;
			cc->parameters[0]["phi_2"] = 0.0;
		}

		int ZONE_now = 1;  // В какой зоне сейчас считаем


		int step = 0;
		double dphi_ = 100.0;
		while (true)
		{
			step++;
			if (step % 100 == 0)
			{
				cout << "step = " << step << "  nevyazka = " << dphi_ << endl;
			}

			//cout << "A" << endl;

			// Пробегаемся по граням считаем нужные потоки
#pragma omp parallel for schedule(dynamic)
			for (auto& gr : this->All_Gran)
			{
				if (gr->cells.size() != 2 && gr->cells[0]->MK_zone != ZONE_now) continue;

				if (gr->cells.size() != 2)
				{
					Cell* A;
					A = gr->cells[0];
					double dl = norm2(A->center[0][0] - gr->center[0][0], A->center[0][1] - gr->center[0][1], A->center[0][2] - gr->center[0][2]);
					gr->parameters["dPhi"] = (A->parameters[0]["Bx"] * gr->normal[0][0] + A->parameters[0]["By"] * gr->normal[0][1] +
						A->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
					gr->parameters["dSS"] = gr->area[0] / dl;
					continue;
				}
				

				if (gr->cells[0]->MK_zone != ZONE_now && gr->cells[1]->MK_zone != ZONE_now) continue;

				Cell* A, * B;
				A = gr->cells[0];
				B = gr->cells[1];


				if (gr->cells[0]->MK_zone == ZONE_now && gr->cells[1]->MK_zone == ZONE_now)
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					gr->parameters["dPhi"] = (B->parameters[0]["phi_1"] - A->parameters[0]["phi_1"]) / dl * gr->area[0];
					gr->parameters["dSS"] = gr->area[0] / dl;
				}
				else
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					if (gr->cells[0]->MK_zone == ZONE_now)
					{
						gr->parameters["dPhi"] = (A->parameters[0]["Bx"] * gr->normal[0][0] + A->parameters[0]["By"] * gr->normal[0][1] +
							A->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
					}
					else
					{
						gr->parameters["dPhi"] = (B->parameters[0]["Bx"] * gr->normal[0][0] + B->parameters[0]["By"] * gr->normal[0][1] +
							B->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
					}
					gr->parameters["dSS"] = gr->area[0] / dl;
				}
			}

			//cout << "B" << endl;

			// Пробегаемся по ячейкам
#pragma omp parallel for schedule(dynamic)
			for (auto& cc : this->All_Cell)
			{
				if (cc->MK_zone != ZONE_now) continue;
				double dPhi = 0.0;
				double dSS = 0.0;

				for (auto& gr : cc->grans)
				{
					dSS += gr->parameters["dSS"];
					if (gr->cells[0]->number == cc->number)
					{
						dPhi += gr->parameters["dPhi"];
					}
					else
					{
						dPhi -= gr->parameters["dPhi"];
					}
				}

				cc->parameters[0]["phi_2"] = cc->parameters[0]["phi_1"] + 1.0 / dSS * dPhi;

			}

			dphi_ = 0.0;

			//cout << "C" << endl;
#pragma omp parallel for schedule(dynamic)
			for (auto& cc : this->All_Cell)
			{
				if (cc->MK_zone != ZONE_now) continue;

				double ddd = fabs(cc->parameters[0]["phi_1"] - cc->parameters[0]["phi_2"]);
				cc->parameters[0]["phi_1"] = cc->parameters[0]["phi_2"];
				if (ddd > dphi_)
				{
#pragma omp critical (sdfs)
					{
						if (ddd > dphi_)
						{
							dphi_ = ddd;
						}
					}
				}
			}

			/*if (dphi_ < 0.1)
			{
				cout << "Yspex" << endl;
				break;
			}*/
			if (step > 300000)
			{
				cout << "Yspex" << endl;
				break;
			}

			if (step % 1000 == 0)
			{
				double Volume = 0.0;
				double E_B = 0.0;
				double E_B_pot = 0.0;
				double E_int = 0.0;
				double E_kin = 0.0;

				double Volume2 = 0.0;
				double E_B2 = 0.0;
				double E_B_pot2 = 0.0;
				double E_int2 = 0.0;
				double E_kin2 = 0.0;

				for (auto& cc : this->All_Cell)
				{
					if (cc->MK_zone != ZONE_now) continue;

					double gr_x = 0.0;
					double gr_y = 0.0;
					double gr_z = 0.0;
					double BB = 0.0;
					bool bnb = false;

					for (auto& gr : cc->grans)
					{
						if (gr->cells.size() < 2)
						{
							bnb = true;
							break;
						}

						if (gr->cells[0]->MK_zone != ZONE_now) bnb = true;
						if (gr->cells[1]->MK_zone != ZONE_now) bnb = true;

						double d1 = norm2(gr->cells[0]->center[0][0] - gr->center[0][0], gr->cells[0]->center[0][1] - gr->center[0][1],
							gr->cells[0]->center[0][2] - gr->center[0][2]);
						double d2 = norm2(gr->cells[1]->center[0][0] - gr->center[0][0], gr->cells[1]->center[0][1] - gr->center[0][1],
							gr->cells[1]->center[0][2] - gr->center[0][2]);

						BB = (gr->cells[0]->parameters[0]["phi_1"] * d2 + gr->cells[1]->parameters[0]["phi_1"] * d1) / (d1 + d2);
						if (gr->cells[0]->number == cc->number)
						{
							gr_x += BB * gr->normal[0][0] * gr->area[0];
							gr_y += BB * gr->normal[0][1] * gr->area[0];
							gr_z += BB * gr->normal[0][2] * gr->area[0];
						}
						else
						{
							gr_x -= BB * gr->normal[0][0] * gr->area[0];
							gr_y -= BB * gr->normal[0][1] * gr->area[0];
							gr_z -= BB * gr->normal[0][2] * gr->area[0];
						}
					}

					gr_x = gr_x / cc->volume[0];
					gr_y = gr_y / cc->volume[0];
					gr_z = gr_z / cc->volume[0];

					if (bnb == true)
					{
						gr_x = cc->parameters[0]["Bx"];
						gr_y = cc->parameters[0]["By"];
						gr_z = cc->parameters[0]["Bz"];
					}

					double r = norm2(cc->center[0][0], cc->center[0][1], cc->center[0][2]);

					if (r * 4.21132 <= 20.0)
					{
						Volume += cc->volume[0];
						E_B += kvv(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) / (8.0 * const_pi) * cc->volume[0];
						E_B_pot += kvv(gr_x, gr_y, gr_z) / (8.0 * const_pi) * cc->volume[0];
						E_int += (cc->parameters[0]["p"]) / (this->phys_param->gamma - 1.0) * cc->volume[0];
						E_kin += cc->parameters[0]["rho"] * kvv(cc->parameters[0]["Vx"], cc->parameters[0]["Vy"], cc->parameters[0]["Vz"]) / (2.0) * cc->volume[0];
					}
					else
					{
						Volume2 += cc->volume[0];
						E_B2 += kvv(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) / (8.0 * const_pi) * cc->volume[0];
						E_B_pot2 += kvv(gr_x, gr_y, gr_z) / (8.0 * const_pi) * cc->volume[0];
						E_int2 += (cc->parameters[0]["p"]) / (this->phys_param->gamma - 1.0) * cc->volume[0];
						E_kin2 += cc->parameters[0]["rho"] * kvv(cc->parameters[0]["Vx"], cc->parameters[0]["Vy"], cc->parameters[0]["Vz"]) / (2.0) * cc->volume[0];
					}
				}

				cout << "----------------------------------------" << endl;
				cout << "rho E = " << (E_B + E_int + E_kin) / Volume << endl;
				cout << "E_B % = " << E_B / (E_B + E_int + E_kin) * 100.0 << endl;
				cout << "E_B_pot % = " << E_B_pot / (E_B + E_int + E_kin) * 100.0 << endl;
				cout << "E_free % = " << (E_B - E_B_pot) / (E_B + E_int + E_kin) * 100.0 << endl;
				cout << "E_int % = " << E_int / (E_B + E_int + E_kin) * 100.0 << endl;
				cout << "E_kin % = " << E_kin / (E_B + E_int + E_kin) * 100.0 << endl;
				cout << "Volume = " << Volume << endl;
				cout << "----------------------------------------" << endl;
				cout << "rho E = " << (E_B2 + E_int2 + E_kin2) / Volume2 << endl;
				cout << "E_B % = " << E_B2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
				cout << "E_B_pot % = " << E_B_pot2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
				cout << "E_free % = " << (E_B2 - E_B_pot2) / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
				cout << "E_int % = " << E_int2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
				cout << "E_kin % = " << E_kin2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
				cout << "Volume = " << Volume2 << endl;
				cout << "----------------------------------------" << endl;
			}


		}



		}
	else if (alg == 21)
	{
		this->Set_MK_Zone();

		for (auto& cc : this->All_Cell)
		{
			if (cc->MK_zone == 5)
			{
				if (cc->center[0][0] > this->geo->L6) cc->MK_zone = 4;
			}
		}

		this->phys_param->param_names.push_back("gr_x");
		this->phys_param->param_names.push_back("gr_y");
		this->phys_param->param_names.push_back("gr_z");
		
		for (auto& cc : this->All_Cell)
		{
			cc->parameters[0]["phi_1"] = 0.0;
			cc->parameters[0]["phi_2"] = 0.0;

			cc->parameters[0]["gr_x"] = 0.0;
			cc->parameters[0]["gr_y"] = 0.0;
			cc->parameters[0]["gr_z"] = 0.0;
		}

		int ZONE_now = 4;  // В какой зоне сейчас считаем

		// Найдём номер ячейки в которой надо проверить BЖ

		int NMNM; 
		if (this->Gran_BS[1]->cells_TVD[0]->MK_zone == ZONE_now) NMNM = this->Gran_BS[1]->cells_TVD[0]->number;
		if (this->Gran_BS[1]->cells_TVD[1]->MK_zone == ZONE_now) NMNM = this->Gran_BS[1]->cells_TVD[1]->number;

		int NMNM2;
		int jk = 0;
		for (auto& cc : this->All_Cell)
		{
			if (cc->MK_zone != ZONE_now) continue;
			jk++;
			if (jk == 1000)
			{
				NMNM2 = cc->number;
				break;
			}
		}


		int step = 0;
		double dphi_ = 100.0;
		double x_max = 0.0;
		double y_max = 0.0;
		double z_max = 0.0;

		while (true)
		{
			step++;
			if (step % 100 == 0)
			{
				cout << "step = " << step << "  nevyazka = " << dphi_ << endl;
				cout << "x_max = " << x_max << "  y_max = " << y_max << "  z_max = " << z_max << endl;
			}

			//cout << "A" << endl;

			// Пробегаемся по граням считаем нужные потоки
#pragma omp parallel for schedule(dynamic)
			for (auto& gr : this->All_Gran)
			{
				if (gr->cells.size() != 2) continue;

				if (gr->cells[0]->MK_zone != ZONE_now && gr->cells[1]->MK_zone != ZONE_now) continue;

				Cell* A, * B;
				A = gr->cells[0];
				B = gr->cells[1];

				if (gr->type2 == Type_Gran_surf::HP)
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					gr->parameters["dPhi"] = 0.0;
					gr->parameters["dSS"] = gr->area[0] / dl;
					continue;
				}


				if (gr->cells[0]->MK_zone == ZONE_now && gr->cells[1]->MK_zone == ZONE_now)
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					double dl1 = norm2(A->center[0][0] - gr->center[0][0], A->center[0][1] - gr->center[0][1], A->center[0][2] - gr->center[0][2]);
					double dl2 = norm2(B->center[0][0] - gr->center[0][0], B->center[0][1] - gr->center[0][1], B->center[0][2] - gr->center[0][2]);

					if (step < 1000)
					{
						gr->parameters["dPhi"] = (B->parameters[0]["phi_1"] - A->parameters[0]["phi_1"]) / dl * gr->area[0];
					}
					else
					{
						gr->parameters["dPhi"] = ((A->parameters[0]["gr_x"] * dl2 + B->parameters[0]["gr_x"] * dl1) / (dl1 + dl2) * gr->normal[0][0] +
							(A->parameters[0]["gr_y"] * dl2 + B->parameters[0]["gr_y"] * dl1) / (dl1 + dl2) * gr->normal[0][1] +
							(A->parameters[0]["gr_z"] * dl2 + B->parameters[0]["gr_z"] * dl1) / (dl1 + dl2) * gr->normal[0][2]) * gr->area[0];
					}


					gr->parameters["dSS"] = gr->area[0] / dl;
				}
				else
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					if (gr->cells[0]->MK_zone == ZONE_now)
					{
						gr->parameters["dPhi"] = -(A->parameters[0]["Bx"] * gr->normal[0][0] + A->parameters[0]["By"] * gr->normal[0][1] +
							A->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
						//gr->parameters["dPhi"] = -(1.0 * gr->normal[0][0]) * gr->area[0];
					}
					else
					{
						gr->parameters["dPhi"] = -(B->parameters[0]["Bx"] * gr->normal[0][0] + B->parameters[0]["By"] * gr->normal[0][1] +
							B->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
						//gr->parameters["dPhi"] = -(1.0 * gr->normal[0][0]) * gr->area[0];
					}
					gr->parameters["dSS"] = gr->area[0] / dl;
				}
			}

			//cout << "B" << endl;

			// Пробегаемся по ячейкам
#pragma omp parallel for schedule(dynamic)
			for (auto& cc : this->All_Cell)
			{
				if (cc->MK_zone != ZONE_now) continue;
				double dPhi = 0.0;
				double dSS = 0.0;

				for (auto& gr : cc->grans)
				{
					dSS += gr->parameters["dSS"];
					if (gr->cells[0]->number == cc->number)
					{
						dPhi += gr->parameters["dPhi"];
					}
					else
					{
						dPhi -= gr->parameters["dPhi"];
					}
				}

				cc->parameters[0]["phi_2"] = cc->parameters[0]["phi_1"] + 1.0 / dSS * dPhi;

			}

			dphi_ = 0.0;

			//cout << "C" << endl;
#pragma omp parallel for schedule(dynamic)
			for (auto& cc : this->All_Cell)
			{
				if (cc->MK_zone != ZONE_now) continue;

				double ddd = fabs(cc->parameters[0]["phi_1"] - cc->parameters[0]["phi_2"]);
				cc->parameters[0]["phi_1"] = cc->parameters[0]["phi_2"];
				if (ddd > dphi_)
				{
#pragma omp critical (sdfs)
					{
						if (ddd > dphi_)
						{
							dphi_ = ddd;
							x_max = cc->center[0][0];
							y_max = cc->center[0][1];
							z_max = cc->center[0][2];
						}
					}
				}
			}

			/*if (dphi_ < 0.1)
			{
				cout << "Yspex" << endl;
				break;
			}*/
			if (step > 7000)
			{
				cout << "Yspex" << endl;
				break;
			}

			if (true)
			{
				double Volume = 0.0;
				double E_B = 0.0;
				double E_B_pot = 0.0;
				double dE_B = 0.0;
				double E_int = 0.0;
				double E_kin = 0.0;

#pragma omp parallel for schedule(dynamic)
				for (auto& cc : this->All_Cell)
				{
					if (cc->MK_zone != ZONE_now) continue;

					double gr_x = 0.0;
					double gr_y = 0.0;
					double gr_z = 0.0;
					double BB = 0.0;
					bool bnb = false;

					for (auto& gr : cc->grans)
					{
						if (gr->cells[0]->MK_zone != ZONE_now) bnb = true;
						if (gr->cells[1]->MK_zone != ZONE_now) bnb = true;

						double d1 = norm2(gr->cells[0]->center[0][0] - gr->center[0][0], gr->cells[0]->center[0][1] - gr->center[0][1],
							gr->cells[0]->center[0][2] - gr->center[0][2]);
						double d2 = norm2(gr->cells[1]->center[0][0] - gr->center[0][0], gr->cells[1]->center[0][1] - gr->center[0][1],
							gr->cells[1]->center[0][2] - gr->center[0][2]);

						BB = (gr->cells[0]->parameters[0]["phi_1"] * d2 + gr->cells[1]->parameters[0]["phi_1"] * d1) / (d1 + d2);
						if (gr->cells[0]->number == cc->number)
						{
							gr_x += BB * gr->normal[0][0] * gr->area[0];
							gr_y += BB * gr->normal[0][1] * gr->area[0];
							gr_z += BB * gr->normal[0][2] * gr->area[0];
						}
						else
						{
							gr_x -= BB * gr->normal[0][0] * gr->area[0];
							gr_y -= BB * gr->normal[0][1] * gr->area[0];
							gr_z -= BB * gr->normal[0][2] * gr->area[0];
						}
					}

					gr_x = gr_x / cc->volume[0];
					gr_y = gr_y / cc->volume[0];
					gr_z = gr_z / cc->volume[0];

					if (bnb == true)
					{
						//gr_x = -1.0; //-cc->parameters[0]["Bx"];
						//gr_y = 0.0;  //-cc->parameters[0]["By"];
						//gr_z = 0.0;  // -cc->parameters[0]["Bz"];

						gr_x = -cc->parameters[0]["Bx"];
						gr_y = -cc->parameters[0]["By"];
						gr_z = -cc->parameters[0]["Bz"];
					}

					cc->parameters[0]["gr_x"] = gr_x;
					cc->parameters[0]["gr_y"] = gr_y;
					cc->parameters[0]["gr_z"] = gr_z;

					if (step % 1000 == 0)
					{
						if (cc->number == NMNM)
						{
							cout << cc->parameters[0]["Bx"] << " " << cc->parameters[0]["By"] << " " << cc->parameters[0]["Bz"] << endl;
							cout << -gr_x << " " << -gr_y << " " << -gr_z << endl;
						}

						if (cc->number == NMNM2)
						{
							cout << " = " << cc->parameters[0]["Bx"] << " " << cc->parameters[0]["By"] << " " << cc->parameters[0]["Bz"] << endl;
							cout << -gr_x << " " << -gr_y << " " << -gr_z << endl;
						}
						
						#pragma omp critical (dwedewfwef)
						{
							Volume += cc->volume[0];
							E_B += kvv(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) / (8.0 * const_pi) * cc->volume[0];
							E_B_pot += kvv(gr_x, gr_y, gr_z) / (8.0 * const_pi) * cc->volume[0];
							E_int += (cc->parameters[0]["p"]) / (this->phys_param->gamma - 1.0) * cc->volume[0];
							E_kin += cc->parameters[0]["rho"] * kvv(cc->parameters[0]["Vx"], cc->parameters[0]["Vy"], cc->parameters[0]["Vz"]) / (2.0) * cc->volume[0];
							dE_B += kv(norm2(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) -
								norm2(gr_x, gr_y, gr_z)) / (8.0 * const_pi);
						}
					}
					
				}


				if (step % 1000 == 0)
				{
					cout << "----------------------------------------" << endl;
					cout << "rho E = " << (E_B + E_int + E_kin) / Volume << endl;
					cout << "dE_B % = " << dE_B / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_B % = " << E_B / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_B_pot % = " << E_B_pot / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_free % = " << (E_B - E_B_pot) / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_int % = " << E_int / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_kin % = " << E_kin / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "Volume = " << Volume << endl;
					cout << "----------------------------------------" << endl;
				}
			}


		}

		for (auto& cc : this->All_Cell)
		{
			cc->parameters[0]["gr_x"] *= -1.0;
			cc->parameters[0]["gr_y"] *= -1.0;
			cc->parameters[0]["gr_z"] *= -1.0;
		}


		// Запись значений в бинарный файл
		std::ofstream out_file("cell_data.bin", std::ios::binary);

		if (out_file.is_open()) 
		{
			for (auto& cc : this->All_Cell) 
			{
				// Записываем в файл
				double gr_x = cc->parameters[0]["gr_x"];
				double gr_y = cc->parameters[0]["gr_y"];
				double gr_z = cc->parameters[0]["gr_z"];

				out_file.write(reinterpret_cast<const char*>(&gr_x), sizeof(double));
				out_file.write(reinterpret_cast<const char*>(&gr_y), sizeof(double));
				out_file.write(reinterpret_cast<const char*>(&gr_z), sizeof(double));
			}
			out_file.close();
		}

		this->Save_for_interpolate("For_intertpolate_0059-.bin", false);
		Interpol SS = Interpol("For_intertpolate_0059-.bin");
		this->Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_2d_(0, 0, 1, 0)_");

	}
	else if (alg == 22)
	{
		this->Set_MK_Zone();

		this->phys_param->param_names.push_back("gr_x");
		this->phys_param->param_names.push_back("gr_y");
		this->phys_param->param_names.push_back("gr_z");

		for (auto& cc : this->All_Cell)
		{
			cc->parameters[0]["phi_1"] = 0.0;
			cc->parameters[0]["phi_2"] = 0.0;

			cc->parameters[0]["gr_x"] = 0.0;
			cc->parameters[0]["gr_y"] = 0.0;
			cc->parameters[0]["gr_z"] = 0.0;
		}

		std::ifstream in_file("cell_data.bin", std::ios::binary);

		if (in_file.is_open())
		{
			for (auto& cc : this->All_Cell)
			{
				double gr_x, gr_y, gr_z;

				in_file.read(reinterpret_cast<char*>(&gr_x), sizeof(double));
				in_file.read(reinterpret_cast<char*>(&gr_y), sizeof(double));
				in_file.read(reinterpret_cast<char*>(&gr_z), sizeof(double));

				cc->parameters[0]["gr_x"] = gr_x;
				cc->parameters[0]["gr_y"] = gr_y;
				cc->parameters[0]["gr_z"] = gr_z;
			}
			in_file.close();
		}


		int ZONE_now = 1;  // В какой зоне сейчас считаем


		int step = 0;
		double dphi_ = 100.0;
		double x_max = 0.0;
		double y_max = 0.0;
		double z_max = 0.0;

		while (true)
		{
			step++;
			if (step % 100 == 0)
			{
				cout << "step = " << step << "  nevyazka = " << dphi_ << endl;
				cout << "x_max = " << x_max << "  y_max = " << y_max << "  z_max = " << z_max << endl;
			}

			// Пробегаемся по граням считаем нужные потоки
#pragma omp parallel for schedule(dynamic)
			for (auto& gr : this->All_Gran)
			{
				if (gr->cells.size() != 2 && gr->cells[0]->MK_zone != ZONE_now) continue;

				if (gr->cells.size() != 2)
				{
					Cell* A;
					A = gr->cells[0];
					double dl = norm2(A->center[0][0] - gr->center[0][0], A->center[0][1] - gr->center[0][1], A->center[0][2] - gr->center[0][2]);
					gr->parameters["dPhi"] = -(A->parameters[0]["Bx"] * gr->normal[0][0] + A->parameters[0]["By"] * gr->normal[0][1] +
						A->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
					gr->parameters["dSS"] = gr->area[0] / dl;
					continue;
				}


				if (gr->cells[0]->MK_zone != ZONE_now && gr->cells[1]->MK_zone != ZONE_now) continue;

				Cell* A, * B;
				A = gr->cells[0];
				B = gr->cells[1];

				if (gr->type2 == Type_Gran_surf::HP)
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					gr->parameters["dPhi"] = 0.0;
					gr->parameters["dSS"] = gr->area[0] / dl;
					continue;
				}


				if (gr->cells[0]->MK_zone == ZONE_now && gr->cells[1]->MK_zone == ZONE_now)
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					double dl1 = norm2(A->center[0][0] - gr->center[0][0], A->center[0][1] - gr->center[0][1], A->center[0][2] - gr->center[0][2]);
					double dl2 = norm2(B->center[0][0] - gr->center[0][0], B->center[0][1] - gr->center[0][1], B->center[0][2] - gr->center[0][2]);

					if (step < 1000)
					{
						gr->parameters["dPhi"] = (B->parameters[0]["phi_1"] - A->parameters[0]["phi_1"]) / dl * gr->area[0];
					}
					else
					{
						gr->parameters["dPhi"] = ((A->parameters[0]["gr_x"] * dl2 + B->parameters[0]["gr_x"] * dl1) / (dl1 + dl2) * gr->normal[0][0] +
							(A->parameters[0]["gr_y"] * dl2 + B->parameters[0]["gr_y"] * dl1) / (dl1 + dl2) * gr->normal[0][1] +
							(A->parameters[0]["gr_z"] * dl2 + B->parameters[0]["gr_z"] * dl1) / (dl1 + dl2) * gr->normal[0][2]) * gr->area[0];
					}

					//gr->parameters["dPhi"] = (B->parameters[0]["phi_1"] - A->parameters[0]["phi_1"]) / dl * gr->area[0];
					gr->parameters["dSS"] = gr->area[0] / dl;
				}
				else
				{
					double dl = norm2(A->center[0][0] - B->center[0][0], A->center[0][1] - B->center[0][1], A->center[0][2] - B->center[0][2]);
					if (gr->cells[0]->MK_zone == ZONE_now)
					{
						gr->parameters["dPhi"] = -(A->parameters[0]["Bx"] * gr->normal[0][0] + A->parameters[0]["By"] * gr->normal[0][1] +
							A->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
					}
					else
					{
						gr->parameters["dPhi"] = -(B->parameters[0]["Bx"] * gr->normal[0][0] + B->parameters[0]["By"] * gr->normal[0][1] +
							B->parameters[0]["Bz"] * gr->normal[0][2]) * gr->area[0];
					}
					gr->parameters["dSS"] = gr->area[0] / dl;
				}
			}

			// Пробегаемся по ячейкам
#pragma omp parallel for schedule(dynamic)
			for (auto& cc : this->All_Cell)
			{
				if (cc->MK_zone != ZONE_now) continue;
				double dPhi = 0.0;
				double dSS = 0.0;

				for (auto& gr : cc->grans)
				{
					dSS += gr->parameters["dSS"];
					if (gr->cells[0]->number == cc->number)
					{
						dPhi += gr->parameters["dPhi"];
					}
					else
					{
						dPhi -= gr->parameters["dPhi"];
					}
				}

				cc->parameters[0]["phi_2"] = cc->parameters[0]["phi_1"] + 1.0 / dSS * dPhi;

			}

			dphi_ = 0.0;


#pragma omp parallel for schedule(dynamic)
			for (auto& cc : this->All_Cell)
			{
				if (cc->MK_zone != ZONE_now) continue;

				double ddd = fabs(cc->parameters[0]["phi_1"] - cc->parameters[0]["phi_2"]);
				cc->parameters[0]["phi_1"] = cc->parameters[0]["phi_2"];
				if (ddd > dphi_)
				{
#pragma omp critical (sdfs)
					{
						if (ddd > dphi_)
						{
							dphi_ = ddd;
							x_max = cc->center[0][0];
							y_max = cc->center[0][1];
							z_max = cc->center[0][2];
						}
					}
				}
			}

			/*if (dphi_ < 0.1)
			{
				cout << "Yspex" << endl;
				break;
			}*/
			if (step > 10000)
			{
				cout << "Yspex" << endl;
				break;
			}

			if (true)
			{
				double Volume = 0.0;
				double E_B = 0.0;
				double E_B_pot = 0.0;
				double dE_B = 0.0;
				double E_int = 0.0;
				double E_kin = 0.0;

				double Volume2 = 0.0;
				double E_B2 = 0.0;
				double E_B_pot2 = 0.0;
				double dE_B2 = 0.0;
				double E_int2 = 0.0;
				double E_kin2 = 0.0;

#pragma omp parallel for schedule(dynamic)
				for (auto& cc : this->All_Cell)
				{
					if (cc->MK_zone != ZONE_now) continue;

					double gr_x = 0.0;
					double gr_y = 0.0;
					double gr_z = 0.0;
					double BB = 0.0;
					bool bnb = false;

					for (auto& gr : cc->grans)
					{
						if (gr->cells.size() < 2)
						{
							bnb = true;
							break;
						}

						if (gr->cells[0]->MK_zone != ZONE_now) bnb = true;
						if (gr->cells[1]->MK_zone != ZONE_now) bnb = true;

						double d1 = norm2(gr->cells[0]->center[0][0] - gr->center[0][0], gr->cells[0]->center[0][1] - gr->center[0][1],
							gr->cells[0]->center[0][2] - gr->center[0][2]);
						double d2 = norm2(gr->cells[1]->center[0][0] - gr->center[0][0], gr->cells[1]->center[0][1] - gr->center[0][1],
							gr->cells[1]->center[0][2] - gr->center[0][2]);

						BB = (gr->cells[0]->parameters[0]["phi_1"] * d2 + gr->cells[1]->parameters[0]["phi_1"] * d1) / (d1 + d2);
						if (gr->cells[0]->number == cc->number)
						{
							gr_x += BB * gr->normal[0][0] * gr->area[0];
							gr_y += BB * gr->normal[0][1] * gr->area[0];
							gr_z += BB * gr->normal[0][2] * gr->area[0];
						}
						else
						{
							gr_x -= BB * gr->normal[0][0] * gr->area[0];
							gr_y -= BB * gr->normal[0][1] * gr->area[0];
							gr_z -= BB * gr->normal[0][2] * gr->area[0];
						}
					}

					gr_x = gr_x / cc->volume[0];
					gr_y = gr_y / cc->volume[0];
					gr_z = gr_z / cc->volume[0];

					if (bnb == true)
					{
						//gr_x = -1.0; //-cc->parameters[0]["Bx"];
						//gr_y = 0.0;  //-cc->parameters[0]["By"];
						//gr_z = 0.0;  // -cc->parameters[0]["Bz"];

						gr_x = -cc->parameters[0]["Bx"];
						gr_y = -cc->parameters[0]["By"];
						gr_z = -cc->parameters[0]["Bz"];
					}

					cc->parameters[0]["gr_x"] = gr_x;
					cc->parameters[0]["gr_y"] = gr_y;
					cc->parameters[0]["gr_z"] = gr_z;

					if (step % 1000 == 0)
					{

						double r = norm2(cc->center[0][0], cc->center[0][1], cc->center[0][2]);
						if (r * 4.21132 <= 20.0)
						{
#pragma omp critical (dwedewfwef)
							{
								Volume += cc->volume[0];
								E_B += kvv(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) / (8.0 * const_pi) * cc->volume[0];
								E_B_pot += kvv(gr_x, gr_y, gr_z) / (8.0 * const_pi) * cc->volume[0];
								E_int += (cc->parameters[0]["p"]) / (this->phys_param->gamma - 1.0) * cc->volume[0];
								E_kin += cc->parameters[0]["rho"] * kvv(cc->parameters[0]["Vx"], cc->parameters[0]["Vy"], cc->parameters[0]["Vz"]) / (2.0) * cc->volume[0];
								dE_B += kv(norm2(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) -
									norm2(gr_x, gr_y, gr_z)) / (8.0 * const_pi);
							}
						}
						else
						{
#pragma omp critical (dwedewfwef2)
							{
								Volume2 += cc->volume[0];
								E_B2 += kvv(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) / (8.0 * const_pi) * cc->volume[0];
								E_B_pot2 += kvv(gr_x, gr_y, gr_z) / (8.0 * const_pi) * cc->volume[0];
								E_int2 += (cc->parameters[0]["p"]) / (this->phys_param->gamma - 1.0) * cc->volume[0];
								E_kin2 += cc->parameters[0]["rho"] * kvv(cc->parameters[0]["Vx"], cc->parameters[0]["Vy"], cc->parameters[0]["Vz"]) / (2.0) * cc->volume[0];
								dE_B2 += kv(norm2(cc->parameters[0]["Bx"], cc->parameters[0]["By"], cc->parameters[0]["Bz"]) -
									norm2(gr_x, gr_y, gr_z)) / (8.0 * const_pi);
							}
						}
					}
				}

				if (step % 1000 == 0)
				{
					cout << "----------------------------------------" << endl;
					cout << "rho E = " << (E_B + E_int + E_kin) / Volume << endl;
					cout << "E_B % = " << E_B / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "dE_B % = " << dE_B / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_B_pot % = " << E_B_pot / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_free % = " << (E_B - E_B_pot) / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_int % = " << E_int / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "E_kin % = " << E_kin / (E_B + E_int + E_kin) * 100.0 << endl;
					cout << "Volume = " << Volume << endl;
					cout << "----------------------------------------" << endl;
					cout << "rho E = " << (E_B2 + E_int2 + E_kin2) / Volume2 << endl;
					cout << "E_B % = " << E_B2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "dE_B % = " << dE_B2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "E_B_pot % = " << E_B_pot2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "E_free % = " << (E_B2 - E_B_pot2) / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "E_int % = " << E_int2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "E_kin % = " << E_kin2 / (E_B2 + E_int2 + E_kin2) * 100.0 << endl;
					cout << "Volume = " << Volume2 << endl;
					cout << "----------------------------------------" << endl;
				}
			}


		}

		// Запись значений в бинарный файл
		if (true)
		{
			std::ofstream out_file("cell_data.bin", std::ios::binary);

			if (out_file.is_open())
			{
				for (auto& cc : this->All_Cell)
				{
					// Записываем в файл
					double gr_x = cc->parameters[0]["gr_x"];
					double gr_y = cc->parameters[0]["gr_y"];
					double gr_z = cc->parameters[0]["gr_z"];

					out_file.write(reinterpret_cast<const char*>(&gr_x), sizeof(double));
					out_file.write(reinterpret_cast<const char*>(&gr_y), sizeof(double));
					out_file.write(reinterpret_cast<const char*>(&gr_z), sizeof(double));
				}
				out_file.close();
			}
		}

		this->Save_for_interpolate("For_intertpolate_0059-.bin", false);
		Interpol SS = Interpol("For_intertpolate_0059-.bin");

		this->Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_2d_(0, 0, 1, 0)_");

		this->Tecplot_print_2D(&SS, 0.0891029508867553, 0.7044237408557898, -0.7041646522383865, -0.00001, "_IHG_polar_", false,
			Eigen::Vector3d(-0.9958639688067080, 0.0756169508599243, -0.0503689624193315),
			Eigen::Vector3d(0.0177656909751554, 0.7057402284561816, 0.7082479157489927),
			Eigen::Vector3d(0.0, 0.0, 0.0));

			}
	else if (alg == 23)
	{
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		cout << "Create Setka Smc" << endl;
		// Создаём вспомогательную Монте-Карло сетку из файлов вспомогательных сеток
		Setka Smc = Setka("SDK_40_2D_Setka.bin", "SDK_40_krug_setka.bin", 40);

		cout << "Create SI_main" << endl;
		// Из основной сетки создаём интерполяционную сетку
		this->Save_for_interpolate("For_intertpolate_work.bin", false);
		Interpol SI_main = Interpol("For_intertpolate_work.bin");

		cout << "Move Setka Smc" << endl;
		// Двигаем поверхности вспомогательной сетки к поверхностям основной
		Smc.Move_to_surf(&SI_main);
		// Точно задаём положение внутренней границы сетки
		Smc.geo->R0 = Smc.phys_param->R_0;

		// Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
		Smc.auto_set_luch_geo_parameter(0, true);
		// Настраиваем новую сетку (также как и основную)   [обязательно]
		if (true)
		{
			// Считаем объёмы, площади и другие геометрические характеристики
			Smc.Calculating_measure(0);
			Smc.Calculating_measure(1);

			// Задаём граничные грани
			Smc.Init_boundary_grans();
		}

		// В сетке для MK очистим ненужные имена переменных 
		if (false)
		{
			Smc.phys_param->param_names.assign(Smc.phys_param->MK_param.begin(), Smc.phys_param->MK_param.end());
		}

		// Заполним сетку МК значениями плазмы из основной сетки (чтобы вместо интерполяции в МК использовать значения в центрах ячеек - так быстрее)
		// переинтерполяция
		if (true)
		{
			Smc.PereInterpolate(&SI_main, false);
		}

		Smc.Test_geometr();

		std::ofstream out("Sp_Sm_for_work_MK.bin", std::ios::binary);
		if (!out.is_open()) 
		{
			throw std::runtime_error("Cannot open file for writing");
		}

		size_t num_cells = Smc.All_Cell.size();
		out.write(reinterpret_cast<const char*>(&num_cells), sizeof(size_t));

		size_t num_ = Smc.phys_param->pui_nW;
		out.write(reinterpret_cast<const char*>(&num_), sizeof(size_t));

		// Загружаем S+ S- для всей сетки
		for (auto& A : Smc.All_Cell)
		{
			A->Init_S(2, Smc.phys_param->pui_nW);
			A->read_S_FromFile(Smc.phys_param->par_n_H_LISM);

			out.write(reinterpret_cast<const char*>(A->pui_Sm.data()), num_ * sizeof(double));

			// Всегда записываем как матрицу 2 x n
			int rows = A->pui_Sp.rows();

			// Записываем первую строку
			for (size_t j = 0; j < num_; ++j)
			{
				double val = A->pui_Sp(0, j);
				out.write(reinterpret_cast<const char*>(&val), sizeof(double));
			}

			// Первая строка матрицы
			if (rows == 2) 
			{
				for (size_t j = 0; j < num_; ++j)
				{
					double val = A->pui_Sp(1, j);
					out.write(reinterpret_cast<const char*>(&val), sizeof(double));
				}
			}
			else 
			{
				for (size_t j = 0; j < num_; ++j)
				{
					double zero = 0.0;
					out.write(reinterpret_cast<const char*>(&zero), sizeof(double));
				}
			}

			A->pui_Sm.resize(0);
			A->pui_Sp.resize(0, 0);
		}

		Smc.Print_SpSm(17.0, 0.0, 0.0);
		Smc.Print_SpSm(20.0, 0.0, 0.0);
		Smc.Print_SpSm(25.0, 0.0, 0.0);
		Smc.Print_SpSm(1.0, 0.0, 0.0);
		Smc.Print_SpSm(5.0, 0.0, 0.0);
		Smc.Print_SpSm(10.0, 0.0, 0.0);
		Smc.Print_SpSm(15.0, 0.0, 0.0);

		int vall = 148;
		out.write(reinterpret_cast<const char*>(&vall), sizeof(int));

		out.close();

		Smc.Save_for_interpolate("For_intertpolate_work_MK.bin", false);

		cout << "End - proverka" << endl;
		Interpol SS = Interpol("For_intertpolate_work_MK.bin");
		SS.Read_Sp_Sm("Sp_Sm_for_work_MK.bin");
		Cell_handle prev_cell = Cell_handle();
		Cell_handle next_cell;
		vector<double> mas_Sm(SS.pui_nW);
		vector<double> mas_Sp1(SS.pui_nW);
		vector<double> mas_Sp2(SS.pui_nW);
		bool b = SS.Get_Source(25.0, 0.0, 0.0, prev_cell, next_cell,
			mas_Sm, mas_Sp1, mas_Sp2);
		if (b == true)
		{
			for (int i = 0; i < SS.pui_nW; i++)
			{
				cout << i << " " << mas_Sm[i] << " " << mas_Sp1[i] << endl;
			}
		}
	}
	else if (alg == 24)
	{
		this->Save_for_interpolate("For_intertpolate_work3.bin", false);
		Interpol SS = Interpol("For_intertpolate_work3.bin");

		double dX = 2.0;  // минимум по 0.5, но лучше меньше
		double dY = 2.0;  // минимум по 0.5, но лучше меньше

		//double dX = 1.0;  // минимум по 0.5, но лучше меньше
		//double dY = 1.0;  // минимум по 0.5, но лучше меньше

		//double dX = 0.25;  // минимум по 0.5, но лучше меньше
		//double dY = 0.25;  // минимум по 0.5, но лучше меньше


		double dZ = 0.05;


		// Перевод в размерные единицы
		double ch_r = 0.0140668;          // расстояние для переведа в парсеки
		double ch_T = 39111.5;            // температура для перевода в Кельвины
		double ch_Halpha = 378663.0 * 0.6106;            // температура для перевода в Кельвины
		// последний коэффициент от того, что там    =   ne nH

		double XL = -285.0;
		double XR = 100.0;
		double YL = -213.0;
		double YR = 213.0;

		//double XL = -71.0;
		//double XR = 71.0;
		//double YL = -142.0;
		//double YR = 142.0;

		//double XL = 0.0;
		//double XR = 57.0;
		//double YL = -71.0;
		//double YR = 71.0;
		

		ofstream fout;
		fout.open("H_alpha_X_Y_2.12-2-v0.txt");

		ofstream fout2;
		fout2.open("soft_X-ray_X_Y.txt");

		ofstream fout3;
		fout3.open("hard_X-ray_X_Y.txt");

		const int NX = static_cast<int>(65.0 / dX + 0.5); // ~130


//#pragma omp parallel for schedule(dynamic)
		for (double X = XR; X > XL; X = X - dX)
		//for (double X = 120.2; X > -120.7; X = X - dX)
		//for (int iX = 0; iX < NX; ++iX)
		{
			//double X = 65.0 - iX * dX;

			//#pragma omp critical 
			//{
				cout << "Culk for X = " << X << endl;
			//}

			double ne, T;

			for (double Y = YL; Y < YR; Y = Y + dY)
			//for (double Y = 0.0; Y < 0.000001; Y = Y + dY)
			//for (double Y = -225.5; Y < 225.5; Y = Y + dY)
			{
				std::array<Cell_handle, 6> prev_cell;
				std::array<Cell_handle, 6> next_cell;
				for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();
				std::unordered_map<string, double> parameters;
				bool fine_int;

				double IH = 0.0;
				double I_soft_X_ray = 0.0;
				double I_hard_X_ray = 0.0;
				double a1, a2;
				fine_int = SS.Get_param(X, Y, 0.0, parameters, prev_cell, next_cell);
				//fine_int = SS.Get_param(X, 0.0, Y, parameters, prev_cell, next_cell);
				//fine_int = SS.Get_param(0.0, X, Y, parameters, prev_cell, next_cell);
				if (fine_int != false)
				{
					for (double Z = -250.0; Z < 250.0; Z = Z + dZ)
					//for (double Z = 0.0; Z < 0.00001; Z = Z + dZ)
					{
						fine_int = SS.Get_param(X, Y, Z + dZ / 2.0, parameters, prev_cell, next_cell);
						//fine_int = SS.Get_param(X, Z + dZ / 2.0, Y, parameters, prev_cell, next_cell);
						//fine_int = SS.Get_param(Z + dZ / 2.0, X, Y, parameters, prev_cell, next_cell);
						if (fine_int == false) continue;
						for (short int i = 0; i < 6; i++) next_cell[i] = prev_cell[i];

						ne = parameters["rho"];
						T = parameters["p"] / parameters["rho"];
						IH += kv(ne) * this->phys_param->interpolate_alpha_eff_Ha(T * ch_T) * dZ;
						this->phys_param->interpolate_Xray(T * ch_T, a1, a2);
						I_soft_X_ray += kv(ne) * a1 * dZ;
						I_hard_X_ray += kv(ne) * a2 * dZ;

						//cout << ne << " " << this->phys_param->interpolate_alpha_eff_Ha(T * ch_T) << " " << T << " " << T * ch_T << endl;
					}
				}

				//#pragma omp critical 
				//{

				fout << X * ch_r << " " << Y * ch_r << " " << IH * ch_Halpha << endl;  // 2.91892 * 1E10
				fout2 << X * ch_r << " " << Y * ch_r << " " << I_soft_X_ray * 33.2504 << endl;
				fout3 << X * ch_r << " " << Y * ch_r << " " << I_hard_X_ray * 173.48 << endl;
				//}
			}

		}

		fout.close();
		fout2.close();
		fout3.close();
	}
	else if (alg == 25)
	{
		
		// Перевод в размерные единицы
		double ch_r = 0.0140668;          // расстояние для переведа в парсеки
		double ch_T = 39111.5;            // температура для перевода в Кельвины
		double ch_r_0 = 0.0000230386;     // безразмерный радиус звезды


		//Dust_spectra DDD = Dust_spectra();

		//Считаем интегралл Kabs(lambda) * F_lambda(lambda)
		if (true)
		{
			double S = 0.0;
			double k, f;
			double dlambda = 1E-9;
			double lambda = 1E-9;
			int ki = 0;
			while (lambda < 0.01)
			{
				ki++;
				k = DDD->interpolate_K_abs(lambda + dlambda / 2.0);
				f = DDD->interpolate_F_kurucz(lambda + dlambda / 2.0);
				S = S + k * f * dlambda;
				lambda = lambda + dlambda;
				if (ki > 500000)
				{
					//cout << "1 proverca:  lambda = " << lambda << "   S = " << S << endl;
					ki = 0;
				}
			}

			DDD->Int1 = S;
			cout << "DDD->Int1 = " << DDD->Int1 << endl;
		}

		ofstream fout;
		fout.open("Kabs(lambda)_planck_b_lambda.txt");

		//Считаем интегралл Kabs(lambda) * planck_b_lambda
		int NN = DDD->L_em_NN;
		double TL = DDD->L_em_TL;
		double TR = DDD->L_em_TR;
		if (true)
		{
			DDD->L_em.resize(NN);
			for (int i = 0; i < NN; i++)
			{
				DDD->L_em[i] = 0.0;
			}

			double k, f;
			double dlambda = 1E-6;
			double lambda = 1E-8;
			int ki = 0;
			while (lambda < 5.0)
			{
				k = DDD->interpolate_K_abs(lambda + dlambda / 2.0);
				for (int i = 0; i < NN; i++)
				{
					double T = TL + (i + 0.5) * (TR - TL) / NN;
					f = DDD->planck_b_lambda(lambda, T);
					DDD->L_em[i] += k * f * dlambda;
				}
				lambda = lambda + dlambda;
			}
		}

		for (int i = 0; i < NN; i++)
		{
			double T = TL + (i + 0.5) * (TR - TL) / NN;
			fout << T << " " << DDD->L_em[i] << endl;
		}
		fout.close();
		DDD->buildInverseTable(500); 


		// Вычисляем температуру пыли в каждой точке
		for (auto& A : this->All_Cell)
		{
			double r = norm2(A->center[0][0], A->center[0][1], A->center[0][2]);
			double aa = kv(ch_r_0 / r) * DDD->Int1 / (4.0 * const_pi);
			double Tdust = DDD->getTemperatureFromL(aa);
			A->parameters[0]["Tdust"] = Tdust;
			A->parameters[1]["Tdust"] = Tdust;

			if (A->parameters[0]["p"] / A->parameters[0]["rho"] * ch_T > 100000)
			{
				A->parameters[0]["rhodust"] = 0.0;
				A->parameters[1]["rhodust"] = 0.0;
				A->parameters[0]["Tdust"] = 0.0;
				A->parameters[1]["Tdust"] = 0.0;
			}
			else
			{
				A->parameters[0]["rhodust"] = A->parameters[0]["rho"] / 165.0;
				A->parameters[1]["rhodust"] = A->parameters[0]["rho"] / 165.0;
			}
		}
		this->phys_param->param_names.push_back("Tdust");
		this->phys_param->param_names.push_back("rhodust");

		// Надо сохранить температуру и плотность пыли в файл

		//std::string filename = "dust_paremeter_0.bin";
		std::string filename = "2.12-2.dust_parameter.bin";
		std::ofstream file(filename, std::ios::binary);
		if (!file.is_open()) 
		{
			std::cerr << "Error gjiuerhgyh7845ygfudhger " << filename << std::endl;
			exit(-1);
		}

		for (auto& A : this->All_Cell)
		{
			double a1 = A->parameters[0]["rhodust"];
			double a2 = A->parameters[0]["Tdust"];
			file.write(reinterpret_cast<const char*>(&a1), sizeof(a1));
			file.write(reinterpret_cast<const char*>(&a2), sizeof(a2));
		}

		// Надо считать температуру пыли из файла
		if (false)
		{
			std::string filename = "dust_paremeter_1-rho10.bin";
			std::ifstream file(filename, std::ios::binary);
			if (!file.is_open())
			{
				std::cerr << "Error gjiuerhgyh7845ygfudhger " << filename << std::endl;
				exit(-1);
			}

			for (auto& A : this->All_Cell)
			{
				double a1, a2;

				file.read(reinterpret_cast<char*>(&a1), sizeof(a1));
				file.read(reinterpret_cast<char*>(&a2), sizeof(a2));
				A->parameters[0]["Tdust"] = A->parameters[1]["Tdust"] = a2;
			}
		}

		
	}
	else if (alg == 26)
	{

		// Перевод в размерные единицы
		double ch_r = 0.0140668;          // расстояние для переведа в парсеки
		double ch_T = 39111.5;            // температура для перевода в Кельвины
		double ch_r_0 = 0.0000230386;     // безразмерный радиус звезды
		double ch_rho_r = 4.37058E-7;     // характерная плотность * размер

		// Надо считать температуру пыли из файла
		if (false)
		{
			std::string filename = "2.8-2.dust_parameter.bin";
			std::ifstream file(filename, std::ios::binary);
			if (!file.is_open())
			{
				std::cerr << "Error gjiuerhgyh7845ygfudhger " << filename << std::endl;
				exit(-1);
			}

			for (auto& A : this->All_Cell)
			{
				double a1, a2;

				file.read(reinterpret_cast<char*>(&a1), sizeof(a1));
				file.read(reinterpret_cast<char*>(&a2), sizeof(a2));
				A->parameters[0]["Tdust"] = A->parameters[1]["Tdust"] = a2;
			}
		}


		// Теперь рисуем сами карты
		this->Save_for_interpolate("For_intertpolate_work3.bin", false);
		Interpol SS = Interpol("For_intertpolate_work3.bin");


		double dX = 2.0;  // минимум по 0.5, но лучше меньше
		double dY = 2.0;  // минимум по 0.5, но лучше меньше

		//double dX = 1.0;  // минимум по 0.5, но лучше меньше
		//double dY = 1.0;  // минимум по 0.5, но лучше меньше

		//double dX = 0.25;  // минимум по 0.5, но лучше меньше
		//double dY = 0.25;  // минимум по 0.5, но лучше меньше

		double XL = -285.0;
		double XR = 100.0;
		double YL = -213.0;
		double YR = 213.0;

		/*double XL = -71.0;
		double XR = 71.0;
		double YL = -142.0;
		double YR = 142.0;*/

		/*double XL = 0.0;
		double XR = 57.0;
		double YL = -71.0;
		double YR = 71.0;*/


		double dZ = 0.05;
		double lambda_0 = 24E-4;
		double kk_abs = DDD->interpolate_K_abs(lambda_0);
		double kk_sca = DDD->interpolate_K_sca(lambda_0);
		double kk_ext = kk_sca + kk_abs;
		double lambda_1 = 70E-4;
		double kk_abs2 = DDD->interpolate_K_abs(lambda_1);
		double kk_sca2 = DDD->interpolate_K_sca(lambda_1);
		double kk_ext2 = kk_sca2 + kk_abs2;

		
		ofstream fout;
		fout.open("2.12-2-infrared-v0.txt");
		// -0, -1, -2 - три разрешения в порядке уменьтшения размера между точками


		const int NX = static_cast<int>(65.0 / dX + 0.5); // ~130


		//for (double X = 185.0; X > -285.0; X = X - dX)
		//for (double X = 120.0; X > -100.0; X = X - dX)
		for (double X = XR; X > XL; X = X - dX)
			//for (int iX = 0; iX < NX; ++iX)
		{
			cout << "Culk for X = " << X << endl;

			double rhodust, Tdust;

			//for (double Y = -250.0; Y < 250.0; Y = Y + dY)
			//for (double Y = -225.0; Y < 225.0; Y = Y + dY)
			//for (double Y = -65.0; Y < 65.0; Y = Y + dY)
			for (double Y = YL; Y < YR; Y = Y + dY)
			{
				std::array<Cell_handle, 6> prev_cell;
				std::array<Cell_handle, 6> next_cell;
				for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();
				std::unordered_map<string, double> parameters;
				bool fine_int;

				double IH = 0.0;
				double IH2 = 0.0;
				double qq = 0.0;
				double qq2 = 0.0;
				double a1, a2;
				fine_int = SS.Get_param(X, Y, 0.0, parameters, prev_cell, next_cell);

				if (fine_int != false)
				{
					for (double Z = 250.0; Z > -250.0; Z = Z - dZ)
					{
						if (norm2(X, Y, Z + dZ / 2.0) < 17.0) continue;
						fine_int = SS.Get_param(X, Y, Z + dZ / 2.0, parameters, prev_cell, next_cell);
						if (fine_int == false) continue;
						for (short int i = 0; i < 6; i++) next_cell[i] = prev_cell[i];

						rhodust = parameters["rhodust"];
						Tdust = parameters["Tdust"];

						qq += rhodust * kk_ext * dZ * ch_rho_r;
						qq2 += rhodust * kk_ext2 * dZ * ch_rho_r;
						//cout << "qq = " << qq << "  " << rhodust * kk_ext * dZ * 3.35078E-7 << endl;
						IH += rhodust * kk_abs * DDD->planck_b_lambda(lambda_0, Tdust) * dZ * exp(-qq);
						IH2 += rhodust * kk_abs2 * DDD->planck_b_lambda(lambda_1, Tdust) * dZ * exp(-qq2);
					}
				}


				fout << X * ch_r << " " << Y * ch_r << " " << IH * ch_rho_r << " " << IH2 * ch_rho_r << endl;  // 2.91892 * 1E10
			}

		}

		fout.close();
	}
	else if (alg == 27)
	{
		this->MK_go_dust();
	}
	
		

	cout << "End Algoritm " << alg << endl;
}
