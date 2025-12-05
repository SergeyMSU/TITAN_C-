#include "Setka.h"
#include <algorithm>
#include <filesystem> // Для работы с файловой системой
namespace fs = std::filesystem; // Создаем псевдоним для удобства

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
	// 11 - расчёт поверхностных токов на разрывах
	// 12 - расчёт объёмных токов
	// 13 - просмотр источников S+/S- и сравнение их с флюидными источниками

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


		for (int i = 1; i <= 6 * 6; i++) // 6 * 2   12 * 5
		{
			auto start = std::chrono::high_resolution_clock::now();
			cout << "IIIII = " << i << endl;

			//S1.Go(true, 600, 1); // 400   1
			cout << "All time = " << this->phys_param->ALL_Time << endl;
			cout << "All time (in days) = " << this->phys_param->ALL_Time / 0.00142358 << endl;
			cout << "All time (in years) = " << this->phys_param->ALL_Time / 0.519607 << endl;
			this->Go(false, 400, 1); // 400   1
			if (i % 100000 == 0)
			{
				this->Go(true, 1000, 1); // 400   1 
			}
			else
			{
				this->Go(true, 100, 1); // 400   1 
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

			if (i % 12 == 0)
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
		// Из основной сетки создаём интерполяционную сетку
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
		if (true)
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
		if (Smc.phys_param->culc_cell_moments == true)
		{
			if (file_exists(Smc.phys_param->MK_file))
			{
				Smc.Download_cell_MK_parameters(Smc.phys_param->MK_file, 1000);
			}
		}

		cout << "2 Sum = " << A->mas_pogl.sum() << endl;

		cout << "Arrays read successfully" << endl;

		Smc.mas_pogl_Culc(1.0, 0.0, 0.0, "upwind");
		Smc.mas_pogl_Culc(1.0, 0.1, 0.0, "sim_upwind");
		Smc.mas_pogl_Culc(0.0, 1.0, 0.0, "crosswind1");
		Smc.mas_pogl_Culc(0.0, 1.0, 1.0, "crosswind2");
		Smc.mas_pogl_Culc(0.0, 0.0, 1.0, "crosswind3");
		Smc.mas_pogl_Culc(-1.0, 0.0, 0.0, "downwind");
		Smc.mas_pogl_Culc(-1.0, 1.0, 0.0, "tail1");
		Smc.mas_pogl_Culc(-1.0, 0.70710678, 0.70710678, "tail2");
		Smc.mas_pogl_Culc(-1.0, 0.0, 1.0, "tail3");

		cout << "Removing arrays" << endl;

		for (size_t idx = 0; idx < Smc.All_Cell.size(); ++idx)
		{
			auto A = Smc.All_Cell[idx];
			A->Delete_mas_pogl();
		}

		Smc.Print_f_proect_in_gran(1);
		Smc.Print_f_proect_in_gran(2);
		Smc.Print_f_proect_in_gran(3);
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
						to_string(5) + "_" + to_string(gr->number) + ".bin";

					if (std::filesystem::exists(name_f))
					{
						std::filesystem::remove(name_f);
					}
				}
			}

			for (auto& gr : Smc.All_Gran)
			{
				for (int ii = 0; ii <= 1; ii++)
				{
					string name_f = Smc.phys_param->AMR_folder + "/" + "func_grans_AMR_" + to_string(ii) + "_H" +
						to_string(6) + "_" + to_string(gr->number) + ".bin";

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
		//zones_number.push_back(4); zones_n_koeff.push_back(1.0);

		zones_number.push_back(2); zones_n_koeff.push_back(1.0);
		zones_number.push_back(2); zones_n_koeff.push_back(1.0);
		


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
	else if (alg == 11)
	{
		ofstream fout;
		string name_f;

		//HP
		if (true)
		{
			name_f = "HP_J.txt";
			fout.open(name_f);
			fout << "TITLE = HP  VARIABLES = x, y, z, phi, the, Jx, Jy, Jz, |J|, J2x, J2y, J2z, |J2|, Bx_L, By_L, Bz_L, Bx_R, By_R, Bz_R" << endl;
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


				B1[0] = A->parameters[0]["Bx"];
				B1[1] = A->parameters[0]["By"];
				B1[2] = A->parameters[0]["Bz"];

				B2[0] = B->parameters[0]["Bx"];
				B2[1] = B->parameters[0]["By"];
				B2[2] = B->parameters[0]["Bz"];


				BB1 = B2 - B1;
				BB2 = B2 + B1;

				//cout << "do = " << B1[0] << endl;
				J = n.cross(BB1);
				J2 = n.cross(BB2);
				//cout << "posle = " << B1[0] << endl;

				J = J / (4.0 * const_pi);
				J2 = J2 / (4.0 * const_pi);

				for (auto& j : i->yzels)
				{
					cc[0] = j->coord[0][0];
					cc[1] = j->coord[0][1];
					cc[2] = j->coord[0][2];

					fout << cc[0] << " " << cc[1] << " " << cc[2] << " " <<
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
		if (true)
		{
			name_f = "TS_J.txt";
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

				J = J / (4.0 * const_pi);

				for (auto& j : i->yzels)
				{
					cc[0] = j->coord[0][0];
					cc[1] = j->coord[0][1];
					cc[2] = j->coord[0][2];

					fout << cc[0] << " " << cc[1] << " " << cc[2] << " " <<
						polar_angle(cc[1], cc[2]) << " " << polar_angle(cc[0], norm2(0.0, cc[1], cc[2])) << " " <<
						J[0] << " " << J[1] << " " << J[2] << " " << J.norm() << endl;
				}
			}


			for (int k = 0; k < this->Gran_TS.size(); k++)
			{
				fout << 4 * k + 1 << " " << 4 * k + 2 << " " << 4 * k + 3 << " " << 4 * k + 4 << endl;
			}


			fout.close();
		}

		// BS
		if (true)
		{
			name_f = "BS_J.txt";
			fout.open(name_f);
			fout << "TITLE = HP  VARIABLES = x, y, z, phi, the, Jx, Jy, Jz, |J|" << endl;
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

				J = J / (4.0 * const_pi);

				for (auto& j : i->yzels)
				{
					cc[0] = j->coord[0][0];
					cc[1] = j->coord[0][1];
					cc[2] = j->coord[0][2];

					fout << cc[0] << " " << cc[1] << " " << cc[2] << " " <<
						polar_angle(cc[1], cc[2]) << " " << polar_angle(cc[0], norm2(0.0, cc[1], cc[2])) << " " <<
						J[0] << " " << J[1] << " " << J[2] << " " << J.norm() << endl;
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
		this->Culc_usual_rotors_in_cell();


		// Надо улучшить ротеры вблизи разрывов
		for (auto& gr : this->Gran_TS)
		{
			auto C1 = gr->cells[1];
			auto C2 = gr->cells_TVD[1];

			C1->parameters[0]["rotB_x"] = C2->parameters[0]["rotB_x"];
			C1->parameters[0]["rotB_y"] = C2->parameters[0]["rotB_y"];
			C1->parameters[0]["rotB_z"] = C2->parameters[0]["rotB_z"];
		}

		for (auto& gr : this->Gran_HP)
		{
			auto C1 = gr->cells[0];
			auto C2 = gr->cells_TVD[0];

			C1->parameters[0]["rotB_x"] = C2->parameters[0]["rotB_x"];
			C1->parameters[0]["rotB_y"] = C2->parameters[0]["rotB_y"];
			C1->parameters[0]["rotB_z"] = C2->parameters[0]["rotB_z"];
		}



		// Рисует тетраэдры в текплот
		if (true)
		{


			this->Save_for_interpolate_one_zone_only("For_intertpolate_work.bin", Type_cell::Zone_2);
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
					params->parameters["rotB_x"] << " " << params->parameters["rotB_y"] << " " << params->parameters["rotB_z"] << " " <<
					kvv(params->parameters["rotB_x"], params->parameters["rotB_y"], params->parameters["rotB_z"]) << endl;
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

		// Трассируем линии тока
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
			for (auto& gr : this->Gran_TS)
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
					double nnn = norm2(parameters["rotB_x"], parameters["rotB_y"], parameters["rotB_z"]);
					line.push_back({ x, y, z, nnn });

					x = x + 0.02 * v * parameters["rotB_x"] / nnn;
					y = y + 0.02 * v * parameters["rotB_y"] / nnn;
					z = z + 0.02 * v * parameters["rotB_z"] / nnn;

					CC = nullptr;
					CC = Find_cell_point(x, y, z, 0, prev);
					if (CC == nullptr) break;
					if (CC->type != Type_cell::Zone_2) break;

					bb = SS.Get_param(x, y, z, parameters);
					if (bb == false) break;
				}

#pragma omp critical (dsds2) 
				{
					file << "ZONE T=\"Line" << line_count++ << "\" I=" << line.size() << " F=POINT" << std::endl;
					for (const auto& point : line)
					{
						file << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << std::endl;
					}
				}

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

	cout << "End Algoritm " << alg << endl;
}
