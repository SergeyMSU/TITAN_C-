#pragma once
#include "Header.h"
using namespace std;


class Dust_spectra
{
public:

	double Int1 = 0.0;  // интегралл Kabs(lambda) * F_lambda(lambda)  dlambda
	double E_esc = 0.0;   // Суммарная энергия ушедших из области пакетов
	mutex mut;

	std::vector<double> L_em;   // Вектор по температуре - интегралл от Kabs * Bplanka
	int L_em_NN = 500;
	double L_em_TL = 1.0;
	double L_em_TR = 1000.0;

	// Обратная таблица: T(L_em) с равномерным шагом по ln(L_em)
	std::vector<double> inv_logL;   // узлы ln(L_em)
	std::vector<double> inv_logT;   // узлы ln(T)
	int inv_N = 0;                  // размер обратной таблицы
	double inv_logL_min = 0.0;      // минимальный ln(L_em)
	double inv_logL_max = 0.0;      // максимальный ln(L_em)
	double inv_logL_step = 0.0;     // шаг по ln(L_em)

	// Спектр звезды
	std::vector<double> lambda_F_kurucz;      // исходные длины волн (см)
	std::vector<double> value_F_kurucz;       // исходные значения K_abs
	double Itot = 0.0;                        // Интегральный поток на поверхности звезды

	// Для быстрой выборки (inverse CDF)
	std::vector<double> cum_int_F_kurucz;   // массив значений первообразной в узлах
	std::vector<double> inv_cdf_lambda_F;   // массив длин волн для равномерных вероятностей
	int inv_cdf_size = 200000;               // размер сетки (можно менять)


	// Свойства пыли
	std::vector<double> lambda_K_abs;      // исходные длины волн (см)
	std::vector<double> value_K_abs;       // исходные значения K_abs
	std::vector<double> lnLambda_K_abs;    // логарифмы длин волн
	std::vector<double> lnValue_K_abs;     // логарифмы значений


	std::vector<double> lambda_K_sca;      // исходные длины волн (см)
	std::vector<double> value_K_sca;       // исходные значения K_abs
	std::vector<double> lnLambda_K_sca;    // логарифмы длин волн
	std::vector<double> lnValue_K_sca;     // логарифмы значений


	std::vector<double> lambda_g;      // длины волн (см)
	std::vector<double> g_vec;         // значения коэффициента g
	std::vector<double> lnLambda_g;    // логарифмы длин волн для быстрого поиска

	// ----------------------------------------------------------------------------------------------------------------------
	// Далее всё для метода переизлучения пакеты из ячейки с заданной температурой
	// ---------- сетки и таблица (переименованы с _sca) ----------
	std::vector<double> T_grid_sca;              // температуры [K]
	std::vector<double> lambda_grid_sca;         // длины волн [см] (логарифмическая сетка)
	std::vector<std::vector<double>> cdf_table_sca; // CDF для каждой температуры
	std::vector<std::vector<double>> inv_lambda_table_sca;
	int N_prob_sca = 5000;  // число точек в равномерной сетке по вероятности

	// Подготовить таблицы CDF для диапазона температур и длин волн
	void prepare_sca(double T_min, double T_max, int N_temp,
		double lambda_min, double lambda_max, int M_lambda);

	// Сохранить все таблицы в файл
	void save_sca(const std::string& filename) const;
	// Загрузить таблицы из файла
	void load_sca(const std::string& filename);

	// Основная функция: возвращает длину волны в СГС по случайному числу u ? [0,1]
	// и температуре T_K [K]. Значение T_K может находиться между узлами сетки.
	double sample_frequency_sca(double uniform_rand, double T_K) const;

	// ---------- вспомогательные методы ----------
	void find_temp_interval_sca(double T_K, size_t& idx_low, size_t& idx_high, double& frac) const;
	void prepare_inverse_sca(int N_prob);



	// -----------------------------------------------------------------------------------------------------------------------
	Dust_spectra();


	void Read_F_kurucz();
	double interpolate_F_kurucz(double lambda);
	void integrate_F_kurucz();                     // Вычисляем интегральный поток на поверхности звезды
	void build_inverse_cdf_F_kurucz();
	double sample_F_kurucz(double xi) const;

	void Read_K_abs();
	double interpolate_K_abs(double lambda);

	void Read_K_sca();
	double interpolate_K_sca(double lambda);

	void Read_g();
	double interpolate_g(double lambda);

	double planck_b_lambda(double lambda_cm, double T_K);


	double temperatureAtIndex(int i) const 
	{
		// исходная сетка T равномерна от TL до TR с точками в центре интервалов
		// T_i = TL + (i + 0.5) * (TR - TL) / NN
		return L_em_TL + (i + 0.5) * (L_em_TR - L_em_TL) / L_em_NN;
	}
	
	// Построение обратной функцйии для L_em для быстрой интерполяции
	void buildInverseTable(int nPoints = 500);

	// Возвращает температуру T по заданному значению L = L_em(T)
	double getTemperatureFromL(double L) const;


};

