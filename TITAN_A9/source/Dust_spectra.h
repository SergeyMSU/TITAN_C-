#pragma once
#include "Header.h"
using namespace std;


class Dust_spectra
{
public:

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

};

