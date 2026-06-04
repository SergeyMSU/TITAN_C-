#include "Dust_spectra.h"


Dust_spectra::Dust_spectra()
{
    this->Read_K_abs();
    this->Read_K_sca();
    this->Read_g();
    this->Read_F_kurucz();
    this->integrate_F_kurucz();
    this->build_inverse_cdf_F_kurucz();

    if (false)
    {
        cout << "Test: F_kurucz sample" << endl;
        cout << 0.1 << "  " << sample_F_kurucz(0.1) << endl;
        cout << 0.2 << "  " << sample_F_kurucz(0.2) << endl;
        cout << 0.3 << "  " << sample_F_kurucz(0.3) << endl;
        cout << 0.4 << "  " << sample_F_kurucz(0.4) << endl;
        cout << 0.5 << "  " << sample_F_kurucz(0.5) << endl;
        cout << 0.6 << "  " << sample_F_kurucz(0.6) << endl;
        cout << 0.7 << "  " << sample_F_kurucz(0.7) << endl;
        cout << 0.8 << "  " << sample_F_kurucz(0.8) << endl;
        cout << 0.9 << "  " << sample_F_kurucz(0.9) << endl;
    }
}

void Dust_spectra::Read_F_kurucz()
{
    //std::string filename = "solar_spectrum_kurucz_cgs.txt";
    std::string filename = "star_spectrum_kurucz_cgs.txt";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        exit(-1);
    }

    std::string line;
    // Пропускаем заголовок (если есть)
    std::getline(file, line);

    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        std::istringstream iss(line);
        double lambda, value;
        if (!(iss >> lambda >> value))
        {
            std::cerr << "Warning: skipped invalid line: " << line << std::endl;
            continue;
        }
        lambda_F_kurucz.push_back(lambda);
        value_F_kurucz.push_back(value);
    }
    file.close();

    if (lambda_F_kurucz.empty())
    {
        std::cerr << "Error: no data read from file." << std::endl;
        exit(-1);
    }

    // Проверка монотонности
    if (!std::is_sorted(lambda_F_kurucz.begin(), lambda_F_kurucz.end()))
    {
        std::cerr << "Warning: lambda_F_kurucz are not in increasing order." << std::endl;
    }

    // ------------------- ПРОВЕРОЧНАЯ ПЕЧАТЬ -------------------
    std::ofstream check_file("check_F_kurucz_interp.txt");
    if (!check_file.is_open())
    {
        std::cerr << "Warning: could not create check file." << std::endl;
        return;
    }

    double lam_min = lambda_F_kurucz.front();
    double lam_max = lambda_F_kurucz.back();
    int N_steps = 200000;
    double step = (lam_max - lam_min) / (N_steps - 1);

    check_file << "# lambda_cm    F_kurucz_interpolated\n";
    for (int i = 0; i < N_steps; ++i)
    {
        double lam = lam_min + i * step;
        double f_interp = interpolate_F_kurucz(lam);
        check_file << std::scientific << std::setprecision(8)
            << lam << "  " << f_interp << "\n";
    }
    check_file.close();

    std::cout << "Check file written: check_F_kurucz_interp.txt" << std::endl;
}

double Dust_spectra::interpolate_F_kurucz(double lambda)
{
    // Граничные условия
    if (lambda <= lambda_F_kurucz.front())
        return 0.0;
    if (lambda >= lambda_F_kurucz.back())
        return 0.0;

    // Поиск интервала, содержащего lambda
    auto it = std::lower_bound(lambda_F_kurucz.begin(), lambda_F_kurucz.end(), lambda);
    size_t idx = it - lambda_F_kurucz.begin();

    // Точное совпадение
    if (std::fabs(lambda_F_kurucz[idx] - lambda) < 1e-12)
        return value_F_kurucz[idx];

    // Линейная интерполяция
    const double& x0 = lambda_F_kurucz[idx - 1];
    const double& x1 = lambda_F_kurucz[idx];
    const double& y0 = value_F_kurucz[idx - 1];
    const double& y1 = value_F_kurucz[idx];

    double t = (lambda - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}

void Dust_spectra::Read_K_abs()
{
    std::string filename = "lambda_K_abs.txt";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        exit(-1);
    }

    std::string line;
    // Пропускаем заголовок
    std::getline(file, line);

    // Чтение данных
    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        std::istringstream iss(line);
        double lambda, value;
        if (!(iss >> lambda >> value))
        {
            std::cerr << "Warning: skipped invalid line: " << line << std::endl;
            continue;
        }
        lambda_K_abs.push_back(lambda);
        value_K_abs.push_back(value);
    }
    file.close();

    if (lambda_K_abs.empty())
    {
        std::cerr << "Error: no data read from file." << std::endl;
        exit(-1);
    }

    // Предварительно вычисляем логарифмы
    lnLambda_K_abs.reserve(lambda_K_abs.size());
    lnValue_K_abs.reserve(value_K_abs.size());
    for (size_t i = 0; i < lambda_K_abs.size(); ++i)
    {
        lnLambda_K_abs.push_back(std::log(lambda_K_abs[i]));
        lnValue_K_abs.push_back(std::log(value_K_abs[i]));
    }

    // Проверка монотонности (по логарифмам)
    if (!std::is_sorted(lnLambda_K_abs.begin(), lnLambda_K_abs.end()))
    {
        std::cerr << "Warning: lambda values are not in increasing order." << std::endl;
    }

    // --- ПЕЧАТЬ ДЛЯ ПРОВЕРКИ ---
    std::ofstream check_file("check_K_abs_interp.txt");
    if (!check_file.is_open())
    {
        std::cerr << "Warning: could not create check file." << std::endl;
        return;
    }

    // Логарифмически равномерный шаг
    double lam_min = lambda_K_abs.front();
    double lam_max = lambda_K_abs.back();
    double log_lam_min = std::log(lam_min);
    double log_lam_max = std::log(lam_max);
    int N_steps = 10000;
    double log_step = (log_lam_max - log_lam_min) / (N_steps - 1);

    check_file << "# lambda_cm    K_abs_interpolated\n";
    for (int i = 0; i < N_steps; ++i)
    {
        double lam = std::exp(log_lam_min + i * log_step);
        double k_interp = interpolate_K_abs(lam);
        check_file << std::scientific << std::setprecision(8)
            << lam << "  " << k_interp << "\n";
    }

    std::cout << "Check file written: check_K_abs_interp.txt" << std::endl;
}

double Dust_spectra::interpolate_K_abs(double lambda)
{
    // Если лямбда меньше минимальной – экстраполяция по первым двум точкам
    if (lambda <= lambda_K_abs.front())
    {
        return 0.0;
    }

    // Если лямбда больше максимальной – экстраполяция по последним двум точкам
    if (lambda >= lambda_K_abs.back())
    {
        return 0.0;
    }

    // Поиск интервала, содержащего lambda (в исходном, не логарифмическом масштабе)
    // Но поскольку данные монотонны, можно использовать логарифмический поиск
    double lnLambda = std::log(lambda);
    auto it = std::lower_bound(lnLambda_K_abs.begin(), lnLambda_K_abs.end(), lnLambda);
    size_t idx = it - lnLambda_K_abs.begin();

    // Точное совпадение с узлом
    if (std::fabs(*it - lnLambda) < 1e-12)
    {
        return value_K_abs[idx];
    }

    // Линейная интерполяция в логарифмическом масштабе
    const double& lnL0 = lnLambda_K_abs[idx - 1];
    const double& lnL1 = lnLambda_K_abs[idx];
    const double& lnV0 = lnValue_K_abs[idx - 1];
    const double& lnV1 = lnValue_K_abs[idx];

    double lnVal = lnV0 + (lnLambda - lnL0) * (lnV1 - lnV0) / (lnL1 - lnL0);
    return std::exp(lnVal);
}

void Dust_spectra::Read_K_sca()
{
    std::string filename = "lambda_K_sca.txt";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        exit(-1);
    }

    std::string line;
    // Пропускаем заголовок
    std::getline(file, line);

    // Чтение данных
    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        std::istringstream iss(line);
        double lambda, value;
        if (!(iss >> lambda >> value))
        {
            std::cerr << "Warning: skipped invalid line: " << line << std::endl;
            continue;
        }
        lambda_K_sca.push_back(lambda);
        value_K_sca.push_back(value);
    }
    file.close();

    if (lambda_K_sca.empty())
    {
        std::cerr << "Error: no data read from file." << std::endl;
        exit(-1);
    }

    // Предварительно вычисляем логарифмы
    lnLambda_K_sca.reserve(lambda_K_sca.size());
    lnValue_K_sca.reserve(value_K_sca.size());
    for (size_t i = 0; i < lambda_K_sca.size(); ++i)
    {
        lnLambda_K_sca.push_back(std::log(lambda_K_sca[i]));
        lnValue_K_sca.push_back(std::log(value_K_sca[i]));
    }

    // Проверка монотонности (по логарифмам)
    if (!std::is_sorted(lnLambda_K_sca.begin(), lnLambda_K_sca.end()))
    {
        std::cerr << "Warning: lambda values are not in increasing order." << std::endl;
    }

    // --- ПЕЧАТЬ ДЛЯ ПРОВЕРКИ ---
    std::ofstream check_file("check_K_sca_interp.txt");
    if (!check_file.is_open())
    {
        std::cerr << "Warning: could not create check file." << std::endl;
        return;
    }

    // Логарифмически равномерный шаг
    double lam_min = lambda_K_sca.front();
    double lam_max = lambda_K_sca.back();
    double log_lam_min = std::log(lam_min);
    double log_lam_max = std::log(lam_max);
    int N_steps = 10000;
    double log_step = (log_lam_max - log_lam_min) / (N_steps - 1);

    check_file << "# lambda_cm    K_sca_interpolated\n";
    for (int i = 0; i < N_steps; ++i)
    {
        double lam = std::exp(log_lam_min + i * log_step);
        double k_interp = interpolate_K_sca(lam);
        check_file << std::scientific << std::setprecision(8)
            << lam << "  " << k_interp << "\n";
    }

    std::cout << "Check file written: check_K_sca_interp.txt" << std::endl;
}

double Dust_spectra::interpolate_K_sca(double lambda)
{
    // Если лямбда меньше минимальной – экстраполяция по первым двум точкам
    if (lambda <= lambda_K_sca.front())
    {
        return 0.0;
    }

    // Если лямбда больше максимальной – экстраполяция по последним двум точкам
    if (lambda >= lambda_K_sca.back())
    {
        return 0.0;
    }

    // Поиск интервала, содержащего lambda (в исходном, не логарифмическом масштабе)
    // Но поскольку данные монотонны, можно использовать логарифмический поиск
    double lnLambda = std::log(lambda);
    auto it = std::lower_bound(lnLambda_K_sca.begin(), lnLambda_K_sca.end(), lnLambda);
    size_t idx = it - lnLambda_K_sca.begin();

    // Точное совпадение с узлом
    if (std::fabs(*it - lnLambda) < 1e-12)
    {
        return value_K_sca[idx];
    }

    // Линейная интерполяция в логарифмическом масштабе
    const double& lnL0 = lnLambda_K_sca[idx - 1];
    const double& lnL1 = lnLambda_K_sca[idx];
    const double& lnV0 = lnValue_K_sca[idx - 1];
    const double& lnV1 = lnValue_K_sca[idx];

    double lnVal = lnV0 + (lnLambda - lnL0) * (lnV1 - lnV0) / (lnL1 - lnL0);
    return std::exp(lnVal);
}

void Dust_spectra::Read_g()
{
    std::string filename = "lambda_g.txt";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        exit(-1);
    }

    std::string line;
    // Пропускаем заголовок (первая строка)
    std::getline(file, line);

    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        std::istringstream iss(line);
        double lambda, g;
        if (!(iss >> lambda >> g))
        {
            std::cerr << "Warning: skipped invalid line: " << line << std::endl;
            continue;
        }
        lambda_g.push_back(lambda);
        g_vec.push_back(g);
    }
    file.close();

    if (lambda_g.empty())
    {
        std::cerr << "Error: no data read from file." << std::endl;
        exit(-1);
    }

    // Вычисляем логарифмы длин волн для интерполяции
    lnLambda_g.reserve(lambda_g.size());
    for (size_t i = 0; i < lambda_g.size(); ++i)
        lnLambda_g.push_back(std::log(lambda_g[i]));

    if (!std::is_sorted(lnLambda_g.begin(), lnLambda_g.end()))
        std::cerr << "Warning: lambda_g are not in increasing order." << std::endl;

    // ------------------- ПРОВЕРОЧНАЯ ПЕЧАТЬ -------------------
    std::ofstream check_file("check_g_interp.txt");
    if (!check_file.is_open())
    {
        std::cerr << "Warning: could not create check file." << std::endl;
        return;
    }

    double lam_min = lambda_g.front();
    double lam_max = lambda_g.back();
    double log_lam_min = std::log(lam_min);
    double log_lam_max = std::log(lam_max);
    int N_steps = 10000;
    double log_step = (log_lam_max - log_lam_min) / (N_steps - 1);

    check_file << "# lambda_cm    g_interpolated\n";
    for (int i = 0; i < N_steps; ++i)
    {
        double lam = std::exp(log_lam_min + i * log_step);
        double g_interp = interpolate_g(lam);
        check_file << std::scientific << std::setprecision(8)
            << lam << "  " << g_interp << "\n";
    }
    check_file.close();

    std::cout << "Check file written: check_g_interp.txt" << std::endl;
}

double Dust_spectra::interpolate_g(double lambda)
{
    // Левая граница
    if (lambda <= lambda_g.front())
        return 1.0;

    // Правая граница
    if (lambda >= lambda_g.back())
        return 0.0;

    // Логарифм запрашиваемой длины волны
    double lnLambda = std::log(lambda);

    // Поиск интервала в lnLambda_g
    auto it = std::lower_bound(lnLambda_g.begin(), lnLambda_g.end(), lnLambda);
    size_t idx = it - lnLambda_g.begin();

    // Точное совпадение с узлом
    if (std::fabs(*it - lnLambda) < 1e-12)
        return g_vec[idx];

    // Линейная интерполяция g (не логарифмируя) по lnLambda
    const double& lnL0 = lnLambda_g[idx - 1];
    const double& lnL1 = lnLambda_g[idx];
    const double& g0 = g_vec[idx - 1];
    const double& g1 = g_vec[idx];

    double t = (lnLambda - lnL0) / (lnL1 - lnL0);
    return g0 + t * (g1 - g0);
}

void Dust_spectra::integrate_F_kurucz()
{
    if (lambda_F_kurucz.size() < 2)
    {
        cout << "Error iewjrfgiuegrfheirf34" << endl;
        exit(-1);
    }
        
    cum_int_F_kurucz.clear();
    cum_int_F_kurucz.reserve(lambda_F_kurucz.size());
    // Первообразная в начальной точке = 0
    cum_int_F_kurucz.push_back(0.0);

    double integral = 0.0;
    for (size_t i = 1; i < lambda_F_kurucz.size(); ++i)
    {
        double dx = lambda_F_kurucz[i] - lambda_F_kurucz[i - 1];
        double avg_value = (value_F_kurucz[i - 1] + value_F_kurucz[i]) * 0.5;
        integral += avg_value * dx;
        double inc = avg_value * dx;
        cum_int_F_kurucz.push_back(cum_int_F_kurucz.back() + inc);
    }

    this->Itot = integral;
}

void Dust_spectra::build_inverse_cdf_F_kurucz()
{
    inv_cdf_lambda_F.resize(this->inv_cdf_size);
    double total = this->Itot;
    int N = this->inv_cdf_size;

    for (int i = 0; i < N; ++i)
    {
        double prob = static_cast<double>(i) / (N - 1);   // 0..1
        double target = prob * total;

        // Найти интервал в cum_int_F_kurucz, содержащий target
        auto it = std::lower_bound(cum_int_F_kurucz.begin(), cum_int_F_kurucz.end(), target);
        size_t idx = it - cum_int_F_kurucz.begin();

        if (idx == 0)
        {
            inv_cdf_lambda_F[i] = lambda_F_kurucz[0];
            continue;
        }
        if (idx >= cum_int_F_kurucz.size())
        {
            inv_cdf_lambda_F[i] = lambda_F_kurucz.back();
            continue;
        }

        // Линейная интерполяция по кумулянте
        const double& cum0 = cum_int_F_kurucz[idx - 1];
        const double& cum1 = cum_int_F_kurucz[idx];
        const double& lam0 = lambda_F_kurucz[idx - 1];
        const double& lam1 = lambda_F_kurucz[idx];

        double t = (target - cum0) / (cum1 - cum0);
        inv_cdf_lambda_F[i] = lam0 + t * (lam1 - lam0);
    }
}

double Dust_spectra::sample_F_kurucz(double xi) const
{
    // xi ? [0,1]
    if (xi <= 0.0) return lambda_F_kurucz.front();
    if (xi >= 1.0) return lambda_F_kurucz.back();

    double idx_real = xi * (inv_cdf_size - 1);
    int idx = static_cast<int>(idx_real);
    double frac = idx_real - idx;

    if (frac < 1e-12) return inv_cdf_lambda_F[idx];
    if (idx + 1 >= inv_cdf_size) return inv_cdf_lambda_F.back();

    // Линейная интерполяция между соседними точками обратной таблицы
    return inv_cdf_lambda_F[idx] + frac * (inv_cdf_lambda_F[idx + 1] - inv_cdf_lambda_F[idx]);
}

double Dust_spectra::planck_b_lambda(double lambda_cm, double T_K)
{

    // Защита от нефизичных входных данных
    if (lambda_cm <= 0.0 || T_K <= 0.0) 
    {
        return 0.0;
    }

    // Предварительные вычисления
    const double two_h_c2 = 0.0000119106;          // 2·h·c?

    const double lambda5 = lambda_cm * lambda_cm *
        lambda_cm * lambda_cm *
        lambda_cm;                // ??
    const double exponent_arg = 1.43948 / (lambda_cm * T_K);

    // Предотвращение переполнения exp при больших аргументах
    if (exponent_arg > 700.0) 
    {
        return 0.0;   // интенсивность пренебрежимо мала
    }

    const double exp_val = std::exp(exponent_arg);
    const double denominator = exp_val - 1.0;

    // Обработка случая, когда знаменатель почти нулевой (малые аргументы)
    // Используем разложение 1/(exp(x)-1) ? 1/x при x ? 0
    if (std::abs(denominator) < 1e-12) 
    {
        // При очень малых x: B_? ? (2·h·c? / ??) * (k_B·T·? / (h·c))
        return two_h_c2 / lambda5 * (T_K * lambda_cm) / (1.43948);
    }

    return two_h_c2 / lambda5 / denominator;
}

// Построение обратной таблицы с равномерным шагом по ln(L_em)
void Dust_spectra::buildInverseTable(int nPoints)
{
    if (L_em.empty() || L_em_NN < 2) {
        throw std::runtime_error("Исходная таблица L_em пуста или слишком мала");
    }

    // Определяем диапазон ln(L_em)
    double logL_min = std::log(L_em.front());
    double logL_max = std::log(L_em.back());
    inv_N = nPoints;
    inv_logL_min = logL_min;
    inv_logL_max = logL_max;
    inv_logL_step = (logL_max - logL_min) / (inv_N - 1);

    inv_logL.clear();
    inv_logT.clear();
    inv_logL.reserve(inv_N);
    inv_logT.reserve(inv_N);

    // Подготовим вспомогательные векторы логарифмов исходной таблицы
    std::vector<double> logL_em(L_em_NN);
    std::vector<double> logT(L_em_NN);
    for (int i = 0; i < L_em_NN; ++i) {
        logL_em[i] = std::log(L_em[i]);
        logT[i] = std::log(temperatureAtIndex(i));
    }

    // Для каждой точки обратной сетки выполняем интерполяцию в исходной таблице
    for (int i = 0; i < inv_N; ++i) {
        double logL_target = inv_logL_min + i * inv_logL_step;
        double L_target = std::exp(logL_target);

        // Находим интервал в исходной таблице, содержащий L_target
        // Используем бинарный поиск (L_em монотонно возрастает)
        auto it = std::upper_bound(L_em.begin(), L_em.end(), L_target);
        if (it == L_em.begin()) {
            // меньше минимального – берём крайнюю левую точку
            inv_logT.push_back(logT.front());
            continue;
        }
        if (it == L_em.end()) {
            // больше максимального – берём крайнюю правую точку
            inv_logT.push_back(logT.back());
            continue;
        }

        int idx_right = std::distance(L_em.begin(), it);
        int idx_left = idx_right - 1;

        double logL_left = logL_em[idx_left];
        double logL_right = logL_em[idx_right];
        double logT_left = logT[idx_left];
        double logT_right = logT[idx_right];

        // Линейная интерполяция в логарифмических координатах
        double weight = (logL_target - logL_left) / (logL_right - logL_left);
        double logT_interp = logT_left + weight * (logT_right - logT_left);
        inv_logT.push_back(logT_interp);
    }

    // Сохраняем узлы inv_logL (только для удобства отладки, но для интерполяции они не нужны,
    // так как сетка равномерна). Заполним для полноты.
    inv_logL.resize(inv_N);
    for (int i = 0; i < inv_N; ++i) {
        inv_logL[i] = inv_logL_min + i * inv_logL_step;
    }
}

// Возвращает температуру T по заданному значению L = L_em(T)
double Dust_spectra::getTemperatureFromL(double L) const
{
    // Защита от выхода за границы исходного диапазона L_em
    if (L <= L_em.front()) {
        return L_em_TL;
    }
    if (L >= L_em.back()) {
        return L_em_TR;
    }

    double logL = std::log(L);

    // Защита от выхода за границы обратной таблицы (на случай ошибок округления)
    if (logL <= inv_logL_min) {
        return std::exp(inv_logT.front());
    }
    if (logL >= inv_logL_max) {
        return std::exp(inv_logT.back());
    }

    // Прямой расчёт индекса в равномерной сетке
    double frac = (logL - inv_logL_min) / inv_logL_step;
    int idx = static_cast<int>(frac);
    // Ограничиваем, чтобы idx был в [0, inv_N-2]
    if (idx < 0) idx = 0;
    if (idx >= inv_N - 1) idx = inv_N - 2;

    double t = frac - idx;   // вес для интерполяции (0..1)
    double logT_interp = inv_logT[idx] + t * (inv_logT[idx + 1] - inv_logT[idx]);
    return std::exp(logT_interp);
}

void Dust_spectra::prepare_sca(double T_min, double T_max, int N_temp,
    double lambda_min, double lambda_max, int M_lambda)
{
    cout << "Start: prepare_sca" << endl;
    // 1. Линейная сетка температур
    T_grid_sca.resize(N_temp);
    for (int i = 0; i < N_temp; ++i) 
    {
        T_grid_sca[i] = T_min + (T_max - T_min) * i / (N_temp - 1);
    }

    // 2. Логарифмическая сетка длин волн
    lambda_grid_sca.resize(M_lambda);
    double log_lmin = std::log(lambda_min);
    double log_lmax = std::log(lambda_max);
    for (int j = 0; j < M_lambda; ++j) {
        double t = static_cast<double>(j) / (M_lambda - 1);
        lambda_grid_sca[j] = std::exp(log_lmin + t * (log_lmax - log_lmin));
    }

    // 3. Построение CDF для каждой температуры
    cdf_table_sca.assign(N_temp, std::vector<double>(M_lambda, 0.0));

    for (int i = 0; i < N_temp; ++i) 
    {
        if (i % 10 == 0) cout << "i = " << i << "   from: " << N_temp << endl;
        double T = T_grid_sca[i];
        std::vector<double> pdf(M_lambda, 0.0);

        // Вычисляем ненормированную PDF: K_abs * B_lambda
        for (int j = 0; j < M_lambda; ++j) 
        {
            double lam = lambda_grid_sca[j];
            double k_abs = interpolate_K_abs(lam);
            double B = planck_b_lambda(lam, T);
            pdf[j] = k_abs * B;
        }

        // Интегрируем методом трапеций (сетка неравномерная)
        double integral = 0.0;
        cdf_table_sca[i][0] = 0.0;
        for (int j = 1; j < M_lambda; ++j) 
        {
            double dlam = lambda_grid_sca[j] - lambda_grid_sca[j - 1];
            integral += 0.5 * (pdf[j - 1] + pdf[j]) * dlam;
            cdf_table_sca[i][j] = integral;
        }

        // Нормировка (последний элемент должен быть 1)
        double total = integral;
        if (total <= 0.0) 
        {
            throw std::runtime_error("Zero total probability for T = " + std::to_string(T));
        }
        for (int j = 0; j < M_lambda; ++j) 
        {
            cdf_table_sca[i][j] /= total;
        }
    }

    cout << "END: prepare_sca" << endl;
}

void Dust_spectra::prepare_inverse_sca(int N_prob) 
{
    cout << "Start: prepare_inverse_sca " << endl;
    if (cdf_table_sca.empty())
        throw std::runtime_error("CDF table must be prepared first (call prepare_sca)");

    N_prob_sca = N_prob;
    size_t N_temp = T_grid_sca.size();
    size_t M = lambda_grid_sca.size();

    inv_lambda_table_sca.assign(N_temp, std::vector<double>(N_prob));

    for (size_t i = 0; i < N_temp; ++i) 
    {
        if (i % 10 == 0) cout << "i = " << i << "   from: " << N_temp << endl;
        const auto& cdf = cdf_table_sca[i];
        for (int k = 0; k < N_prob; ++k) 
        {
            double p = (k + 0.5) / N_prob;   // центр интервала

            // Бинарный поиск интервала в CDF (выполняется один раз при подготовке)
            auto it = std::upper_bound(cdf.begin(), cdf.end(), p);
            size_t j = std::distance(cdf.begin(), it);
            if (j == 0) 
            {
                inv_lambda_table_sca[i][k] = lambda_grid_sca.front();
                continue;
            }
            if (j >= M) 
            {
                inv_lambda_table_sca[i][k] = lambda_grid_sca.back();
                continue;
            }
            double p_low = cdf[j - 1];
            double p_high = cdf[j];
            double t = (p - p_low) / (p_high - p_low);

            double lam_low = lambda_grid_sca[j - 1];
            double lam_high = lambda_grid_sca[j];
            // Логарифмическая интерполяция
            inv_lambda_table_sca[i][k] = lam_low * std::pow(lam_high / lam_low, t);
        }
    }

    cout << "End: prepare_inverse_sca " << endl;
}

void Dust_spectra::find_temp_interval_sca(double T_K, size_t& idx_low, size_t& idx_high, double& frac) const 
{
    if (T_K <= T_grid_sca.front()) {
        idx_low = idx_high = 0;
        frac = 0.0;
        return;
    }
    if (T_K >= T_grid_sca.back()) {
        idx_low = idx_high = T_grid_sca.size() - 1;
        frac = 0.0;
        return;
    }
    auto it = std::upper_bound(T_grid_sca.begin(), T_grid_sca.end(), T_K);
    idx_high = std::distance(T_grid_sca.begin(), it);
    idx_low = idx_high - 1;
    frac = (T_K - T_grid_sca[idx_low]) / (T_grid_sca[idx_high] - T_grid_sca[idx_low]);
}

double Dust_spectra::sample_frequency_sca(double uniform_rand, double T_K) const
{
    if (uniform_rand < 0.0) uniform_rand = 0.0;
    if (uniform_rand > 1.0) uniform_rand = 1.0;

    // Если быстрая таблица не подготовлена – используем обычный бинарный поиск
    if (inv_lambda_table_sca.empty()) 
    {
        throw std::runtime_error("Fast inverse table not prepared. Call prepare_inverse_sca first.");
    }

    size_t i_low, i_high;
    double frac;
    find_temp_interval_sca(T_K, i_low, i_high, frac);

    // Функция получения ? по u для одной температуры (индекс i)
    auto sample_lambda_fast = [&](size_t temp_idx, double u) -> double {
        const auto& table = inv_lambda_table_sca[temp_idx]; // N_prob элементов
        // Позиция на равномерной сетке [0, 1]
        double pos = u * (N_prob_sca - 1);
        int k = static_cast<int>(pos);
        if (k < 0) k = 0;
        if (k >= N_prob_sca - 1) k = N_prob_sca - 2;
        double t = pos - k;  // 0..1
        double lam_low = table[k];
        double lam_high = table[k + 1];
        // Логарифмическая интерполяция внутри интервала
        return lam_low * std::pow(lam_high / lam_low, t);
        };

    double lam_low_T = sample_lambda_fast(i_low, uniform_rand);
    double lam_high_T = sample_lambda_fast(i_high, uniform_rand);

    double lambda_cm = lam_low_T + frac * (lam_high_T - lam_low_T);
    if (lambda_cm <= 0.0) lambda_cm = lambda_grid_sca.front();

    return lambda_cm;
}

void Dust_spectra::save_sca(const std::string& filename) const 
{
    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("Cannot open file for writing: " + filename);

    size_t N = T_grid_sca.size();
    size_t M = lambda_grid_sca.size();
    ofs << N << " " << M << "\n";

    for (size_t i = 0; i < N; ++i) ofs << T_grid_sca[i] << (i + 1 == N ? "\n" : " ");
    for (size_t j = 0; j < M; ++j) ofs << lambda_grid_sca[j] << (j + 1 == M ? "\n" : " ");

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            ofs << cdf_table_sca[i][j] << (j + 1 == M ? "\n" : " ");
        }
    }

    // Сохраняем обратную таблицу, если она есть
    int has_inv = inv_lambda_table_sca.empty() ? 0 : 1;
    ofs << has_inv << "\n";
    if (has_inv) {
        ofs << N_prob_sca << "\n";
        for (size_t i = 0; i < T_grid_sca.size(); ++i) {
            for (int k = 0; k < N_prob_sca; ++k) {
                ofs << inv_lambda_table_sca[i][k] << (k + 1 == N_prob_sca ? "\n" : " ");
            }
        }
    }
}

void Dust_spectra::load_sca(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs) throw std::runtime_error("Cannot open file for reading: " + filename);

    size_t N, M;
    ifs >> N >> M;

    T_grid_sca.resize(N);
    for (size_t i = 0; i < N; ++i) ifs >> T_grid_sca[i];

    lambda_grid_sca.resize(M);
    for (size_t j = 0; j < M; ++j) ifs >> lambda_grid_sca[j];

    cdf_table_sca.assign(N, std::vector<double>(M));
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            ifs >> cdf_table_sca[i][j];

    int has_inv;
    ifs >> has_inv;
    if (has_inv) {
        ifs >> N_prob_sca;
        inv_lambda_table_sca.assign(T_grid_sca.size(), std::vector<double>(N_prob_sca));
        for (size_t i = 0; i < T_grid_sca.size(); ++i) {
            for (int k = 0; k < N_prob_sca; ++k) {
                ifs >> inv_lambda_table_sca[i][k];
            }
        }
    }
    else {
        inv_lambda_table_sca.clear();
        N_prob_sca = 0;
    }
}