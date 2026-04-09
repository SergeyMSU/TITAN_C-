#include "CoolingFunction.h"


// Чтение данных из файла
void CoolingFunction::ReadCoolingFunction(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        exit(-1);
    }

    std::string line;
    // Пропускаем заголовок (первая строка)
    std::getline(file, line);

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        double T, Lambda;
        if (!(iss >> T >> Lambda)) {
            std::cerr << "Warning: skipped invalid line: " << line << std::endl;
            continue;
        }
        T_values.push_back(T);
        Lambda_values.push_back(Lambda);
    }
    file.close();

    if (T_values.empty()) {
        std::cerr << "Error: no data read from file." << std::endl;
        exit(-1);
    }

    // Предварительное вычисление логарифмов
    lnT.reserve(T_values.size());
    lnLambda.reserve(Lambda_values.size());
    for (size_t i = 0; i < T_values.size(); ++i) {
        lnT.push_back(std::log(T_values[i]));
        lnLambda.push_back(std::log(Lambda_values[i]));
    }

    // Проверка монотонности по логарифму T
    if (!std::is_sorted(lnT.begin(), lnT.end())) {
        std::cerr << "Warning: T values are not in increasing order." << std::endl;
    }

    // Запись проверочного файла
    WriteCheckFile(filename);
}

// Интерполяция Lambda/n_H^2 по температуре T
double CoolingFunction::InterpolateCooling(double T) const
{
    // Если T вне диапазона данных — возвращаем 0
    if (T <= T_values.front() || T >= T_values.back()) {
        return 0.0;
    }

    double logT = std::log(T);
    auto it = std::lower_bound(lnT.begin(), lnT.end(), logT);
    size_t idx = it - lnT.begin();

    // Точное совпадение с узлом (с учётом погрешности)
    if (idx < lnT.size() && std::fabs(*it - logT) < 1e-12) {
        return Lambda_values[idx];
    }

    // Линейная интерполяция в log-log пространстве
    const double& lnT0 = lnT[idx - 1];
    const double& lnT1 = lnT[idx];
    const double& lnL0 = lnLambda[idx - 1];
    const double& lnL1 = lnLambda[idx];

    double lnL = lnL0 + (logT - lnT0) * (lnL1 - lnL0) / (lnT1 - lnT0);
    return std::exp(lnL);
}


// Создание проверочного файла с интерполяцией на равномерной логарифмической сетке
void CoolingFunction::WriteCheckFile(const string& name) const
{
    std::ofstream check_file(name + "_check.txt");
    if (!check_file.is_open()) {
        std::cerr << "Warning: could not create check file." << std::endl;
        return;
    }

    double T_min = T_values.front();
    double T_max = T_values.back();
    double logT_min = std::log(T_min);
    double logT_max = std::log(T_max);
    int N_steps = 10000;
    double log_step = (logT_max - logT_min) / (N_steps - 1);

    check_file << "# T [K]    Lambda_interpolated [erg cm^3/s]\n";
    check_file << std::scientific << std::setprecision(8);
    for (int i = 0; i < N_steps; ++i) {
        double T = std::exp(logT_min + i * log_step);
        double lambda_interp = InterpolateCooling(T);
        check_file << T << "  " << lambda_interp << "\n";
    }

    std::cout << "Check file written: check_cooling_interp.txt" << std::endl;
}