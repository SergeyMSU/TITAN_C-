#include "Phys_param.h"

Phys_param::Phys_param()
{
	this->Matr << -0.9958639688067077,  0.01776569097515556,  0.08910295088675518,
		           0.07561695085992419, 0.7057402284561812,   0.7044237408557894,
		          -0.05036896241933166, 0.7082479157489926,  -0.7041646522383864;

	this->Matr2 << -0.9958639688067080, 0.0756169508599243, -0.0503689624193315,
	            	0.0177656909751554, 0.7057402284561816,  0.7082479157489927,
		            0.0891029508867553, 0.7044237408557898, -0.7041646522383865;

    // Открываем файл для чтения
    std::ifstream file("nVT1au_new_2000-2022av.dat");
    if (!file.is_open()) {
        std::cerr << "Error 6564397608" << std::endl;
        exit(-1);
    }

    std::string line;
    bool is_first_line = true;

    // Читаем файл построчно
    while (std::getline(file, line)) {
        // Пропускаем первую строку (заголовок)
        if (is_first_line) {
            is_first_line = false;
            continue;
        }

        // Разбиваем строку на отдельные значения
        std::istringstream iss(line);
        double hlat, np, v, t;

        // Если строка содержит данные, записываем их в векторы
        if (iss >> hlat >> np >> v >> t) {
            this->heliolat_deg.push_back(hlat);
            this->n_p_cm3.push_back(np);
            this->V_kms.push_back(v);
            this->T_K.push_back(t);
        }
    }

    file.close();
}

double Phys_param::Get_rho_0(const double& the)
{
    const auto& angles = this->heliolat_deg;
    const auto& np = this->n_p_cm3;

    // Проверка на выход за пределы диапазона
    if (the < angles.front() || the > angles.back()) {
        cout << "Ygol " + std::to_string(the) + " out of diapazon!" << endl;
        cout << "ERROR  0989767567" << endl;
    }

    // Поиск ближайших точек
    size_t i = 0;
    while (i < angles.size() - 1 && angles[i + 1] < the) {
        i++;
    }

    // Если угол точно совпадает с одним из значений в файле
    if (std::abs(the - angles[i]) < 1e-6) {
        return np[i];
    }

    // Линейная интерполяция: n_p = np1 + (theta - theta1) * (np2 - np1) / (theta2 - theta1)
    double theta1 = angles[i];
    double theta2 = angles[i + 1];
    double np1 = np[i];
    double np2 = np[i + 1];

    return (np1 + (the - theta1) * (np2 - np1) / (theta2 - theta1)) / this->char_rho;
}

double Phys_param::Get_v_0(const double& the)
{
    const auto& angles = this->heliolat_deg;
    const auto& np = this->V_kms;

    // Проверка на выход за пределы диапазона
    if (the < angles.front() || the > angles.back()) {
        cout << "Ygol " + std::to_string(the) + " out of diapazon!" << endl;
        cout << "ERROR  0989767567" << endl;
    }

    // Поиск ближайших точек
    size_t i = 0;
    while (i < angles.size() - 1 && angles[i + 1] < the) {
        i++;
    }

    // Если угол точно совпадает с одним из значений в файле
    if (std::abs(the - angles[i]) < 1e-6) {
        return np[i];
    }

    // Линейная интерполяция: n_p = np1 + (theta - theta1) * (np2 - np1) / (theta2 - theta1)
    double theta1 = angles[i];
    double theta2 = angles[i + 1];
    double np1 = np[i];
    double np2 = np[i + 1];

    return (np1 + (the - theta1) * (np2 - np1) / (theta2 - theta1)) / this->char_v;

    return 0.0;
}
