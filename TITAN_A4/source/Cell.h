#pragma once
#include"Header.h"

// Ячейкам задаются зоны в функции New_append_surfaces
enum class Type_cell {
	none,   // 0
	Zone_1,  // 1   
	Zone_2,  // 2  
	Zone_3,  // 3   
	Zone_4  // 4   
};


class Cell
{
public:
	vector<Yzel*> yzels;
	vector<Gran*> grans;
	//vector<Edge*> edges;          // Вектор рёбер нужен только для правильного вычисления ротеров
	

	int number = 0;                // номера начинаются с единицы
	bool is_inner = false;         // Является ли ячейка внутренней (которые считаются отдельно)
	bool is_TVD = true;         // Можно ли делать TVD в этой ячейке
	// Фактически равно false только для первых ячеек  вблизи нуля

	mutex mut;
	short int is_need = 0;          // параметр, который можно исползовать для разных нужд

	Type_cell type = Type_cell::none;  // по умолчанию создаём обычный узел

	short unsigned int MK_zone_r = 0;
	short unsigned int MK_zone_phi = 0;
	short unsigned int MK_zone = 0;

	std::array < unordered_map<string, double>, 2> parameters;   // Параметры в ячейке (тут могут быть значения
	// плазменных полей, значения дивергенций и т.д.
	// "rho", 'p', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz'
	// 'Q' - маячок для переноса и определения HP

	unordered_map<string, double> geo_parameters;  // Геометрические параметры для 
	// Удобства и ускорения расчёта Монте-Карло
	// l_size  характерный размер ячейки

	unordered_map<string, Eigen::VectorXd> interpolate_alpha;

	// Источники рождения пикапов (их должно быть в ячейке столько же, сколько у нас есть сортов пикапов в этой зоне)
	vector<double> pui_Sm;   // (n)
	Eigen::MatrixXd pui_Sp;   // (2, n)


	Eigen::MatrixXd mas_pogl;      // (sort, n)    Массив поглощения

	// Функции распределения пикапов
	vector<double> f_pui_1;   // (n)
	vector<double> f_pui_2;   // (n)


	vector<double> F_integr_pui_1;   // (pui_F_n)
	vector<double> nu_integr_pui_1;   // (pui_F_n)
	vector<double> Mz_integr_pui_1;   // (pui_F_n)
	vector<double> E_integr_pui_1;   // (pui_F_n)
	vector<double> F_integr_pui_2;   // (pui_F_n)
	vector<double> nu_integr_pui_2;   // (pui_F_n)
	vector<double> Mz_integr_pui_2;   // (pui_F_n)
	vector<double> E_integr_pui_2;   // (pui_F_n)

	//vector<double> pui_Sm;
	//vector<double> pui_Sp;

	void Init_mas_pogl(short int n, short int nH);  // Инициализация массива поглощения
	void Delete_mas_pogl(); 
	void write_mas_pogl_ToFile(Phys_param* phys_param);
	void read_mas_pogl_FromFile(Phys_param* phys_param);
	int pogl_mas_number(const double& Ve, const double& L, const double& R, short int n);
	


	void Init_pui_integral(short int n, short int zone);  // Инициализация интеграллов для розыгрыша pui
	void Delete_pui_integral(void);
	void write_pui_integral_ToFile(void);
	void read_pui_integral_FromFile(Phys_param*& Phys_param);
	void pui_integral_Culc(Phys_param* phys_param);
	void print_nu_integr_pui(Phys_param* phys_param, string name = "");
	void print_F_integr_pui(string name = "");


	void Init_f_pui(short int n, short int zone);  // Инициализация f_pui, заполняет нулями
	void Delete_f_pui(void);
	void write_pui_ToFile(void);
	void read_pui_FromFile(void);
	void print_pui(double Wmax, string nam);
	void culc_pui_n_T(const double& pui_wR);
	
	double pui_get_f(const double& w, short int ii, const double& Wmax);  // Возвращает значение f_pui в данной точке // TODO!!!
	/// @brief Возвращает значение функции распределения пикапов в заданной точке
	/// 
	/// Функция выполняет интерполяцию функции распределения пикапов (pickup ions) 
	/// для получения значения в произвольной точке энергетического спектра.
	/// Используется для расчетов методом Монте-Карло при моделировании движения 
	/// и взаимодействия пикап-ионов в плазме.
	///
	/// **Принцип работы:**
	/// - Функция работает с дискретными массивами f_pui_1 и f_pui_2, которые содержат
	///   значения функции распределения на сетке по энергии/скорости
	/// - Выполняет интерполяцию между соседними узлами сетки для получения
	///   промежуточных значений
	/// - Нормировка на максимальную энергию Wmax обеспечивает корректное масштабирование
	///
	/// **Область применения:**
	/// - Генерация случайных скоростей пикап-ионов по заданному распределению
	/// - Вычисление моментов функции распределения
	/// - Расчет источников и стоков в уравнениях переноса пикапов
	/// - Моделирование взаимодействия пикапов с фоновой плазмой
	///
	/// @param w Энергия или скорость, для которой нужно получить значение функции распределения
	///          Должна быть в диапазоне [0, Wmax] для корректной интерполяции
	/// @param ii Индекс типа пикап-иона (1 или 2):
	///           - 1: использует массив f_pui_1 (первый тип пикапов)  
	///           - 2: использует массив f_pui_2 (второй тип пикапов)
	///           Соответствует различным сортам ионов или разным энергетическим группам
	/// @param Wmax Максимальная энергия/скорость в сетке функции распределения
	///             Используется для нормировки и определения шага сетки
	///             Должна соответствовать верхней границе массивов f_pui_1/f_pui_2
	///
	/// @return Интерполированное значение функции распределения f(w)
	///         Возвращает 0.0 если w выходит за пределы допустимого диапазона
	///         или если массивы функций распределения не инициализированы
	///
	/// @note Функция требует предварительной инициализации массивов f_pui_1 и f_pui_2
	///       через вызов Init_f_pui() или загрузки данных из файла read_pui_FromFile()
	/// @note Для корректной работы массивы должны быть заполнены актуальными значениями
	///       функции распределения, полученными из расчетов или экспериментальных данных
	/// @note TODO комментарий указывает на необходимость доработки или оптимизации функции
	///
	/// @see Init_f_pui() - инициализация массивов функций распределения
	/// @see read_pui_FromFile() - загрузка функций распределения из файла  
	/// @see culc_pui_n_T() - вычисление моментов функции распределения
	/// @see MK_Add_pui_source() - добавление источника пикапов при моделировании МК
	
	
	double pui_get_nu(const double& w, short int ii, const double& Wmax);
	double PUI_get_F_integer(const double& ksi, short int ii);

	void MK_pui_charge_exchange_velocity(Sensor* sens, Setka* SS, Phys_param* Phys,
		const double& Upx, const double& Upy,
		const double& Upz, const double& UHx, const double& UHy, const double& UHz,
		double& VHx, double& VHy, double& VHz, short int num_pui);
	
	
	void Init_S(short int k, short int n);   // Инициализация S+ S-, заполняет нулями
	void write_S_ToFile(void);
	void read_S_FromFile(const double& n_H_lism);
	void print_SmSp(double Wmax, string nam, Phys_param*& phys);

	void Get_RBF_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par);
	void Get_IDW_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par, Phys_param* phys_param);

	Cell* Get_Sosed(Gran* gr); // Получить соседа через данную грань
	// возвращает nullptr, если соседа через данную грань нет, т.е. грань граничная
	// Для работы функции номера ячеек должны быть актуальны

	friend Cell* Get_Sosed(Cell* C, Gran* gr);
	// Как предыдущая, но сравнение идёт по указателю, не требует номеров ячеек

	bool Belong_point(const double& x, const double& y, const double& z, short int now, bool fast, Cell*& Next);
	/// @brief Проверяет принадлежность точки ячейке
	/// 
	/// Функция определяет, находится ли заданная точка с координатами (x, y, z) 
	/// внутри данной ячейки. Может работать в двух режимах - быстром и точном.
	///
	/// **Точный режим (fast = false):**
	/// - Выполняет полную геометрическую проверку принадлежности точки ячейке
	/// - Разбивает каждую грань ячейки на 2 треугольника 
	/// - Для каждого треугольника строит тетраэдр с центром ячейки как вершиной
	/// - Проверяет принадлежность точки каждому тетраэдру методом барицентрических координат
	/// - Возвращает true, если точка принадлежит хотя бы одному тетраэдру
	/// - Более точный, но медленный алгоритм
	///
	/// **Быстрый режим (fast = true):**
	/// - Выполняет приближенную проверку на основе положения точки относительно граней
	/// - Для каждой грани вычисляет скалярное произведение нормали на вектор от центра грани к точке
	/// - Если точка находится "снаружи" относительно какой-либо грани (скалярное произведение < -1e-10),
	///   то устанавливает Next на соседнюю ячейку через эту грань и возвращает false
	/// - Если точка "внутри" относительно всех граней, возвращает true
	/// - Быстрый, но менее точный алгоритм (имеет погрешность)
	///
	/// @param x Координата точки по оси X
	/// @param y Координата точки по оси Y  
	/// @param z Координата точки по оси Z
	/// @param now Индекс временного слоя (0 или 1) для выбора актуальных координат узлов
	/// @param fast Режим проверки: true - быстрая приближенная проверка, false - точная геометрическая
	/// @param Next Ссылка на указатель соседней ячейки для продолжения поиска
	///            - Устанавливается только в быстром режиме при выходе точки за пределы ячейки
	///            - В точном режиме всегда остается nullptr
	///            - Указывает на соседнюю ячейку через грань, которую "пересекла" точка
	///
	/// @return true если точка принадлежит ячейке, false в противном случае
	///
	/// @note Быстрый режим предназначен для алгоритмов трассировки траекторий частиц
	/// @note Точный режим используется для надежного определения принадлежности в сложных случаях
	/// @note При fast=true параметр Next помогает направить поиск в нужную соседнюю ячейку
	/// @see Set_Cell_Geo_for_MK() - функция должна быть вызвана для корректной работы


	void Set_Cell_Geo_for_MK(void);
	/// @brief Устанавливает геометрические параметры ячейки для метода Монте-Карло
	/// 
	/// Функция вычисляет и сохраняет в geo_parameters основные геометрические характеристики 
	/// ячейки, необходимые для эффективного выполнения расчетов методом Монте-Карло.
	/// Рассчитываются граничные координаты (bounding box) ячейки по всем трем осям.
	///
	/// Вычисляемые параметры:
	/// - "x_min", "x_max" - минимальная и максимальная координаты по оси X
	/// - "y_min", "y_max" - минимальная и максимальная координаты по оси Y  
	/// - "z_min", "z_max" - минимальная и максимальная координаты по оси Z
	/// - "l_size" - характерный размер ячейки (для оценки шага интегрирования)
	///
	/// Эти параметры используются для:
	/// - Быстрой проверки принадлежности точки ячейке
	/// - Оптимизации алгоритмов трассировки лучей
	/// - Генерации случайных точек внутри ячейки
	/// - Ускорения поиска пересечений с границами ячейки
	///
	/// @note Функция должна вызываться после инициализации узлов ячейки (yzels)
	/// @note Результаты сохраняются в geo_parameters для последующего использования
	/// @see Gran::Set_Gran_Geo_for_MK() - аналогичная функция для граней
	/// @see Belong_point() - использует вычисленные параметры для проверки принадлежности



	double center[2][3];           // Центр грани (также в предыдущий и следуюoий момент времени)
	double volume[2];


	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 1);

	double func_R(unsigned short int i_time); // Расстояние от центра до начала координат



	void Tecplot_print_cell(void);

	void MK_Add_moment(MK_particle& P, const double& cp, const double& u, const double& mu_ex,
		const double& u1, const double& u2, const double& u3, const double& skalar, Phys_param* phys_param);

	void MK_Add_particle(MK_particle& P, const double& time, Phys_param* phys_param);
	void MK_Add_pui_source(MK_particle& P, const double& wr, const double& nu_ex, const double& mu,
		const double& time, Phys_param* phys_param, short int zone, short int parent);

	void MK_calc_Sm(Phys_param* phys_param);
	void MK_normir_Moments(Phys_param* phys_param);

};

