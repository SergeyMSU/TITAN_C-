#pragma once
#include "Header.h"

using namespace std;



class Setka
{
public:
	class Geo_param* geo;
	class Phys_param* phys_param;
	


	bool regim_otladki = true;


	vector<Yzel*> All_Yzel;
	vector <vector<Yzel*>> Yzel_2D;
	vector <vector<Yzel*>> Krug_Yzel;         // Узлы в кругу (в головной области)  
	// первая координата - на какой сфере точке (0 - для ближник точек), вторая координата - номер точки на луче
	// Krug_Yzel и Yzel_2D пересекаются (крайтими узлами). Однако A_Luch и A2_Luch не пересекаются
	vector <vector<Yzel*>> Krug_Yzel_2;       // Узлы в кругу (в хвосте)


	vector<Gran*> All_Gran;
	vector<Gran*> All_boundary_Gran;   // Список граничных граней (какой именно тип границы
	// и его обработка указывается в самой грани
	// определяется функцией  Init_boundary_grans()

	vector<Gran*> Gran_inner_area;
	vector<Gran*> Gran_outer_area;
	// Грани для отдельно счёта внутренней и внешней области, у этих множеств есть пересечение

	vector<Gran*> Gran_TS;  // Сюда добавляются только те грани, которые реально выделяются 
	// (невыделяемое продолжение поверхностей сюда не добавляется)
	// Порядок этих граней произвольный!
	vector<Gran*> Gran_HP;  // В HP вошла только часть D лучей, координаты которых > this->geo->L6
	vector<Gran*> Gran_BS;


	Cell* Cell_Center;               // Фиктивная ячейка в центре системы координат
	// На неё не будут ссылаться никакие грани чтобы не нарушать алгоритмы
	// В ней будет рассчитываться H2, H3, H4

	vector<Cell*> All_Cell;           // Все ячейки сетки
	vector <vector<Cell*>> Cell_2D;   // Ячейки в части сетки, которая получена вращением 2Д сетки
	vector <vector<Cell*>> Cell_layer_head;   // Ячейки в головной части сетки (вблизи оси)
	vector <vector<Cell*>> Cell_layer_tail;   // Ячейки в хвостовой части сетки (вблизи оси)

	vector<Cell*> Cell_inner_area;          
	vector<Cell*> Cell_outer_area;          
	// У этих множеств, в отличие от граней, нет пересечения


	vector<Luch*> All_Luch;
	// Лучи для конкретной реализации сетки (не универсальный блок программы)
	vector <vector<Luch*>> A_Luch;   //  [Plane][Number]  сначала номер плоскости, в которой лежит луч, затем номер самого луча
	vector <vector<Luch*>> B_Luch;
	vector <vector<Luch*>> C_Luch;
	vector <vector<Luch*>> D_Luch;
	vector <vector<Luch*>> E_Luch;
	vector <vector<Luch*>> H_Luch;
	vector <vector<Luch*>> G_Luch;
	map <string, vector<vector<Luch*>>*> All_name_luch; 
	vector<string> name_luch;

	vector<Luch*> A2_Luch;
	vector<Luch*> C2_Luch;


	Surfaces* Surf1;  // Для считывания поверхностей разрыва и движения сетки к ним
	

	vector<Sensor*> Sensors;

	// Параметры для Монте-Карло
	vector < vector<Gran*>> MK_Grans;
	vector <double> MK_Potoks;  // Потоки через зоны (через вышеопределённые наборы граней)

	Setka();
	~Setka();

	void Algoritm(short int alg);

	void Winslow_method(void);  // Реализован winslow метод для триангуляции круга (можно для любой фигуры адаптировать)

	// ****************************************************************************
	///// БЛОК ГЕОМЕТРИИ ****************************************************************************
	//****************************************************************************

	void New_initial();  // Функция начального построения сетки из файлов 2Д сетки и триангуляции круга
	void New_connect();  // Функция начального построения сетки: создания граней и связывания ячеек
	
	void New_append_surfaces();  // Определяем грани на поверхностях разрыва
	// Создание граней на TS, HP, BS 
	// связывание граничных граней с соседями
	// Здесь также задаётся тип узлов (на TS или между разрывами)
	// Здесь также задаётся тиа ячеек (к какой зоне относится)


	void Calculating_measure(unsigned short int st_time);
	// Вычисление площадей, центров и нормалей граней, объёмов и центров ячеек

	Cell* Find_cell_point(const double& x, const double& y, const double& z, short int now, Cell*& previos);
	// Если нет previos, то он должен быть = nullptr

	bool Time_to_vilet(const MK_particle& P, double& time, Gran*& gran);

	void Renumerate(void); // Перенумерует все элементы сетки (необходимо в случае добавления новых или удаления)
	bool Test_geometr(void);  // Тестируем геометрию сетки (что всё связано и работает корректно)

	void Set_luch_parametr();
	// Добавляет лучам необходимые параметры (например, углы для радиальных лучей) - это ускорит расчёты
	// Нужна для начального построения сетки

	void auto_set_luch_geo_parameter(int for_new);
	// Автоматическая настройка геометрических параметров сетки:
	// настройка для каждого луча (особенно сплайн-лучей) длины начального и конечного вектора
	// работа функции может быть долгой, должна вызываться в программе не часто (не на каждом шаге)

	void Read_old_surface(string name);
	// Считывает файл поверхностей (в формате старой программы на фортране)

	void Move_to_surf(Surfaces* Surf);

	void Smooth_head_TS(void);
	// Сглаживание TS в головной области.
	// головной участок TS приблежается функцией  a x^2 + b y^2 + c z^2 + d xy + e xz + f yz = 0
	// используется метод наименьших квадратов с весами.
	// далее точки TS двигаются к аппроксимированному "эллипсойду"
	// Это надо, чтобы вручную исправить TS, если неустойчивость развилась на столько, 
	// что поверхностное натяжение уже не спасает
	void Smooth_head_HP(void);  // То же самое, но для контакта

	void Smooth_head_TS2(void); // Попытка сгладить выпирающие и вдавленные узлы
	void Smooth_head_HP2(void);
	// Не получилось, так как если поверхность не радиальная, то и сглаживание почти не работает
	// Возможно будет работать в случае сильного разваливания 


	// ****************************************************************************
	// Интерполяция ****************************************************************************
	// ****************************************************************************

	// Заполнение значений на гранях
	void Set_Gran_par_for_interpolate(void);


	void Tecplot_print_plane_interpolation(Eigen::Vector3d A,
		Eigen::Vector3d v1, Eigen::Vector3d v2, int l1, int r1, int l2, int r2);

	// ****************************************************************************
	// ФИЗИКА ****************************************************************************
	// ****************************************************************************

	void Init_boundary_grans(void); // Объявляет какие грани являются граничными
	// Создаёт список граничных граней
	// Выдяляет ячейки во внутренней области, которые считаются отдельно,
	// Также создаётся список граней для внутренней области, которые тоже считаются отдельно 
	// А также грани на границе внётренней области, которые считаются также отдельно

	int determ_zone(Cell* C, short int now); // Определить зону, в которой находится ячейка
	// Зоны есть  1, 2, 3, 4

	void Set_MK_Zone(void);

	void Calc_sourse_MF(Cell* C, boost::multi_array<double, 2>& SOURSE,
		short int now, short int zone);

	void Snos_on_Gran(Gran* gr, unordered_map<string, double>& par_left,
		unordered_map<string, double>& par_right, short int now);

	void Culc_Velocity_surface(short int now, const double& time, short int metod = 1);

	void Init_physics(void); // Заполняет начальные значения параметров в ячейках и граничные на гранях
	// Предлагается задавать граничные условия на гранях (должно быть проще это обрабатывать)

	void Init_TVD(void);
	// Инициализация соседей грани для ТВД процедуры

	void Go(bool is_inner_area, size_t steps__, short int metod = 1); // Запуск расчёта

	double Culc_Gran_Potok(Gran* gr, unsigned short int now, 
		short int metod, string& name, const double& time);  // Расчитывает поток через грань
	// все случаи реализуются внутри
	// возвращает шаг по времени

	void Save_cell_parameters(string filename);
	void Download_cell_parameters(string filename);

	void Save_for_interpolate(string filename);

	// Монте карло ***********************************************************
	void MK_prepare(short int zone_MK); // Настройка всего для Монте-Карло
	void MK_delete(short int zone_MK); 
	void MK_go(short int zone_MK);      // Запуск всех частиц
	void MK_fly_immit(MK_particle& P, short int zone_MK, Sensor* Sens);  // Запуск частицы, имитационный метод
	void M_K_Change_Velosity(Sensor* sens, const double& Ur, const double& Uthe,
		const double& Uphi, const double& Vr, const double& Vthe,
		const double& Vphi, double& Wr, double& Wthe, double& Wphi, const double& cp);
	// Сохраняет функции распределения в файлы и освобождает память

	//  ****************************************************************************
	// ВИЗУАЛИЗАЦИЯ ****************************************************************************
	//	****************************************************************************
	// Ниже всё, что касается визуализации сетки

	// Для Tecplot
	void Tecplot_print_all_yzel_in_3D(string name);
	// Печатает все узлы (но по слоям, чтобы можно было удобно их просмотреть)

	void Tecplot_print_cut_plane_parameters(const Eigen::Vector3d & A, 
		const Eigen::Vector3d & v1,
		const Eigen::Vector3d & v2);
	// Лучшая функция 2Д вывода.
	// Разрезает каждую ячейку плоскостью и печатает четырёх угольники
	// со значением параметров в узлах. Текплот мгновенно строит их них карты

	void Tecplot_print_all_yzel_with_condition();
	// Печатает узлы с каким-то условием

	void Tecplot_print_all_lush_in_3D(string name);
	void Tecplot_print_all_cell_in_3D();

	void Tecplot_print_cell_plane_parameters(); 
	// Быстрый просмотр результатов, печатает просто центры ячеек и значения плазмы в них 
	// (для ячеек на первой плоскости вращения)



	void Tecplot_print_krug_yzel_in_3D(int num);
	void Tecplot_print_opor_yzel_in_luchs_3D(string name);
	void Tecplot_print_all_lush_in_2D();
	// Печатает лучи в 2Д плоскости, но есть временной параметр, меняя который получаем лучи в разных плоскостях
	// т.е. просматриваем фактически всю сетку (кроме части в головной и хвостовой области)
	void Tecplot_print_plane_lush(int plane); // как предыдущая, но печатает в конкретной плоскости
	void Tecplot_print_plane_surfase(int plane); // как предыдущая, но печатает в конкретной плоскости
	void Tecplot_print_All_surfase_in_2D(); // как предыдущая, но печатает в конкретной плоскости

	// ГРАНИ
	void Tecplot_print_all_gran_in_cell(); // печатает все грани в 3Д (можно перелистывать их в текплоте)
	// ограничиваем количество граней в выводе (иначе файл весит 200 мб)
	void Tecplot_print_all_gran_in_surface(string name_surf);
	// печатает грани на поверхности

	void Tecplot_print_gran_with_condition(short int iii);

	// *************************************************
	// High_level_of_visualization
	// *************************************************

	void Tecplot_print_1D(Interpol* Int1, const Eigen::Vector3d& Origin,
		const Eigen::Vector3d& vec, string name, const double& leng);

	void Tecplot_print_2D(Interpol* Int1, const double& a, const double& b,
		const double& c, const double& d,
		string name);
	// Плоскость   a x + b y + c z + d = 0;
};

