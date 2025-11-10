#pragma once
#include "Header.h"

// —труктура дл€ данных активных €чеек
struct ActiveCellData 
{
	double f = 0.0;
	double Spotok = 0.0; // Ёто поток в €чейке не умноженный на грань!
	unordered_map<string, double> param;    // ¬спомогательные параметры в €чейке, например дл€ линейного сноса
	// f = A + Bx (Vx - Vx0) + By (Vy - Vy0) + Bz (Vz - Vz0)
	// Bx, By, Bz
	// Max - максимальное значение линейной функции в €чейке (дл€ метода отказов)
	
	struct {
		unsigned is_signif : 1;      // 1 бит // —ущестывенна€ €чейка, та, которую можно делить, если надо 
	// определ€етс€ по процентру плотности к всему объЄму
		unsigned need_devide_x : 1;  // 1 бит  // Ќужно ли еЄ делить вдоль x?
		unsigned need_devide_y : 1;  // 1 бит  // Ќужно ли еЄ делить вдоль y?
		unsigned need_devide_z : 1;  // 1 бит  // Ќужно ли еЄ делить вдоль z?
	} active_flags;
};


class AMR_cell
{
public:
	AMR_cell* I_self;            // ”казатель на себ€

	uint8_t level = 0;
	AMR_cell* parent = nullptr;    // ячейка - родитель
	uint8_t nx = 0;           // Ќомер данной €чейки в €чейке-родителе
	uint8_t ny = 0;
	uint8_t nz = 0;

	struct Flags {
		unsigned is_divided : 1;     // 1 бит  // –азделена ли €чейка
	} flags;  // –азмер: 1 байт

	boost::multi_array<AMR_cell*, 3> cells;  // ячейки - дети

	// ƒанные активных €чеек (может быть nullptr дл€ неактивных €чеек)
	std::unique_ptr<ActiveCellData> active_data;

	AMR_cell();

	// ‘ункции дл€ работы с данными активных €чеек
	void ensure_active_data();  // —оздать active_data если его нет
	void clear_active_data();   // ќчистить active_data дл€ освобождени€ пам€ти
	bool has_active_data() const { return active_data != nullptr; }
	
	// √еттеры и сеттеры дл€ безопасного доступа к данным
	double getF() const 
	{ 
		if (active_data)
		{
			return active_data->f;
		}
		else
		{
			cout << "Error 8934yt807goueyw4hfgiuehgf873egg" << endl;
			return 0.0;
		}
	}
	void setF(double value) { ensure_active_data(); active_data->f = value; }
	
	double getSpotok() const 
	{ 
		if (active_data)
		{
			return active_data->Spotok;
		}
		else
		{
			cout << "Error ergert34t43r5t3345t3" << endl;
			return 0.0;
		}
	}
	void setSpotok(double value) { ensure_active_data(); active_data->Spotok = value; }
	
	bool isSignif() const
	{ 
		if (active_data)
		{
			return active_data->active_flags.is_signif;
		}
		else
		{
			cout << "Error etrhgrtyhr6y456yt46" << endl;
			return false;
		}
	}
	void setIsSignif(bool value) { ensure_active_data(); active_data->active_flags.is_signif = value; }
	
	bool needDevideX() const 
	{ 
		if (active_data)
		{
			return active_data->active_flags.need_devide_x;
		}
		else
		{
			cout << "Error etrhgrety45t45y56yh" << endl;
			return false;
		}
	}
	void setNeedDevideX(bool value) { ensure_active_data(); active_data->active_flags.need_devide_x = value; }
	
	bool needDevideY() const { return active_data ? active_data->active_flags.need_devide_y : false; }
	void setNeedDevideY(bool value) { ensure_active_data(); active_data->active_flags.need_devide_y = value; }
	
	bool needDevideZ() const { return active_data ? active_data->active_flags.need_devide_z : false; }
	void setNeedDevideZ(bool value) { ensure_active_data(); active_data->active_flags.need_devide_z = value; }
	
	bool isDivided() const { return flags.is_divided; }
	void setIsDivided(bool value) { flags.is_divided = value; }
	
	unordered_map<string, double>& getParam() 
	{ 
		if (active_data)
		{
			return active_data->param;
		}
		else
		{
			cout << "Error 564y56y3gregegr" << endl;
			exit(-1);
		}
	}
	/*const unordered_map<string, double>& getParam() const 
	{ 
		static unordered_map<string, double> empty_map;
		return active_data ? active_data->param : empty_map; 
	}*/

	double Get_SpotokV(void);
	void Get_Moment(AMR_f* AMR, double & m, double& mu, double& mux, double& muu);
	void Get_f(AMR_f* AMR, double& S);


	void Cell_partially_free_space(void);
	// ќсвобождает место дл€ €чеек, которые €вл€ютс€ неактивными

	void Re_Cell_partially_free_space(void);
	

	void Culc_gradients(AMR_f* AMR);
	// ¬ычисл€ет градианты дл€ данной €чейки, использу€ еЄ соседей
	// задаютс€ Bx, By, Bz

	void divide(AMR_f* AMR, unsigned short int n1, unsigned short int n2, unsigned short int n3); // –азделить €чейку

	AMR_cell* find_cell(const double& x, const double& y, const double& z, const double& xL, 
		const double& xR, const double& yL, const double& yR, const double& zL, const double& zR);
	// »щет €чейку по еЄ сосед€м

	AMR_cell* get_sosed(AMR_f* AMR, short int nn);
	// nn = 0, 1, 2, 3, 4, 5
	//     по х вперЄд - назад, по y ...

	void Print_info(void);

	void Get_random_velosity_in_cell(AMR_f* AMR, const double& ksi, const double& Squ, Eigen::Vector3d& Vel, Sensor* Sens);

	void Get_index(std::vector<std::array<unsigned int, 3>>& numbers);
	// ѕолучить индекс €чейки

	void Get_Center(AMR_f* AMR, std::array<double, 3>& center); // ѕолучить центр €чейки (даже если она разбита)
	void Get_Center(AMR_f* AMR, std::array<double, 3>& center, std::array<double, 3>& razmer); // ѕолучить центр €чейки (даже если она разбита)

	void Get_Centers(AMR_f* AMR, std::vector<std::array<double, 3>>& centers); // ѕолучить центры €чейки (включа€ центры подъечеек)
	void Get_all_cells(std::vector< AMR_cell*>& cells); // ѕолучить список действительных €чеек (неразделЄнных)

	void Slice_plane(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d, std::vector < std::vector<std::array<double, 3>>>& points);
	// –азрезать €чейку плоскостью

	void Save_cell(std::ofstream& out);
	void Read_cell(std::ifstream& in);

	void Delete(void); 
	// ”дал€ет текущую €чейку и все еЄ соседние
};

