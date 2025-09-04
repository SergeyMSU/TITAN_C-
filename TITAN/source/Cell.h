#pragma once
#include"Header.h"

// ячейкам задаютс€ зоны в функции New_append_surfaces
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
	//vector<Edge*> edges;          // ¬ектор рЄбер нужен только дл€ правильного вычислени€ ротеров
	

	int number = 0;                // номера начинаютс€ с единицы
	bool is_inner = false;         // явл€етс€ ли €чейка внутренней (которые считаютс€ отдельно)
	bool is_TVD = true;         // ћожно ли делать TVD в этой €чейке
	// ‘актически равно false только дл€ первых €чеек  вблизи нул€

	mutex mut;
	short int is_need = 0;          // параметр, который можно исползовать дл€ разных нужд

	Type_cell type = Type_cell::none;  // по умолчанию создаЄм обычный узел

	short unsigned int MK_zone_r = 0;
	short unsigned int MK_zone_phi = 0;
	short unsigned int MK_zone = 0;

	std::array < unordered_map<string, double>, 2> parameters;   // ѕараметры в €чейке (тут могут быть значени€
	// плазменных полей, значени€ дивергенций и т.д.
	// "rho", 'p', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz'
	// 'Q' - ма€чок дл€ переноса и определени€ HP

	unordered_map<string, double> geo_parameters;  // √еометрические параметры дл€ 
	// ”добства и ускорени€ расчЄта ћонте- арло
	// l_size  характерный размер €чейки

	unordered_map<string, Eigen::VectorXd> interpolate_alpha;

	// »сточники рождени€ пикапов (их должно быть в €чейке столько же, сколько у нас есть сортов пикапов в этой зоне)
	vector<double> pui_Sm;
	vector<double> pui_Sp;


	void Get_RBF_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par);
	void Get_IDW_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par, Phys_param* phys_param);

	Cell* Get_Sosed(Gran* gr); // ѕолучить соседа через данную грань
	// возвращает nullptr, если соседа через данную грань нет, т.е. грань гранична€
	// ƒл€ работы функции номера €чеек должны быть актуальны

	friend Cell* Get_Sosed(Cell* C, Gran* gr);
	//  ак предыдуща€, но сравнение идЄт по указателю, не требует номеров €чеек

	bool Belong_point(const double& x, const double& y, const double& z, short int now, bool fast, Cell*& Next);
	// ѕринадлежит ли точка €чейке

	void Set_Cell_Geo_for_MK(void);


	double center[2][3];           // ÷ентр грани (также в предыдущий и следуюoий момент времени)
	double volume[2];


	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 1);

	double func_R(unsigned short int i_time); // –ассто€ние от центра до начала координат


	void Tecplot_print_cell(void);

	void MK_Add_particle(MK_particle& P, const double& time);
	void MK_Add_pui_source(const double& wr, const double& nu_ex, const double& mu,
		const double& time, Phys_param* phys_param);

	void MK_calc_Sm(Phys_param* phys_param);
	void MK_normir_Moments(Phys_param* phys_param);

};

