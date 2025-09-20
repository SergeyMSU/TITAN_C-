#pragma once
#include"Header.h"

// ������� �������� ���� � ������� New_append_surfaces
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
	//vector<Edge*> edges;          // ������ ���� ����� ������ ��� ����������� ���������� �������
	

	int number = 0;                // ������ ���������� � �������
	bool is_inner = false;         // �������� �� ������ ���������� (������� ��������� ��������)
	bool is_TVD = true;         // ����� �� ������ TVD � ���� ������
	// ���������� ����� false ������ ��� ������ �����  ������ ����

	mutex mut;
	short int is_need = 0;          // ��������, ������� ����� ����������� ��� ������ ����

	Type_cell type = Type_cell::none;  // �� ��������� ������ ������� ����

	short unsigned int MK_zone_r = 0;
	short unsigned int MK_zone_phi = 0;
	short unsigned int MK_zone = 0;

	std::array < unordered_map<string, double>, 2> parameters;   // ��������� � ������ (��� ����� ���� ��������
	// ���������� �����, �������� ����������� � �.�.
	// "rho", 'p', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz'
	// 'Q' - ������ ��� �������� � ����������� HP

	unordered_map<string, double> geo_parameters;  // �������������� ��������� ��� 
	// �������� � ��������� ������� �����-�����
	// l_size  ����������� ������ ������

	unordered_map<string, Eigen::VectorXd> interpolate_alpha;

	// ��������� �������� ������� (�� ������ ���� � ������ ������� ��, ������� � ��� ���� ������ ������� � ���� ����)
	vector<double> pui_Sm;   // (n)
	Eigen::MatrixXd pui_Sp;   // (2, n)

	// ������� ������������� �������
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

	void Init_pui_integral(short int n, short int zone);  // ������������� f_pui, ��������� ������
	void Delete_pui_integral(void);
	void write_pui_integral_ToFile(void);
	void read_pui_integral_FromFile(void);
	void pui_integral_Culc(Phys_param* phys_param);


	void Init_f_pui(short int n, short int zone);  // ������������� f_pui, ��������� ������
	void Delete_f_pui(void);
	void write_pui_ToFile(void);
	void read_pui_FromFile(void);
	void print_pui(double Wmax, string nam);
	void culc_pui_n_T(const double& pui_wR);
	double pui_get_f(const double& w, short int ii);  // ���������� �������� f_pui � ������ ����� // TODO!!!

	void Init_S(short int k, short int n);   // ������������� S+ S-, ��������� ������
	void write_S_ToFile(void);
	void read_S_FromFile(void);
	void print_SmSp(double Wmax, string nam);

	void Get_RBF_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par);
	void Get_IDW_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par, Phys_param* phys_param);

	Cell* Get_Sosed(Gran* gr); // �������� ������ ����� ������ �����
	// ���������� nullptr, ���� ������ ����� ������ ����� ���, �.�. ����� ���������
	// ��� ������ ������� ������ ����� ������ ���� ���������

	friend Cell* Get_Sosed(Cell* C, Gran* gr);
	// ��� ����������, �� ��������� ��� �� ���������, �� ������� ������� �����

	bool Belong_point(const double& x, const double& y, const double& z, short int now, bool fast, Cell*& Next);
	// ����������� �� ����� ������

	void Set_Cell_Geo_for_MK(void);


	double center[2][3];           // ����� ����� (����� � ���������� � ������o�� ������ �������)
	double volume[2];


	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 1);

	double func_R(unsigned short int i_time); // ���������� �� ������ �� ������ ���������


	void Tecplot_print_cell(void);

	void MK_Add_particle(MK_particle& P, const double& time);
	void MK_Add_pui_source(MK_particle& P, const double& wr, const double& nu_ex, const double& mu,
		const double& time, Phys_param* phys_param, short int zone, short int parent);

	void MK_calc_Sm(Phys_param* phys_param);
	void MK_normir_Moments(Phys_param* phys_param);

};

