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

	// ��������� �������� �������
	vector<double> pui_Sm;
	vector<double> pui_Sp;


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
	void MK_Add_pui_source(const double& wr, const double& nu_ex, const double& mu,
		const double& time, Phys_param* phys_param);

	void MK_calc_Sm(Phys_param* phys_param);
	void MK_normir_Moments(Phys_param* phys_param);

};

