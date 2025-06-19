#pragma once
#include"Header.h"

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
	int number = 0;                // ������ ���������� � �������
	bool is_inner = false;         // �������� �� ������ ���������� (������� ��������� ��������)

	Type_cell type = Type_cell::none;  // �� ��������� ������ ������� ����

	std::array < unordered_map<string, double>, 2> parameters;   // ��������� � ������ (��� ����� ���� ��������
	// ���������� �����, �������� ����������� � �.�.
	// "rho", 'p', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz'
	// 'Q' - ������ ��� �������� � ����������� HP

	unordered_map<string, Eigen::VectorXd> interpolate_alpha;


	void Get_RBF_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par);
	void Get_IDW_interpolation(const double& x, const double& y, const double& z, unordered_map<string, double>& par, Phys_param* phys_param);

	Cell* Get_Sosed(Gran* gr); // �������� ������ ����� ������ �����
	// ���������� nullptr, ���� ������ ����� ������ ����� ���, �.�. ����� ���������
	// ��� ������ ������� ������ ����� ������ ���� ���������

	friend Cell* Get_Sosed(Cell* C, Gran* gr);
	// ��� ����������, �� ��������� ��� �� ���������, �� ������� ������� �����

	bool Belong_point(const double& x, const double& y, const double& z, short int now, bool fast, Cell*& Next);
	// ����������� �� ����� ������

	double center[2][3];           // ����� ����� (����� � ���������� � ������o�� ������ �������)
	double volume[2];


	void Culc_center(unsigned short int st_time);
	void Culc_volume(unsigned short int st_time, unsigned short int method = 1);

	double func_R(unsigned short int i_time); // ���������� �� ������ �� ������ ���������

};

