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
	vector <vector<Yzel*>> Krug_Yzel;         // ���� � ����� (� �������� �������)  
	// ������ ���������� - �� ����� ����� ����� (0 - ��� ������� �����), ������ ���������� - ����� ����� �� ����
	// Krug_Yzel � Yzel_2D ������������ (�������� ������). ������ A_Luch � A2_Luch �� ������������
	vector <vector<Yzel*>> Krug_Yzel_2;       // ���� � ����� (� ������)


	vector<Gran*> All_Gran;
	vector<Gran*> All_boundary_Gran;   // ������ ��������� ������ (����� ������ ��� �������
	// � ��� ��������� ����������� � ����� �����
	// ������������ ��������  Init_boundary_grans()

	vector<Gran*> Gran_inner_area;
	vector<Gran*> Gran_outer_area;
	// ����� ��� �������� ����� ���������� � ������� �������, � ���� �������� ���� �����������

	vector<Gran*> Gran_TS;  // ���� ����������� ������ �� �����, ������� ������� ���������� 
	// (������������ ����������� ������������ ���� �� �����������)
	// ������� ���� ������ ������������!
	vector<Gran*> Gran_HP;  // � HP ����� ������ ����� D �����, ���������� ������� > this->geo->L6
	vector<Gran*> Gran_BS;


	Cell* Cell_Center;               // ��������� ������ � ������ ������� ���������
	// �� �� �� ����� ��������� ������� ����� ����� �� �������� ���������
	// � ��� ����� �������������� H2, H3, H4

	vector<Cell*> All_Cell;           // ��� ������ �����
	vector <vector<Cell*>> Cell_2D;   // ������ � ����� �����, ������� �������� ��������� 2� �����
	vector <vector<Cell*>> Cell_layer_head;   // ������ � �������� ����� ����� (������ ���)
	vector <vector<Cell*>> Cell_layer_tail;   // ������ � ��������� ����� ����� (������ ���)

	vector<Cell*> Cell_inner_area;          
	vector<Cell*> Cell_outer_area;          
	// � ���� ��������, � ������� �� ������, ��� �����������


	vector<Luch*> All_Luch;
	// ���� ��� ���������� ���������� ����� (�� ������������� ���� ���������)
	vector <vector<Luch*>> A_Luch;   //  [Plane][Number]  ������� ����� ���������, � ������� ����� ���, ����� ����� ������ ����
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


	Surfaces* Surf1;  // ��� ���������� ������������ ������� � �������� ����� � ���
	

	vector<Sensor*> Sensors;

	// ��������� ��� �����-�����
	vector < vector<Gran*>> MK_Grans;
	vector <double> MK_Potoks;  // ������ ����� ���� (����� ��������������� ������ ������)

	Setka();
	~Setka();

	void Algoritm(short int alg);

	void Winslow_method(void);  // ���������� winslow ����� ��� ������������ ����� (����� ��� ����� ������ ������������)

	// ****************************************************************************
	///// ���� ��������� ****************************************************************************
	//****************************************************************************

	void New_initial();  // ������� ���������� ���������� ����� �� ������ 2� ����� � ������������ �����
	void New_connect();  // ������� ���������� ���������� �����: �������� ������ � ���������� �����
	
	void New_append_surfaces();  // ���������� ����� �� ������������ �������
	// �������� ������ �� TS, HP, BS 
	// ���������� ��������� ������ � ��������
	// ����� ����� ������� ��� ����� (�� TS ��� ����� ���������)
	// ����� ����� ������� ��� ����� (� ����� ���� ���������)


	void Calculating_measure(unsigned short int st_time);
	// ���������� ��������, ������� � �������� ������, ������� � ������� �����

	Cell* Find_cell_point(const double& x, const double& y, const double& z, short int now, Cell*& previos);
	// ���� ��� previos, �� �� ������ ���� = nullptr

	bool Time_to_vilet(const MK_particle& P, double& time, Gran*& gran);

	void Renumerate(void); // ������������ ��� �������� ����� (���������� � ������ ���������� ����� ��� ��������)
	bool Test_geometr(void);  // ��������� ��������� ����� (��� �� ������� � �������� ���������)

	void Set_luch_parametr();
	// ��������� ����� ����������� ��������� (��������, ���� ��� ���������� �����) - ��� ������� �������
	// ����� ��� ���������� ���������� �����

	void auto_set_luch_geo_parameter(int for_new);
	// �������������� ��������� �������������� ���������� �����:
	// ��������� ��� ������� ���� (�������� ������-�����) ����� ���������� � ��������� �������
	// ������ ������� ����� ���� ������, ������ ���������� � ��������� �� ����� (�� �� ������ ����)

	void Read_old_surface(string name);
	// ��������� ���� ������������ (� ������� ������ ��������� �� ��������)

	void Move_to_surf(Surfaces* Surf);

	void Smooth_head_TS(void);
	// ����������� TS � �������� �������.
	// �������� ������� TS ������������ ��������  a x^2 + b y^2 + c z^2 + d xy + e xz + f yz = 0
	// ������������ ����� ���������� ��������� � ������.
	// ����� ����� TS ��������� � ������������������� "����������"
	// ��� ����, ����� ������� ��������� TS, ���� �������������� ��������� �� �������, 
	// ��� ������������� ��������� ��� �� �������
	void Smooth_head_HP(void);  // �� �� �����, �� ��� ��������

	void Smooth_head_TS2(void); // ������� �������� ���������� � ���������� ����
	void Smooth_head_HP2(void);
	// �� ����������, ��� ��� ���� ����������� �� ����������, �� � ����������� ����� �� ��������
	// �������� ����� �������� � ������ �������� ������������ 


	// ****************************************************************************
	// ������������ ****************************************************************************
	// ****************************************************************************

	// ���������� �������� �� ������
	void Set_Gran_par_for_interpolate(void);


	void Tecplot_print_plane_interpolation(Eigen::Vector3d A,
		Eigen::Vector3d v1, Eigen::Vector3d v2, int l1, int r1, int l2, int r2);

	// ****************************************************************************
	// ������ ****************************************************************************
	// ****************************************************************************

	void Init_boundary_grans(void); // ��������� ����� ����� �������� ����������
	// ������ ������ ��������� ������
	// �������� ������ �� ���������� �������, ������� ��������� ��������,
	// ����� �������� ������ ������ ��� ���������� �������, ������� ���� ��������� �������� 
	// � ����� ����� �� ������� ��������� �������, ������� ��������� ����� ��������

	int determ_zone(Cell* C, short int now); // ���������� ����, � ������� ��������� ������
	// ���� ����  1, 2, 3, 4

	void Set_MK_Zone(void);

	void Calc_sourse_MF(Cell* C, boost::multi_array<double, 2>& SOURSE,
		short int now, short int zone);

	void Snos_on_Gran(Gran* gr, unordered_map<string, double>& par_left,
		unordered_map<string, double>& par_right, short int now);

	void Culc_Velocity_surface(short int now, const double& time, short int metod = 1);

	void Init_physics(void); // ��������� ��������� �������� ���������� � ������� � ��������� �� ������
	// ������������ �������� ��������� ������� �� ������ (������ ���� ����� ��� ������������)

	void Init_TVD(void);
	// ������������� ������� ����� ��� ��� ���������

	void Go(bool is_inner_area, size_t steps__, short int metod = 1); // ������ �������

	double Culc_Gran_Potok(Gran* gr, unsigned short int now, 
		short int metod, string& name, const double& time);  // ����������� ����� ����� �����
	// ��� ������ ����������� ������
	// ���������� ��� �� �������

	void Save_cell_parameters(string filename);
	void Download_cell_parameters(string filename);

	void Save_for_interpolate(string filename);

	// ����� ����� ***********************************************************
	void MK_prepare(short int zone_MK); // ��������� ����� ��� �����-�����
	void MK_delete(short int zone_MK); 
	void MK_go(short int zone_MK);      // ������ ���� ������
	void MK_fly_immit(MK_particle& P, short int zone_MK, Sensor* Sens);  // ������ �������, ������������ �����
	void M_K_Change_Velosity(Sensor* sens, const double& Ur, const double& Uthe,
		const double& Uphi, const double& Vr, const double& Vthe,
		const double& Vphi, double& Wr, double& Wthe, double& Wphi, const double& cp);
	// ��������� ������� ������������� � ����� � ����������� ������

	//  ****************************************************************************
	// ������������ ****************************************************************************
	//	****************************************************************************
	// ���� ��, ��� �������� ������������ �����

	// ��� Tecplot
	void Tecplot_print_all_yzel_in_3D(string name);
	// �������� ��� ���� (�� �� �����, ����� ����� ���� ������ �� �����������)

	void Tecplot_print_cut_plane_parameters(const Eigen::Vector3d & A, 
		const Eigen::Vector3d & v1,
		const Eigen::Vector3d & v2);
	// ������ ������� 2� ������.
	// ��������� ������ ������ ���������� � �������� ������ ���������
	// �� ��������� ���������� � �����. ������� ��������� ������ �� ��� �����

	void Tecplot_print_all_yzel_with_condition();
	// �������� ���� � �����-�� ��������

	void Tecplot_print_all_lush_in_3D(string name);
	void Tecplot_print_all_cell_in_3D();

	void Tecplot_print_cell_plane_parameters(); 
	// ������� �������� �����������, �������� ������ ������ ����� � �������� ������ � ��� 
	// (��� ����� �� ������ ��������� ��������)



	void Tecplot_print_krug_yzel_in_3D(int num);
	void Tecplot_print_opor_yzel_in_luchs_3D(string name);
	void Tecplot_print_all_lush_in_2D();
	// �������� ���� � 2� ���������, �� ���� ��������� ��������, ����� ������� �������� ���� � ������ ����������
	// �.�. ������������� ���������� ��� ����� (����� ����� � �������� � ��������� �������)
	void Tecplot_print_plane_lush(int plane); // ��� ����������, �� �������� � ���������� ���������
	void Tecplot_print_plane_surfase(int plane); // ��� ����������, �� �������� � ���������� ���������
	void Tecplot_print_All_surfase_in_2D(); // ��� ����������, �� �������� � ���������� ���������

	// �����
	void Tecplot_print_all_gran_in_cell(); // �������� ��� ����� � 3� (����� ������������� �� � ��������)
	// ������������ ���������� ������ � ������ (����� ���� ����� 200 ��)
	void Tecplot_print_all_gran_in_surface(string name_surf);
	// �������� ����� �� �����������

	void Tecplot_print_gran_with_condition(short int iii);

	// *************************************************
	// High_level_of_visualization
	// *************************************************

	void Tecplot_print_1D(Interpol* Int1, const Eigen::Vector3d& Origin,
		const Eigen::Vector3d& vec, string name, const double& leng);

	void Tecplot_print_2D(Interpol* Int1, const double& a, const double& b,
		const double& c, const double& d,
		string name);
	// ���������   a x + b y + c z + d = 0;
};

