#pragma once
#include "Header.h"
#include <vector>

class AMR_f
{
public:
	double xL;
	double xR;
	double yL;
	double yR;
	double zL;
	double zR;

	double procent_signif = 0.3;  // 0.3;
	double procent_devide = 1.0;  // 1.0;

	array<double, 3> Vn;
	array<double, 3> Vt;
	array<double, 3> Vm;

	unsigned int xn;
	unsigned int yn;
	unsigned int zn;

	// ��������� ��� ����������� ���� �� ������� �������
	// ���������� �� ������� �������������
	double Sf;
	double Sfu;
	double Sfux;
	double Sfuu;

	unordered_map<string, double> parameters;  // �������������� ��-���������

	double SpotokV = 0.0;  // ����� ����� ������ ����� ������� �������������
	// ���� ����� ��� ������� �� ������� ����� (���� ������ ��� ����� �����)

	AMR_f* AMR_self;
	mutex mut;

	boost::multi_array<AMR_cell*, 3> cells;

	void Culk_SpotokV(const double& Squ);
	// ������� ����� ������� ������������� ����� �������
	// ����� ��������� ����� � ������������ (����������) �������

	void Get_random_velosity(AMR_f* AMR, const double& Squ, Eigen::Vector3d& Vel, Sensor* Sens);

	AMR_f();
	AMR_f(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
		const double& zR, unsigned int xn, unsigned int yn, unsigned int zn);

	void AMR_resize(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
		const double& zR, unsigned int xn, unsigned int yn, unsigned int zn);

	void Get_real_koordinate(const double& x, const double& y, 
		const double& z, double& Vx, double& Vy, double& Vz);
	// �� ��������� ����������� �������� �������� ���������� ����������

	void Get_lokal_koordinate(const double& Vx, const double& Vy, 
		const double& Vz, double& x, double& y, double& z);

	void Add_particle(const double& Vx, const double& Vy,
		const double& Vz, const double& mu);

	void Normir_velocity_volume(const double& squ);
	// ���������� � ����� ������� �����-�����

	void Set_bazis(void);
	// �� ������� ���������� ��� ������ �������

	AMR_cell* find_cell(const double& x, const double& y, const double& z);
	// ���� ������ (��������� �� ��) �� ����������

	void Get_all_cells(std::vector<AMR_cell*>& cells); 
	// �������� ������ �������������� ����� (�.�. ���� ������ ���������, ��� ��
	// ����������, � ���������� � ���� � �.�.).

	void Fill_maxwel_inf(const double& Vinf); 
	// ��������� ����������� �� �������������
	// � �������� ��������� ���������

	double Integrate_Maxwell_V(const double& xL, const double& xR,
		const double& yL, const double& yR,
		const double& zL, const double& zR, const double& Vinf);

	void Fill_test(void);
	// ��������� ������ ����������

	void Fill_null(void);
	// ��������� �������� ������ ����

	unsigned int Refine(void);
	unsigned int de_Refine(void);

	void Save(string namef);
	void Read(string namef);

	unsigned int Size(void);


	void Print_info(void);

	void Print_all_center_Tecplot(AMR_f* AMR, const string& name = "_");
	void Print_slice_Tecplot(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d);
	// ���������  a x + b y + c z + d = 0

	void Print_all_sosed_Tecplot(AMR_f* AMR);

	void Print_1D_Tecplot(AMR_f* AMR, const double& VV);

	void Delete(void);
	// ������� �����, ������� ������ � �.�.
};

