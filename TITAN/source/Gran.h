#pragma once
#include "Header.h"

// ��������� ������������� ��� ���������� ��������� �������
enum class Type_Gran {
	Us,    // 0    ������� �����
	Inner_Hard,      // 1       ��������� ����� �� ���������� ����� � ������� ���������� ���������
	Outer_Hard,      // 2       ��������� ����� �� ������� ������� � ������� ���������� ���������
	Outer_Soft      // 3       ��������� ����� �� ������� ������� � ������� ���������� ���������
};

enum class Type_Gran_surf {
	Us,    // 0    ������� �����
	TS,      // 1       ��������� ����� �� ���������� ����� � ������� ���������� ���������
	HP,      // 2       ��������� ����� �� ������� ������� � ������� ���������� ���������
	BS      // 3       ��������� ����� �� ������� ������� � ������� ���������� ���������
};

class Gran
{
public:
	vector<Yzel*> yzels;  // ���� ����� (���� ������ ���� ����������� �� �����)
	vector<Cell*> cells;  // ��� ������ ����� (������� ������ �������� �� ������ �� ������
	// ������������ ����� �� ����� ���������� ��� ����, ����� ����� ��������� ������ ��� ����� 
	// ������������ ����� �����
	vector<Cell*> cells_TVD;  // ������ ��� ����� ��� ��������� �� �����

	vector<Edge*> edges;


	vector<Gran*> grans_surf;  // �����-������, ������������� ������������ �������
	// ���� ������ �������� ������ ��� ������, ������������� ������������
	// �� ����������� � ������� 


	bool work1 = false;            // ������� ������, ������������ ��� ������ ������ (����� �������������� ���� ������ ��������)
	int number = 0;                // ������ ���������� � �������
	Type_Gran type = Type_Gran::Us;  // �� ��������� ������ �������
	Type_Gran_surf type2 = Type_Gran_surf::Us;  // �� ��������� ������ �������

	vector<short int> MK_type;

	unordered_map<string, double> parameters;   // ��������� �� ����� (��� ����� ���� ��������
	// ���������� ����� (������� �������� � �.�.), ����� ���� �������� �������

	unordered_map<string, double> geo_parameters;

	short int Get_method(); // ����� ������� Laks/HLL/HLLC/HLLD ������� ������ ����� �����

	double normal[2][3];           // ������� ����� (����� � ���������� � ������o�� ������ �������)
	double center[2][3];           // ����� ����� (����� � ���������� � ������o�� ������ �������)
	double area[2];           // ������� ����� (����� � ���������� � ������o�� ������ �������)

	// ��� �����-�����
	vector<array<AMR_f*, 2>> AMR;
	// ������ ���������� �� ������� �����, ������ ������ �������
	double MK_Potok;  // ��������� ����� �� ���� ������ � �����
	unsigned short int N_particle = 0; // ����� ������, �������� � �����
	mutex mut;
	
	void Read_AMR(short int ni, short int nH, bool need_refine = false);

	bool Have_zone_number(short int z);
	// ��������� ���� �� � ������� MK_type ���� � ������� z?

	void Culc_measure(unsigned short int st_time);
	// ��������� normal, center, area
	// ��� ���������� ������ ����������� ������� ������ ����� ������ ���� ���������

	double culc_velosity(short int now1, const double& time);
	double func_R(unsigned short int i_time); // ���������� �� ������ �� ������ ���������

	// ������� ��������� �� ������ �����
	friend bool areCellsEqual(const Gran& cell1, const Gran& cell2);
	friend bool areCellsEqual(const Gran* cell1, const Gran* cell2);
	friend bool areCellsEqual_my(const Gran* cell1, const Gran* cell2);

	// ��� �����-�����
	void Get_Random_pozition(Eigen::Vector3d& poz, Sensor* Sens);
	// �������� ��������� ������� �� �����

	void Set_Gran_Geo_for_MK(void);
	// ��������� ����������� ��������� ��� ��
	// x_min, x_max, y_min, y_max, z_min, z_max - ��� �������� ����������� ���������� � ������
	//

	bool Luch_iz_cross_approx(const Eigen::Vector3d& R, const Eigen::Vector3d& V);
	// ���������� �� ��� ������ ����� (��������������, �� ��������������� ������ �����).

	bool Luch_crossing(const Eigen::Vector3d& R, const Eigen::Vector3d& V, double& time);

	bool rayTriangleIntersect(
		const Eigen::Vector3d& orig, const Eigen::Vector3d& dir,
		const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
		double& t);

	Gran();
};

