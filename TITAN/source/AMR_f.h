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

	array<double, 3> Vn;
	array<double, 3> Vt;
	array<double, 3> Vm;

	unsigned int xn;
	unsigned int yn;
	unsigned int zn;

	double Sf;
	double Sfu;
	double Sfuu;

	AMR_f* AMR_self;

	boost::multi_array<AMR_cell*, 3> cells;

	AMR_f();
	AMR_f(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
		const double& zR, unsigned int xn, unsigned int yn, unsigned int zn);

	void Set_bazis(void);
	// ѕо нормали определ€ет два других вектора

	AMR_cell* find_cell(const double& x, const double& y, const double& z);
	// »щет €чейку (указатель на неЄ) по координате

	void Get_all_cells(std::vector<AMR_cell*>& cells); 
	// ѕолучить список действительных €чеек (т.е. если €чейка разделена, она не
	// включаютс€, а включаютс€ еЄ дети и т.д.).

	void Fill_test(void);
	// «аполнить €чейки максвеллом

	void Refine(void);

	void Save(string namef);
	void Read(string namef);


	void Print_info(void);

	void Print_all_center_Tecplot(AMR_f* AMR);
	void Print_slice_Tecplot(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d);
	// плоскость  a x + b y + c z + d = 0

	void Print_all_sosed_Tecplot(AMR_f* AMR);

	void Print_1D_Tecplot(AMR_f* AMR);

	void Delete(void);
	// ”дал€ет сетку, очищает пам€ть и т.д.
};

