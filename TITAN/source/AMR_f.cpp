#include "AMR_f.h"
using namespace std;

void AMR_f::Culk_SpotokV(const double& Squ)
{
	std::vector<AMR_cell*> cells;
	std::array<double, 3> center;
	std::array<double, 3> razmer;

	this->SpotokV = 0.0;

	this->Get_all_cells(cells);
	double V;
	double Vx, Vy, Vz;

	//cout << "All_cells_do = " << cells.size() << endl;

	for (const auto& i : cells)
	{
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		this->Get_real_koordinate(center[0], center[1], center[2],
			Vx, Vy, Vz);
		i->Spotok = V * i->f * fabs(Vx);
		this->SpotokV += i->Spotok;
	}

	this->SpotokV *= Squ;

	const auto& shape = this->cells.shape();
	const size_t nx = shape[0];
	const size_t ny = shape[1];
	const size_t nz = shape[2];
	for (size_t i = 0; i < nx; ++i)
	{
		for (size_t j = 0; j < ny; ++j)
		{
			for (size_t k = 0; k < nz; ++k)
			{
				AMR_cell* cell = this->cells[i][j][k];
				cell->Spotok = cell->Get_SpotokV();
			}
		}
	}

}

void AMR_f::Get_random_velosity(AMR_f* AMR, const double& Squ, Eigen::Vector3d& Vel, Sensor* Sens)
{
	// Squ - площадь грани
	// Vel - возвращаемая скорость частицы
	// Sens - датчик случайных чисел

	double ksi = Sens->MakeRandom() * this->SpotokV;
	double SS = 0.0;

	const auto& shape = this->cells.shape();
	const size_t nx = shape[0];
	const size_t ny = shape[1];
	const size_t nz = shape[2];
	for (size_t i = 0; i < nx; ++i)
	{
		for (size_t j = 0; j < ny; ++j)
		{
			for (size_t k = 0; k < nz; ++k)
			{
				AMR_cell* cell = this->cells[i][j][k];
				if (SS + cell->Spotok * Squ > ksi)
				{
					// Нам нужна эта ячейка
					cell->Get_random_velosity_in_cell(AMR, ksi - SS, Squ, Vel, Sens);

					double A0 = Vel[0] * this->Vn[0] + Vel[1] * this->Vt[0] +
						Vel[2] * this->Vm[0];
					double A1 = Vel[0] * this->Vn[1] + Vel[1] * this->Vt[1] +
						Vel[2] * this->Vm[1];
					double A2 = Vel[0] * this->Vn[2] + Vel[1] * this->Vt[2] +
						Vel[2] * this->Vm[2];

					Vel[0] = A0;
					Vel[1] = A1;
					Vel[2] = A2;
					return;
				}
				else
				{
					SS += cell->Spotok * Squ;
				}
			}
		}
	}

	cout << "Error 6453211877  " << endl;
	whach(SS);
	whach(ksi);
	whach(this->SpotokV);
	exit(-1);

}

AMR_f::AMR_f()
{
	this->xL = 0.0;
	this->xR = 0.0;

	this->yL = 0.0;
	this->yR = 0.0;

	this->zL = 0.0;
	this->zR = 0.0;

	this->xn = 0;
	this->yn = 0;
	this->zn = 0;

	this->Vn = { 0.0, 0.0, 0.0 };
	this->Vt = { 0.0, 0.0, 0.0 };
	this->Vm = { 0.0, 0.0, 0.0 };
}

AMR_f::AMR_f(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
	const double& zR, unsigned int xn, unsigned int yn, unsigned int zn)
{
	this->xL = xL;
	this->xR = xR;

	this->yL = yL;
	this->yR = yR;

	this->zL = zL;
	this->zR = zR;

	this->xn = xn;
	this->yn = yn;
	this->zn = zn;

	this->Vn = { 0.0, 0.0, 0.0 };
	this->Vt = { 0.0, 0.0, 0.0 };
	this->Vm = { 0.0, 0.0, 0.0 };


	this->cells.resize(boost::extents[xn][yn][zn]);

	for (int i = 0; i < xn; ++i) {
		for (int j = 0; j < yn; ++j) {
			for (int k = 0; k < zn; ++k) {
				auto A = new AMR_cell();
				A->nx = i;
				A->ny = j;
				A->nz = k;
				A->parent = nullptr;
				A->I_self = A;
				A->level = 0;
				this->cells[i][j][k] = A;
			}
		}
	}
}

void AMR_f::Get_real_koordinate(const double& x, const double& y, const double& z, double& Vx, double& Vy, double& Vz)
{
	Eigen::Vector3d ex, ey, ez, ee;
	ex << this->Vn[0], this->Vn[1], this->Vn[2];
	ey << this->Vt[0], this->Vt[1], this->Vt[2];
	ez << this->Vm[0], this->Vm[1], this->Vm[2];

	ee = x * ex + y * ey + z * ez;
	Vx = ee[0];
	Vy = ee[1];
	Vz = ee[2];
	return;
}

void AMR_f::Get_lokal_koordinate(const double& Vx, const double& Vy, const double& Vz, double& x, double& y, double& z)
{
	Eigen::Vector3d ex, ey, ez, ee;
	ex << this->Vn[0], this->Vn[1], this->Vn[2];
	ey << this->Vt[0], this->Vt[1], this->Vt[2];
	ez << this->Vm[0], this->Vm[1], this->Vm[2];
	ee << Vx, Vy, Vz;

	x = ee.dot(ex);
	y = ee.dot(ey);
	z = ee.dot(ez);
	return;
}

void AMR_f::Add_particle(const double& Vx, const double& Vy, const double& Vz, const double& mu)
{
	double x, y, z;
	double Vnn;
	Vnn = fabs(scalarProductFast(Vx, Vy, Vz, this->Vn[0], this->Vn[1], this->Vn[2]));
	if (Vnn < 1e-6) Vnn = 1e-6;
	this->Get_lokal_koordinate(Vx, Vy, Vz, x, y, z);
	unsigned short int kkk = 0;
	a2:
	kkk++;
	AMR_cell* cell = find_cell(x, y, z);

	if (kkk > 100)
	{
		cout << "Error 75463208674" << endl;
		whach(this->Vn[0]);
		whach(this->Vn[1]);
		whach(this->Vn[2]);
		whach(this->Vt[0]);
		whach(this->Vt[1]);
		whach(this->Vt[2]);
		whach(this->Vm[0]);
		whach(this->Vm[1]);
		whach(this->Vm[2]);
		whach(Vx);
		whach(Vy);
		whach(Vz);
		whach(x);
		whach(y);
		whach(z);
		exit(-1);
	}

	if (cell == nullptr)
	{
		x *= 0.995;
		y *= 0.995;
		z *= 0.995;
		goto a2;
	}
	cell->f += mu/Vnn;
}

void AMR_f::Normir_velocity_volume(void)
{
	std::vector<AMR_cell*> cells;
	std::array<double, 3> center;
	std::array<double, 3> razmer;
	this->Get_all_cells(cells);

	for (auto& cel : cells)
	{
		cel->Get_Center(this->AMR_self, center, razmer);
		cel->f /= razmer[0] * razmer[1] * razmer[2];
	}
}

void AMR_f::Set_bazis(void)
{
	Eigen::Vector3d n;
	Eigen::Vector3d t;
	Eigen::Vector3d m;

	n << this->Vn[0], this->Vn[1], this->Vn[2];

	get_bazis(n, t, m);
	this->Vt[0] = t[0];
	this->Vt[1] = t[1];
	this->Vt[2] = t[2];

	this->Vm[0] = m[0];
	this->Vm[1] = m[1];
	this->Vm[2] = m[2];
}

AMR_cell* AMR_f::find_cell(const double& x, const double& y, const double& z)
{
	if(x < this->xL || x > this->xR) return nullptr;
	if(y < this->yL || y > this->yR) return nullptr;
	if(z < this->zL || z > this->zR) return nullptr;

	double dx = (this->xR - this->xL) / this->xn;
	int index1 = static_cast<int>((x - this->xL) / dx);
	if (index1 == this->xn) index1 = this->xn - 1;

	double dy = (this->yR - this->yL) / this->yn;
	int index2 = static_cast<int>((y - this->yL) / dy);
	if (index2 == this->yn) index2 = this->yn - 1;

	double dz = (this->zR - this->zL) / this->zn;
	int index3 = static_cast<int>((z - this->zL) / dz);
	if (index3 == this->zn) index3 = this->zn - 1;

	auto A = this->cells[index1][index2][index3];

	if (A->is_divided == false)
	{
		return A;
	}
	else
	{
		return A->find_cell(x, y, z,
			xL + index1 * dx, xL + (index1 + 1) * dx,
			yL + index2 * dy, yL + (index2 + 1) * dy,
			zL + index3 * dz, zL + (index3 + 1) * dz);
	}
}

void AMR_f::Get_all_cells(vector<AMR_cell*>& cells)
{
	const auto& shape = this->cells.shape();
	const size_t nx = shape[0];
	const size_t ny = shape[1];
	const size_t nz = shape[2];
	for (size_t i = 0; i < nx; ++i)
	{
		for (size_t j = 0; j < ny; ++j)
		{
			for (size_t k = 0; k < nz; ++k)
			{
				AMR_cell* cell = this->cells[i][j][k];
				if (cell->is_divided == false) {
					cells.push_back(cell);
				}
				else
				{
					cell->Get_all_cells(cells);
				}
			}
		}
	}
}

void AMR_f::Fill_maxwel_inf(const double& Vinf)
{
	Eigen::Vector3d ex, ey, ez, ee;
	ex << this->Vn[0], this->Vn[1], this->Vn[2];
	ey << this->Vt[0], this->Vt[1], this->Vt[2];
	ez << this->Vm[0], this->Vm[1], this->Vm[2];

	unsigned int NN = 1;

	while(NN > 0)
	{
		NN = 0;
		std::vector<AMR_cell*> cells;
		std::array<double, 3> center;
		this->Get_all_cells(cells);

		for (auto& cel : cells)
		{
			cel->Get_Center(this->AMR_self, center);
			ee = center[0] * ex + center[1] * ey + center[2] * ez;
			cel->f = maxwell(1.0, 1.0, Vinf, 0.0, 0.0, ee[0], ee[1], ee[2]);
		}

		NN = this->Refine();
	}
}

void AMR_f::Fill_null(void)
{
	std::vector<AMR_cell*> cells;
	this->Get_all_cells(cells);

	for (auto& i : cells)
	{
		i->f = 0.0;
	}
}

void AMR_f::Fill_test(void)
{
	std::vector<AMR_cell*> cells;
	this->Get_all_cells(cells);
	std::array<double, 3> center;
	std::array<double, 3> razmer;

	double cp1 = 0.2;
	double cp2 = 0.8;
	double u1 = -0.3;
	double u2 = +0.5;

	for (auto& i : cells)
	{
		i->Get_Center(this->AMR_self, center, razmer);
		double a, b, c, d, e, f;
		a = center[0] - razmer[0]/2;
		b = center[0] + razmer[0]/2;
		c = center[1] - razmer[1]/2;
		d = center[1] + razmer[1]/2;
		e = center[2] - razmer[2] / 2;
		f = center[2] + razmer[2] / 2;
		double s1 = maxwell(1.0, cp1, u1, 0.0, 0.0, a, c, e) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, a, c, e);
		double s2 = maxwell(1.0, cp1, u1, 0.0, 0.0, a, c, f) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, a, c, f);
		double s3 = maxwell(1.0, cp1, u1, 0.0, 0.0, a, d, e) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, a, d, e);
		double s4 = maxwell(1.0, cp1, u1, 0.0, 0.0, a, d, f) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, a, d, f);
		double s5 = maxwell(1.0, cp1, u1, 0.0, 0.0, b, c, e) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, b, c, e);
		double s6 = maxwell(1.0, cp1, u1, 0.0, 0.0, b, c, f) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, b, c, f);
		double s7 = maxwell(1.0, cp1, u1, 0.0, 0.0, b, d, e) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, b, d, e);
		double s8 = maxwell(1.0, cp1, u1, 0.0, 0.0, b, d, f) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, b, d, f);
		double s9 = maxwell(1.0, cp1, u1, 0.0, 0.0, (a + b)/2, (c + d)/2, (e + f)/2) +
			10 * maxwell(1.0, cp2, u2, 0.0, 0.0, (a + b) / 2, (c + d) / 2, (e + f) / 2);

		i->f = s9;
	}
}

unsigned int AMR_f::Refine(void)
{
	this->Sf = 0.0;
	this->Sfu = 0.0;
	this->Sfuu = 0.0;

	std::vector<AMR_cell*> cells;
	std::array<double, 3> center;
	std::array<double, 3> razmer;

	this->Get_all_cells(cells);
	double V, u, m, mu, muu;

	//cout << "All_cells_do = " << cells.size() << endl;

	for (const auto& i : cells)
	{
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		u = norm2(center[0], center[1], center[2]);
		this->Sf += V * i->f;
		this->Sfu += V * i->f * u;
		this->Sfuu += V * i->f * kv(u);
	}

	if (this->Sf < 1e-8 || this->Sfu < 1e-8 || this->Sfuu < 1e-8) return 0;

	double procent = this->procent_signif;
	for (const auto& i : cells)
	{
		i->is_signif = false;
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		u = norm2(center[0], center[1], center[2]);
		m = V * i->f;
		mu = V * i->f * u;
		muu = V * i->f * kv(u);

		if (m * 100.0 / this->Sf > procent) i->is_signif = true;
		if (mu * 100.0 / this->Sfu > procent) i->is_signif = true;
		if (muu * 100.0 / this->Sfuu > procent) i->is_signif = true;
	}

	procent = this->procent_devide;
	for (const auto& i : cells)
	{
		i->need_devide_x = false;
		i->need_devide_y = false;
		i->need_devide_z = false;

		if (i->is_signif == false) continue;

		auto A = i->get_sosed(this->AMR_self, 0);
		if (A != nullptr) if (fabs(i->f - A->f) * 100.0 / i->f > procent)
		{
			i->need_devide_x = true;
			A->need_devide_x = true;
		}
		A = i->get_sosed(this->AMR_self, 1);
		if (A != nullptr) if (fabs(i->f - A->f) * 100.0 / i->f > procent)
		{
			i->need_devide_x = true;
			A->need_devide_x = true;
		}
		

		A = i->get_sosed(this->AMR_self, 2);
		if (A != nullptr) if (fabs(i->f - A->f) * 100.0 / i->f > procent)
		{
			i->need_devide_y = true;
			A->need_devide_y = true;
		}
		A = i->get_sosed(this->AMR_self, 3);
		if (A != nullptr) if (fabs(i->f - A->f) * 100.0 / i->f > procent)
		{
			i->need_devide_y = true;
			A->need_devide_y = true;
		}

		A = i->get_sosed(this->AMR_self, 4);
		if (A != nullptr) if (fabs(i->f - A->f) * 100.0 / i->f > procent)
		{
			i->need_devide_z = true;
			A->need_devide_z = true;
		}
		A = i->get_sosed(this->AMR_self, 5);
		if (A != nullptr) if (fabs(i->f - A->f) * 100.0 / i->f > procent)
		{
			i->need_devide_z = true;
			A->need_devide_z = true;
		}
	}

	short int k1 = 1;
	short int k2 = 1;
	short int k3 = 1;
	unsigned int NN = 0;
	for (const auto& i : cells)
	{
		if (i->need_devide_x == true)
		{
			k1 = 2;
		}
		else
		{
			k1 = 1;
		}

		if (i->need_devide_y == true)
		{
			k2 = 2;
		}
		else
		{
			k2 = 1;
		}

		if (i->need_devide_z == true)
		{
			k3 = 2;
		}
		else
		{
			k3 = 1;
		}

		if (k1 == 1 && k2 == 1 && k3 == 1) continue;

		i->divide(k1, k2, k3);
		NN++;
	}

	//cout << "Devide " << NN << "  cells" << endl;
	return NN;
}

void AMR_f::Save(string namef)
{
	std::ofstream out(namef, std::ios::binary);
	if (!out) {
		throw std::runtime_error("Cannot open file for writing: " + namef);
	}

	// Сохраняем важные параметры сетки
	double a;
	a = this->xL;
	out.write(reinterpret_cast<const char*>(&a), sizeof(double));
	a = this->xR;
	out.write(reinterpret_cast<const char*>(&a), sizeof(double));
	a = this->yL;
	out.write(reinterpret_cast<const char*>(&a), sizeof(double));
	a = this->yR;
	out.write(reinterpret_cast<const char*>(&a), sizeof(double));
	a = this->zL;
	out.write(reinterpret_cast<const char*>(&a), sizeof(double));
	a = this->zR;
	out.write(reinterpret_cast<const char*>(&a), sizeof(double));

	out.write(reinterpret_cast<const char*>(this->Vn.data()), 3 * sizeof(double));
	out.write(reinterpret_cast<const char*>(this->Vt.data()), 3 * sizeof(double));
	out.write(reinterpret_cast<const char*>(this->Vm.data()), 3 * sizeof(double));


	// Сохраняем размеры основной сетки
	size_t dims[3] = { this->cells.shape()[0], this->cells.shape()[1], this->cells.shape()[2] };
	out.write(reinterpret_cast<const char*>(dims), 3 * sizeof(size_t));

	// Записываем все ячейки
	for (size_t i = 0; i < dims[0]; ++i) 
	{
		for (size_t j = 0; j < dims[1]; ++j) 
		{
			for (size_t k = 0; k < dims[2]; ++k) 
			{
				this->cells[i][j][k]->Save_cell(out);
			}
		}
	}

}

void AMR_f::Read(string namef)
{
	std::ifstream in(namef, std::ios::binary);
	if (!in) {
		throw std::runtime_error("Cannot open file for reading: " + namef);
	}

	// Сохраняем важные параметры сетки
	double a;
	in.read(reinterpret_cast<char*>(&a), sizeof(double));
	this->xL = a;
	in.read(reinterpret_cast<char*>(&a), sizeof(double));
	this->xR = a;
	in.read(reinterpret_cast<char*>(&a), sizeof(double));
	this->yL = a;
	in.read(reinterpret_cast<char*>(&a), sizeof(double));
	this->yR = a;
	in.read(reinterpret_cast<char*>(&a), sizeof(double));
	this->zL = a;
	in.read(reinterpret_cast<char*>(&a), sizeof(double));
	this->zR = a;

	in.read(reinterpret_cast<char*>(this->Vn.data()), 3 * sizeof(double));
	in.read(reinterpret_cast<char*>(this->Vt.data()), 3 * sizeof(double));
	in.read(reinterpret_cast<char*>(this->Vm.data()), 3 * sizeof(double));

	// Читаем размеры multi_array
	size_t dims[3];
	in.read(reinterpret_cast<char*>(dims), 3 * sizeof(size_t));

	// Выделяем память под корневую сетку
	this->cells.resize(boost::extents[dims[0]][dims[1]][dims[2]]);

	for (size_t i = 0; i < dims[0]; ++i) 
	{
		for (size_t j = 0; j < dims[1]; ++j) 
		{
			for (size_t k = 0; k < dims[2]; ++k) 
			{
				this->cells[i][j][k] = new AMR_cell();
				this->cells[i][j][k]->Read_cell(in);
			}
		}
	}
}

void AMR_f::Print_info(void)
{
	const auto& shape = this->cells.shape();
	const size_t nx = shape[0];
	const size_t ny = shape[1];
	const size_t nz = shape[2];
	for (size_t i = 0; i < nx; ++i) 
	{
		for (size_t j = 0; j < ny; ++j) 
		{
			for (size_t k = 0; k < nz; ++k) 
			{
				AMR_cell* cell = cells[i][j][k];
				if (cell != nullptr) {
					cell->Print_info();
				}
			}
		}
	}
}

void AMR_f::Print_all_center_Tecplot(AMR_f* AMR)
{
	ofstream fout;
	string name_f = "Tecplot_AMR_f_print_cell_center_3D.txt";
	//std::vector<std::array<double, 3>> centers;
	std::array<double, 3> center;
	std::vector< AMR_cell*> cells;

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z, f" << endl;

	this->Get_all_cells(cells);

	/*const size_t dim1 = this->cells.shape()[0];
	const size_t dim2 = this->cells.shape()[1];
	const size_t dim3 = this->cells.shape()[2];

	for (size_t i = 0; i < dim1; ++i) {
		for (size_t j = 0; j < dim2; ++j) {
			for (size_t k = 0; k < dim3; ++k) {
				AMR_cell* cell = cells[i][j][k];
				cell->Get_Centers(AMR, centers);
			}
		}
	}*/

	for (auto& j : cells)
	{
		j->Get_Center(this->AMR_self, center);
		fout << center[0] << " " << center[1] << " " << center[2] << " " << j->f << endl;
	}

	fout.close();
	cout << "Print " << cells.size() << " cells" << endl;
}

void AMR_f::Print_slice_Tecplot(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d)
{
	std::vector<AMR_cell*> all_cells;
	this->Get_all_cells(all_cells);
	std::vector < std::vector<std::array<double, 3>>> points;

	for (auto& i : all_cells)
	{
		i->Slice_plane(AMR, a, b, c, d, points);
	}

	unsigned int NN2 = 0;

	for (const auto& i : points)
	{
		NN2 += i.size();
	}

	ofstream fout;
	string name_f = "Tecplot_setka_srez.txt";

	fout.open(name_f);
	fout << "TITLE = HP" << endl;
	fout << "VARIABLES = X, Y, Z" << endl;
	fout << "ZONE T=HP, NODES = " << NN2 << ", ELEMENTS = " << NN2 << ", F = FEPOINT, ET = LINESEG" << endl;

	Eigen::Vector3d C;
	for (const auto& i : points)
	{
		for (const auto& j : i)
		{
			C(0) = j[0];
			C(1) = j[1];
			C(2) = j[2];
			fout << C(0) << " " << C(1) << " " << C(2) << endl;
		}
	}


	size_t all_k1 = 1;
	for (const auto& i : points)
	{
		size_t k1 = i.size();
		for (size_t ii = 0; ii < k1; ii++)
		{
			size_t k2 = ii + 1;
			if (k2 >= k1) k2 = 0;
			fout << all_k1 + ii << " " << all_k1 + k2 << endl;
		}

		all_k1 = all_k1 + k1;
	}

	fout.close();
}

void AMR_f::Print_all_sosed_Tecplot(AMR_f* AMR)
{
	ofstream fout;
	string name_f = "Tecplot_print_sosed_3D.txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	vector< AMR_cell*> all_cells;
	this->Get_all_cells(all_cells);

	fout << "ZONE T=HP, N = " << all_cells.size() * 6 * 2 << ", E = " << all_cells.size() * 6 << ", F=FEPOINT, ET=LINESEG" << endl;

	for (auto& cc : all_cells)
	{
		std::array<double, 3> center;
		cc->Get_Center(AMR, center);
		for (size_t ii = 0; ii < 6; ii++)
		{
			auto ss = cc->get_sosed(AMR, ii);
			if (ss == nullptr) ss = cc;

			std::array<double, 3> center2;
			ss->Get_Center(AMR, center2);

			fout << center[0] << " " << center[1] << " " << center[2] << endl;
			fout << (center[0] + center2[0]) / 2.0 << " " << (center[1] + center2[1]) / 2.0
				<< " " << (center[2] + center2[2]) / 2.0 << endl;
		}
	}

	for (int m = 0; m < all_cells.size() * 6; m++)
	{
		fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
	}


	fout.close();
}

void AMR_f::Print_1D_Tecplot(AMR_f* AMR, const double& VV)
{
	Eigen::Vector3d ex, ey, ez, ee;
	ex << this->Vn[0], this->Vn[1], this->Vn[2];
	ey << this->Vt[0], this->Vt[1], this->Vt[2];
	ez << this->Vm[0], this->Vm[1], this->Vm[2];
	std::array<double, 3> center;
	ofstream fout;
	string name_f = "Tecplot_Print_1D_Tecplot.txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, f, ff" << endl;

	double x = -20;
	double dx = 0.03;

	double cp1 = 1.0;
	double u1 = VV;

	while (x < 20.0)
	{
		x = x + dx;
		ee << x, 0.0, 0.0;

		auto A = this->find_cell(ee.dot(ex), ee.dot(ey), ee.dot(ez));
		double ff = maxwell(1.0, cp1, u1, 0.0, 0.0, x, 0.0, 0.0);
		if (A == nullptr)
		{
			fout << x << " " << 0.0 << " " << ff << endl;
		}
		else
		{
			fout << x << " " << A->f << " " << ff << endl;
		}
	}

	fout.close();

}

void AMR_f::Delete(void)
{
	size_t dims[3] = { this->cells.shape()[0], this->cells.shape()[1], this->cells.shape()[2] };
	
	// Записываем все ячейки
	for (size_t i = 0; i < dims[0]; ++i)
	{
		for (size_t j = 0; j < dims[1]; ++j)
		{
			for (size_t k = 0; k < dims[2]; ++k)
			{
				this->cells[i][j][k]->Delete();
			}
		}
	}
	for (size_t i = 0; i < dims[0]; ++i)
	{
		for (size_t j = 0; j < dims[1]; ++j)
		{
			for (size_t k = 0; k < dims[2]; ++k)
			{
				delete this->cells[i][j][k];
			}
		}
	}

	this->cells.resize(boost::extents[0][0][0]);
}

