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

	//cout << "All_cells_do = " << cells.size() << endl;

	for (const auto& i : cells)
	{
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		//this->Get_real_koordinate(center[0], center[1], center[2],
		//	Vx, Vy, Vz);
		//i->Spotok = V * i->f * fabs(Vx);
		i->setSpotok(V * i->getF() * fabs(center[0]));
		this->SpotokV += i->getSpotok();
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
				cell->setSpotok(cell->Get_SpotokV());
			}
		}
	}



}

void AMR_f::Partially_free_space(void)
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
				cell->Cell_partially_free_space();
			}
		}
	}
}

void AMR_f::Re_Partially_free_space(void)
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
				cell->Re_Cell_partially_free_space();
			}
		}
	}
}

void AMR_f::Get_random_velosity(AMR_f* AMR, const double& Squ, Eigen::Vector3d& Vel, Sensor* Sens)
{
	// Squ - площадь грани
	// Vel - разыгрываемая скорость частицы
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
				if (SS + cell->getSpotok() * Squ > ksi)
				{
					// Эта ячейка нам подходит
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
					SS += cell->getSpotok() * Squ;
				}
			}
		}
	}

	cout << "Error 6453211877  " << endl;
	whach(SS);
	whach(nx);
	whach(ny);
	whach(nz);
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

	this->Sf = 0.0;
	this->Sfu = 0.0;
	this->Sfux = 0.0;
	this->Sfuu = 0.0;
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

	this->Sf = 0.0;
	this->Sfu = 0.0;
	this->Sfux = 0.0;
	this->Sfuu = 0.0;


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


void AMR_f::AMR_resize(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
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

	this->cells.resize(boost::extents[xn][yn][zn]);

	for (int i = 0; i < xn; ++i) 
	{
		for (int j = 0; j < yn; ++j) 
		{
			for (int k = 0; k < zn; ++k) 
			{
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


	this->mut.lock();
	this->parameters["n"] += (mu / Vnn);
	this->parameters["Smu"] += mu;
	this->parameters["nn"] += mu * kv(1.0 / Vnn);
	this->mut.unlock();

	if (Vnn < 1e-6) Vnn = 1e-6;
	this->Get_lokal_koordinate(Vx, Vy, Vz, x, y, z);

	if (x < this->xL) x = this->xL + 1e-7;
	if (x > this->xR) x = this->xR - 1e-7;

	if (y < this->yL) y = this->yL + 1e-7;
	if (y > this->yR) y = this->yR - 1e-7;

	if (z < this->zL) z = this->zL + 1e-7;
	if (z > this->zR) z = this->zR - 1e-7;

	unsigned short int kkk = 0;
	a2:
	kkk++;
	AMR_cell* cell = find_cell(x, y, z);

	if (kkk > 30)
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
		whach(this->xL);
		whach(this->xR);
		whach(this->yL);
		whach(this->yR);
		whach(this->zL);
		whach(this->zR);
		exit(-1);
	}

	if (cell == nullptr)
	{
		if (kkk < 20)
		{
			x *= 0.995;
			y *= 0.995;
			z *= 0.995;
		}
		else
		{
			x *= 0.995;
		}
		goto a2;
	}

	this->mut.lock();
	//cell->f += mu/Vnn;
	cell->setF(cell->getF() + mu);    // Убираем нормировку на Vnn, будем нормировать потоком на границе
	this->mut.unlock();
}

void AMR_f::Normir_velocity_volume(const double& squ)
{
	std::vector<AMR_cell*> cells;
	std::array<double, 3> center;
	std::array<double, 3> razmer;
	this->Get_all_cells(cells);

	if (this->parameters.find("n") != this->parameters.end()) this->parameters["n"] /= squ;
	if (this->parameters.find("Smu") != this->parameters.end()) this->parameters["Smu"] /= squ;
	if (this->parameters.find("nn") != this->parameters.end()) this->parameters["nn"] /= squ;

	for (auto& cel : cells)
	{
		cel->Get_Center(this->AMR_self, center, razmer);
		//cel->f /= (razmer[0] * razmer[1] * razmer[2] * squ);
		cel->setF(cel->getF() / (razmer[0] * razmer[1] * razmer[2] * squ * center[0])); 
		// нормировка еще и на элемент Vx (деления источника заново выведены)
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

	if (A->isDivided() == false)
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
				if (cell->isDivided() == false) {
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

double AMR_f::Integrate_Maxwell_V(const double& xL, const double& xR,
	const double& yL, const double& yR,
	const double& zL, const double& zR, const double& Vinf)
{
	Eigen::Vector3d ex, ey, ez, ee;
	ex << this->Vn[0], this->Vn[1], this->Vn[2];
	ey << this->Vt[0], this->Vt[1], this->Vt[2];
	ez << this->Vm[0], this->Vm[1], this->Vm[2];

	const short unsigned int Nx = 10;
	const short unsigned int Ny = 10;
	const short unsigned int Nz = 10;

	double dx = (xR - xL) / (Nx);
	double dy = (yR - yL) / (Ny);
	double dz = (zR - zL) / (Nz);

	double x, y, z, f;
	double S = 0.0;

	for (short unsigned int i = 0; i < Nx; ++i)
	{
		for (short unsigned int j = 0; j < Ny; ++j)
		{
			for (short unsigned int k = 0; k < Nz; ++k)
			{
				x = xL + i * dx + dx / 2.0;
				y = yL + j * dy + dy / 2.0;
				z = zL + k * dz + dz / 2.0;
				ee = x * ex + y * ey + z * ez;
				f = maxwell(1.0, 1.0, Vinf, 0.0, 0.0, ee[0], ee[1], ee[2]);
				S += f * x;
			}
		}
	}

	S *= (dx * dy * dz);
	return S;
}

void AMR_f::Fill_maxwel_inf(const double& Vinf)
{
	Eigen::Vector3d ex, ey, ez, ee;
	ex << this->Vn[0], this->Vn[1], this->Vn[2];
	ey << this->Vt[0], this->Vt[1], this->Vt[2];
	ez << this->Vm[0], this->Vm[1], this->Vm[2];

	unsigned int NN = 1;
	unsigned int Nall = 0;

	while(NN > 0)
	{
		NN = 0;
		std::vector<AMR_cell*> cells;
		std::array<double, 3> center;
		std::array<double, 3> razmer;
		this->Get_all_cells(cells);

		for (auto& cel : cells)
		{
			cel->Get_Center(this->AMR_self, center, razmer);
			ee = center[0] * ex + center[1] * ey + center[2] * ez;
			double S = Integrate_Maxwell_V(center[0] - razmer[0] / 2.0,
				center[0] + razmer[0] / 2.0, center[1] - razmer[1] / 2.0,
				center[1] + razmer[1] / 2.0, center[2] - razmer[2] / 2.0,
				center[2] + razmer[2] / 2.0, Vinf);
			cel->setF(S / center[0] / (razmer[0] * razmer[1] * razmer[2]));
			//cel->f = maxwell(1.0, 1.0, Vinf, 0.0, 0.0, ee[0], ee[1], ee[2]);
		}

		NN = this->Refine(0);
		Nall = cells.size();
	}

	//cout << "Nall = " << Nall << endl;
}

void AMR_f::Fill_null(void)
{
	std::vector<AMR_cell*> cells;
	this->Get_all_cells(cells);

	for (auto& i : cells)
	{
		i->setF(0.0);
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

		i->setF(s9);
	}
}

unsigned int AMR_f::de_Refine(short int H_n)
{
	// Сначала нужно для всех ячеек (даже разделённых) посчитать значения переменных
	size_t dim1 = this->cells.shape()[0];
	size_t dim2 = this->cells.shape()[1];
	size_t dim3 = this->cells.shape()[2];

	std::array<double, 3> center;
	std::array<double, 3> razmer;

	for (size_t i = 0; i < dim1; ++i)
	{
		for (size_t j = 0; j < dim2; ++j)
		{
			for (size_t k = 0; k < dim3; ++k)
			{
				AMR_cell* cell = this->cells[i][j][k];
				cell->Get_Center(this->AMR_self, center, razmer);
				double SS = 0.0;
				cell->Get_f(this->AMR_self, SS);
				cell->setF(SS / (razmer[0] * razmer[1] * razmer[2] * center[0]));
			}
		}
	}

	// -------------------

	this->Sf = 0.0;
	this->Sfu = 0.0;
	this->Sfux = 0.0;
	this->Sfuu = 0.0;

	std::vector<AMR_cell*> cells;
	std::vector<AMR_cell*> parents;

	this->Get_all_cells(cells);
	double V, u, m = 0.0, mu = 0.0, muu = 0.0, mux = 0.0;

	for (const auto& i : cells)
	{
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		u = norm2(center[0], center[1], center[2]);
		this->Sf += V * i->getF();
		this->Sfu += V * i->getF() * u;
		this->Sfux += V * i->getF() * center[0];
		this->Sfuu += V * i->getF() * kv(u);
	}

	if (this->Sf < 1e-8 || this->Sfu < 1e-8 ||
		this->Sfuu < 1e-8 || this->Sfux < 1e-8) return 0;

	AMR_cell* parent;
	unsigned int N_delete = 0;
	double procent;

	for (const auto& i : cells)
	{
		parent = i->parent;
		if (parent == nullptr) continue;
		parent->setNeedDevideX(false);
	}

	for (const auto& i : cells)
	{
		procent = this->procent_signif / 2.0;
		m = 0.0; 
		mu = 0.0;
		muu = 0.0;
		mux = 0.0;
		parent = i->parent;
		if (parent == nullptr) continue;
		if (parent->needDevideX() == true) continue;
		parent->setIsSignif(false);
		parent->Get_Moment(this->AMR_self, m, mu, mux, muu);
		if (m * 100.0 / this->Sf > procent) parent->setIsSignif(true);
		if (mu * 100.0 / this->Sfu > procent) parent->setIsSignif(true);
		if (mux * 100.0 / this->Sfux > procent) parent->setIsSignif(true);
		if (muu * 100.0 / this->Sfuu > procent) parent->setIsSignif(true);


		// Специальное разбиение вокруг проблемной точки
		parent->Get_Center(this->AMR_self, center, razmer);
		if (false && (H_n == 3 || H_n == 4))
		{
			double LL = center[0] - razmer[0] / 2.0;
			double RR = center[0] + razmer[0] / 2.0;
			if (center[0] < -5.0 && center[0] > -7.0 && razmer[0] > 0.3)
			{
				parent->setNeedDevideX(false);
			}
			else if (LL < -5.0 && LL > -7.0 && razmer[0] > 0.3)
			{
				parent->setNeedDevideX(false);
			}
			else if (RR < -5.0 && RR > -7.0 && razmer[0] > 0.3)
			{
				parent->setNeedDevideX(false);
			}
			continue;
		}

		// Если ячека не пустая, то её размеры не могут быть больше 0.5 (чтобы нормально отделить нулевые области от ненулевых
		if (parent->getF() > 0.0000000001)
		{
			double dd = 2.0;
			if (H_n == 4)  dd = 0.7;
			if (H_n == 3)  dd = 0.7;
			if (H_n == 2)  dd = 3.0;
			if (H_n == 1)  dd = 2.0;
			if (razmer[0] > dd || razmer[1] > dd || razmer[2] > dd)
			{
				parent->setNeedDevideX(false);
				continue;
			}
		}

		if (parent->isSignif() == false)
		{
			parent->setNeedDevideX(true);
			parents.push_back(parent);
			continue;
		}

		// Проверяем если она существенная, но дальше делится не будет
		procent = this->procent_devide/2.0;

		bool bkl = false;
		for (short int il = 0; il < 6; il++)
		{
			auto A = parent->get_sosed(this->AMR_self, il);
			if (A != nullptr) if (fabs(parent->getF() - A->getF()) * 100.0 / parent->getF() > procent)
			{
				bkl = true;
				break;
			}
		}

		if (bkl == true) continue;


		// Если дошли до сюда, то можно удалять дочерние ячейки
		parent->setNeedDevideX(true);
		parents.push_back(parent);
		continue;
	}

	std::sort(parents.begin(), parents.end(), [](AMR_cell* a, AMR_cell* b) {
		return a->level > b->level;  // ">" для сортировки по убыванию
		});


	N_delete += parents.size();

	for (const auto& i : parents)
	{
		if(i == nullptr) continue;
		if(i->isDivided() == false) continue;
		if (i->needDevideX() == false) continue;

		i->setNeedDevideX(false);
		// В этом случае можно удалять дочерние ячейки
		i->setIsDivided(false);
		//i->Get_Center(this->AMR_self, center, razmer);
		//cout << center[0] << " " << center[1] << " " << center[2] << endl;
		dim1 = i->cells.shape()[0];
		dim2 = i->cells.shape()[1];
		dim3 = i->cells.shape()[2];
		for (size_t ii = 0; ii < dim1; ++ii)
		{
			for (size_t j = 0; j < dim2; ++j)
			{
				for (size_t k = 0; k < dim3; ++k)
				{
					AMR_cell* cell = i->cells[ii][j][k];
					cell->Delete();
					delete cell;
					i->cells[ii][j][k] = nullptr;
				}
			}
		}
		i->cells.resize(boost::extents[0][0][0]);
	}
	
	//if(N_delete > 0) cout << "Ydaleno yacheek  = " << N_delete << endl;
	return N_delete;
}

void AMR_f::Copy_and_Refine(std::vector<Int_point*>& Cells, Delaunay* Delone)
{
	std::vector<AMR_cell*> cells; 
	std::array<double, 3> center;

	std::unordered_map < string, double> parameters;
	Cell_handle prev_cell = Cell_handle();
	Cell_handle next_cell;
	unsigned int K2 = 10000;
	unsigned int kk = 0;

	while (K2 > 1)
	{
		cells.clear();
		this->Get_all_cells(cells);
		kk++;
		for (const auto& i : cells)
		{
			i->Get_Center(this->AMR_self, center);
			bool bb = Get_param_amr(center[0], center[1], center[2], parameters, Cells, Delone, prev_cell, next_cell);
			if(bb == true) prev_cell = next_cell;
			i->setF(parameters["f"]);
		}
		K2 = this->Refine(0);
		if (kk > 100)
		{
			cout << "kk > 100  " << kk << " " << K2 << endl;
			break;
		}
	}

	K2 = 10000;
	kk = 0;
	while (K2 > 1)
	{
		cells.clear();
		this->Get_all_cells(cells);
		kk++;
		for (const auto& i : cells)
		{
			i->Get_Center(this->AMR_self, center);
			bool bb = Get_param_amr(center[0], center[1], center[2], parameters, Cells, Delone, prev_cell, next_cell);
			if (bb == true) prev_cell = next_cell;
			i->setF(parameters["f"]);
		}
		K2 = this->de_Refine();
		if (kk > 100)
		{
			cout << "kk > 100  " << kk << " " << K2 << endl;
			break;
		}
	}



}

void AMR_f::Culc_gradients(void)
{
	std::vector<AMR_cell*> cells;
	this->Get_all_cells(cells);
	for (const auto& i : cells)
	{
		i->Culc_gradients(this->AMR_self);
	}
}

void AMR_f::Clean_low()
{
	// Функция чистит значения в ячейках, котрые в 3 раза меньше несущественных значений
	this->Sf = 0.0;
	this->Sfu = 0.0;
	this->Sfux = 0.0;
	this->Sfuu = 0.0;

	std::vector<AMR_cell*> cells;
	std::array<double, 3> center;
	std::array<double, 3> razmer;

	this->Get_all_cells(cells);
	double V, u, m = 0.0, mu = 0.0, muu = 0.0, mux = 0.0;

	//cout << "All_cells_do = " << cells.size() << endl;

	for (const auto& i : cells)
	{
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		u = norm2(center[0], center[1], center[2]);
		this->Sf += V * i->getF();
		this->Sfu += V * i->getF() * u;
		this->Sfux += V * i->getF() * center[0];
		this->Sfuu += V * i->getF() * kv(u);
	}

	double procent = this->procent_signif/100.0;
	for (const auto& i : cells)
	{
		i->setIsSignif(false);
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		u = norm2(center[0], center[1], center[2]);
		m = V * i->getF();
		mu = V * i->getF() * u;
		mux = V * i->getF() * center[0];
		muu = V * i->getF() * kv(u);

		if (m * 100.0 / this->Sf < procent) i->setF(0.0);
		if (mu * 100.0 / this->Sfu < procent) i->setF(0.0);
		if (mux * 100.0 / this->Sfux < procent) i->setF(0.0);
		if (muu * 100.0 / this->Sfuu < procent) i->setF(0.0);
	}
}

unsigned int AMR_f::Refine(short int H_n)
{
	// H_n может понадобиться для особого мельчения разных сортов
	this->Sf = 0.0;
	this->Sfu = 0.0;
	this->Sfux = 0.0;
	this->Sfuu = 0.0;

	std::vector<AMR_cell*> cells;
	std::array<double, 3> center;
	std::array<double, 3> razmer;

	this->Get_all_cells(cells);
	double V, u, m = 0.0, mu = 0.0, muu = 0.0, mux = 0.0;

	//cout << "All_cells_do = " << cells.size() << endl;

	for (const auto& i : cells)
	{
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		u = norm2(center[0], center[1], center[2]);
		this->Sf += V * i->getF();
		this->Sfu += V * i->getF() * u;
		this->Sfux += V * i->getF() * center[0];
		this->Sfuu += V * i->getF() * kv(u);
	}

	if (this->Sf < 1e-8 || this->Sfu < 1e-8 || this->Sfuu < 1e-8 || this->Sfux < 1e-8) return 0;

	for (const auto& i : cells)
	{
		i->setNeedDevideX(false);
		i->setNeedDevideY(false);
		i->setNeedDevideZ(false);
	}

	double procent = this->procent_signif;
	for (const auto& i : cells)
	{
		i->setIsSignif(false);
		i->Get_Center(this->AMR_self, center, razmer);
		V = razmer[0] * razmer[1] * razmer[2];
		u = norm2(center[0], center[1], center[2]);
		m = V * i->getF();
		mu = V * i->getF() * u;
		mux = V * i->getF() * center[0];
		muu = V * i->getF() * kv(u);

		if (m * 100.0 / this->Sf > procent) i->setIsSignif(true);
		if (mu * 100.0 / this->Sfu > procent) i->setIsSignif(true);
		if (mux * 100.0 / this->Sfux > procent) i->setIsSignif(true);
		if (muu * 100.0 / this->Sfuu > procent) i->setIsSignif(true);

		// Специальное разбиение вокруг проблемной точки
		if (false)
		{
			double LL = center[0] - razmer[0] / 2.0;
			double RR = center[0] + razmer[0] / 2.0;
			if (center[0] < -5.0 && center[0] > -7.0 && razmer[0] > 0.3)
			{
				i->setNeedDevideX(true);
			}
			else if(LL < -5.0 && LL > -7.0 && razmer[0] > 0.3)
			{
				i->setNeedDevideX(true);
			}
			else if (RR < -5.0 && RR > -7.0 && razmer[0] > 0.3)
			{
				i->setNeedDevideX(true);
			}
		}

		// Если ячека не пустая, то её размеры не могут быть больне 0.5 (чтобы нормально отделить нулевые области от ненулевых
		if (i->getF() > 0.0000000001)
		{
			double dd = 3.0;
			if (H_n == 9)  dd = 1.5;
			if (H_n == 8)  dd = 1.5;
			if (H_n == 7)  dd = 3.0;
			if (H_n == 6)  dd = 4.5;
			if (H_n == 5)  dd = 4.5;
			if (H_n == 4)  dd = 0.7;
			if (H_n == 3)  dd = 0.7;
			if (H_n == 2)  dd = 3.0;
			if (H_n == 1)  dd = 2.0;
			if(razmer[0] > dd)  i->setNeedDevideX(true);
			if(razmer[1] > dd)  i->setNeedDevideY(true);
			if(razmer[2] > dd)  i->setNeedDevideZ(true);
		}

	}

	procent = this->procent_devide;
	for (const auto& i : cells)
	{
		/*i->flags.need_devide_x = false;
		i->flags.need_devide_y = false;
		i->flags.need_devide_z = false;*/

		if (i->isSignif() == false) continue;

		// Добавил, чтобы ячеки вблизи несущественных имели размер не больше 0.3
		// Это чтобы было хорошо видно границу несущественных

		auto A = i->get_sosed(this->AMR_self, 0);
		if (A != nullptr)
		{
			if (fabs(i->getF() - A->getF()) * 100.0 / i->getF() > procent)
			{
				i->setNeedDevideX(true);
				A->setNeedDevideX(true);
			}
			else if (false)//(A->flags.is_signif == false)
			{
				i->Get_Center(this->AMR_self, center, razmer);
				if (razmer[0] > 0.3)
				{
					i->setNeedDevideX(true);
				}
			}
		}

		A = i->get_sosed(this->AMR_self, 1);
		if (A != nullptr)
		{
			if (fabs(i->getF() - A->getF()) * 100.0 / i->getF() > procent)
			{
				i->setNeedDevideX(true);
				A->setNeedDevideX(true);
			}
			else if (false)//(A->flags.is_signif == false)
			{
				i->Get_Center(this->AMR_self, center, razmer);
				if (razmer[0] > 0.3)
				{
					i->setNeedDevideX(true);
				}
			}
		}
	
		A = i->get_sosed(this->AMR_self, 2);
		if (A != nullptr)
		{
			if (fabs(i->getF() - A->getF()) * 100.0 / i->getF() > procent)
			{
				i->setNeedDevideY(true);
				A->setNeedDevideY(true);
			}
			else if (false)//(A->flags.is_signif == false)
			{
				i->Get_Center(this->AMR_self, center, razmer);
				if (razmer[1] > 0.3)
				{
					i->setNeedDevideX(true);
				}
			}
		}

		A = i->get_sosed(this->AMR_self, 3);
		if (A != nullptr)
		{
			if (fabs(i->getF() - A->getF()) * 100.0 / i->getF() > procent)
			{
				i->setNeedDevideY(true);
				A->setNeedDevideY(true);
			}
			else if (false)//(A->flags.is_signif == false)
			{
				i->Get_Center(this->AMR_self, center, razmer);
				if (razmer[1] > 0.3)
				{
					i->setNeedDevideX(true);
				}
			}
		}

		A = i->get_sosed(this->AMR_self, 4);
		if (A != nullptr)
		{
			if (fabs(i->getF() - A->getF()) * 100.0 / i->getF() > procent)
			{
				i->setNeedDevideZ(true);
				A->setNeedDevideZ(true);
			}
			else if (false)//(A->flags.is_signif == false)
			{
				i->Get_Center(this->AMR_self, center, razmer);
				if (razmer[2] > 0.3)
				{
					i->setNeedDevideX(true);
				}
			}
		}

		A = i->get_sosed(this->AMR_self, 5);
		if (A != nullptr)
		{
			if (fabs(i->getF() - A->getF()) * 100.0 / i->getF() > procent)
			{
				i->setNeedDevideZ(true);
				A->setNeedDevideZ(true);
			}
			else if (false)//(A->flags.is_signif == false)
			{
				i->Get_Center(this->AMR_self, center, razmer);
				if (razmer[2] > 0.3)
				{
					i->setNeedDevideX(true);
				}
			}
		}


		// Слишком маленькие ячейки больше делить не надо
		i->Get_Center(this->AMR_self, center, razmer);
		if (razmer[0] < 0.05) i->setNeedDevideX(false);
		if (razmer[1] < 0.05) i->setNeedDevideY(false);
		if (razmer[2] < 0.05) i->setNeedDevideZ(false);
	}

	

	short int k1 = 1;
	short int k2 = 1;
	short int k3 = 1;
	unsigned int NN = 0;
	for (const auto& i : cells)
	{
		if (i->needDevideX() == true)
		{
			k1 = 2;
		}
		else
		{
			k1 = 1;
		}

		if (i->needDevideY() == true)
		{
			k2 = 2;
		}
		else
		{
			k2 = 1;
		}

		if (i->needDevideZ() == true)
		{
			k3 = 2;
		}
		else
		{
			k3 = 1;
		}

		if (k1 == 1 && k2 == 1 && k3 == 1) continue;

		i->divide(this->AMR_self, k1, k2, k3);
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

	out.close();
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


	if (dims[0] > 1000 || dims[1] > 1000 || dims[2] > 1000)
	{
		cout << "ERROR  ergergergefgrvergheg" << endl;
		cout << namef << endl;
		cout << dims[0] << " " << dims[1] << " " << dims[2] << endl;
		exit(-1);
	}



	// Выделяем память под корневую сетку
	this->cells.resize(boost::extents[dims[0]][dims[1]][dims[2]]);

	this->xn = dims[0];
	this->yn = dims[1];
	this->zn = dims[2];

	for (size_t i = 0; i < dims[0]; ++i) 
	{
		for (size_t j = 0; j < dims[1]; ++j) 
		{
			for (size_t k = 0; k < dims[2]; ++k) 
			{
				this->cells[i][j][k] = new AMR_cell();
				this->cells[i][j][k]->nx = i;
				this->cells[i][j][k]->ny = j;
				this->cells[i][j][k]->nz = k;
				this->cells[i][j][k]->parent = nullptr;
				this->cells[i][j][k]->I_self = this->cells[i][j][k];
				this->cells[i][j][k]->Read_cell(in);
			}
		}
	}

	in.close();
}

unsigned int AMR_f::Size(void)
{
	std::vector<AMR_cell*> cells;
	this->Get_all_cells(cells);
	return cells.size();
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

void AMR_f::Print_all_center_Tecplot(AMR_f* AMR, const string& name)
{
	ofstream fout;
	string name_f = name + "_Tecplot_AMR_f_print_cell_center_3D.txt";
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
		fout << center[0] << " " << center[1] << " " << center[2] << " " << j->getF() << endl;
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
			fout << x << " " << A->getF() << " " << ff << endl;
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

void AMR_f::Analyze_memory_usage(bool detailed_output)
{
	size_t total_cells = 0;
	size_t active_cells = 0;
	size_t divided_cells = 0;
	size_t leaf_cells = 0;
	size_t memory_base_data = 0;
	size_t memory_active_data = 0;
	size_t memory_child_arrays = 0;
	size_t memory_param_maps = 0;
	size_t param_map_entries = 0;
	int max_depth = 0;

	cout << "=== ANALIZ ISPOLZOVANIYA PAMYATI AMR_f ===" << endl;
	cout << "Nachinayu analiz struktury AMR..." << endl;

	// Анализируем корневые ячейки
	const auto& shape = this->cells.shape();
	const size_t nx = shape[0];
	const size_t ny = shape[1];
	const size_t nz = shape[2];

	// Память самого объекта AMR_f
	size_t amr_f_base_size = sizeof(AMR_f);
	size_t amr_f_containers_size = 0;
	
	// Оценка памяти контейнеров AMR_f
	amr_f_containers_size += estimate_multi_array_memory(this->cells);
	amr_f_containers_size += estimate_unordered_map_memory(this->parameters);
	
	for (size_t i = 0; i < nx; ++i)
	{
		for (size_t j = 0; j < ny; ++j)
		{
			for (size_t k = 0; k < nz; ++k)
			{
				if (this->cells[i][j][k] != nullptr) 
				{
					analyze_cell_memory_recursive(this->cells[i][j][k], total_cells, active_cells,
						divided_cells, leaf_cells, memory_base_data, memory_active_data,
						memory_child_arrays, memory_param_maps, param_map_entries, 0, max_depth);
				}
			}
		}
	}

	// Вывод результатов
	cout << "\n=== OBSHCHAYA STATISTIKA ===" << endl;
	cout << "Vsego yacheek: " << total_cells << endl;
	cout << "Aktivnyh yacheek (nerazdelyonnykh): " << leaf_cells << endl;
	cout << "Razdelyonnykh yacheek: " << divided_cells << endl;
	cout << "Yacheek s aktivnymi dannymi: " << active_cells << endl;
	cout << "Maksimalnaya glubina dereva: " << max_depth << endl;

	cout << "\n=== ISPOLZOVANIE PAMYATI ===" << endl;
	
	// Память AMR_f
	cout << "Bazovyy obyekt AMR_f: " << amr_f_base_size << " bayt" << endl;
	cout << "Konteynery AMR_f: " << amr_f_containers_size << " bayt" << endl;
	
	// Память ячеек
	cout << "Bazovye dannye yacheek: " << memory_base_data << " bayt" << endl;
	cout << "Aktivnye dannye yacheek: " << memory_active_data << " bayt" << endl;
	cout << "Massivy dochernikh yacheek: " << memory_child_arrays << " bayt" << endl;
	cout << "Konteynery param (unordered_map): " << memory_param_maps << " bayt" << endl;
	
	size_t total_cells_memory = memory_base_data + memory_active_data + memory_child_arrays + memory_param_maps;
	size_t total_memory = amr_f_base_size + amr_f_containers_size + total_cells_memory;
	
	cout << "\nItogo pamyati yacheek: " << total_cells_memory << " bayt (" << total_cells_memory / 1024.0 / 1024.0 << " MB)" << endl;
	cout << "OBSHCHAYA PAMYAT: " << total_memory << " bayt (" << total_memory / 1024.0 / 1024.0 << " MB)" << endl;

	cout << "\n=== ANALIZ EFFEKTIVNOSTI ===" << endl;
	if (total_cells > 0) {
		cout << "Sredniy razmer yacheyki: " << (double)total_cells_memory / total_cells << " bayt" << endl;
		cout << "Protsent aktivnykh yacheek: " << (double)active_cells / total_cells * 100.0 << "%" << endl;
		cout << "Protsent pamyati na aktivnye dannye: " << (double)memory_active_data / total_cells_memory * 100.0 << "%" << endl;
		cout << "Protsent pamyati na param konteynery: " << (double)memory_param_maps / total_cells_memory * 100.0 << "%" << endl;
	}

	if (param_map_entries > 0) {
		cout << "Srednee kolichestvo parametrov na yacheyku s aktivnymi dannymi: " << (double)param_map_entries / active_cells << endl;
		cout << "Sredniy razmer odnogo param konteyner: " << (double)memory_param_maps / active_cells << " bayt" << endl;
	}

	cout << "\n=== REKOMENDATSII PO OPTIMIZATSII ===" << endl;
	
	if (memory_param_maps > memory_active_data / 2) {
		cout << "!!! KRITICHNO: Konteynery param zanimayut slishkom mnogo pamyati!" << endl;
		cout << "   Rekomendatsiya: Zamenit unordered_map na strukturu s fiksirovannymi polyami" << endl;
		cout << "   Potentsialnaya ekonomiya: ~" << memory_param_maps * 0.8 / 1024.0 / 1024.0 << " MB" << endl;
	}

	if ((double)memory_child_arrays / total_cells_memory > 0.3) {
		cout << "!!! Massivy dochernikh yacheek zanimayut mnogo pamyati" << endl;
		cout << "   Rekomendatsiya: Optimizirovat boost::multi_array ili ispolzovat vector" << endl;
	}

	if ((double)active_cells / total_cells < 0.5) {
		cout << "!!! Nizkiy protsent aktivnykh yacheek" << endl;
		cout << "   Rekomendatsiya: Rassmotret bolee agressivnuyu ochistku neaktivnykh dannykh" << endl;
	}

	if (detailed_output && total_cells > 0) {
		cout << "\n=== DETALNAYA INFORMATSIYA ===" << endl;
		cout << "Razmer bazovoy struktury AMR_cell: " << sizeof(AMR_cell) << " bayt" << endl;
		cout << "Razmer ActiveCellData: " << sizeof(ActiveCellData) << " bayt" << endl;
		cout << "Razmer unique_ptr<ActiveCellData>: " << sizeof(std::unique_ptr<ActiveCellData>) << " bayt" << endl;
		cout << "Razmer boost::multi_array<AMR_cell*, 3>: " << sizeof(boost::multi_array<AMR_cell*, 3>) << " bayt" << endl;
	}

	cout << "\n=== ANALIZ ZAVERSHYON ===" << endl;
}

void AMR_f::analyze_cell_memory_recursive(AMR_cell* cell, size_t& total_cells, size_t& active_cells, 
	size_t& divided_cells, size_t& leaf_cells, size_t& memory_base_data, 
	size_t& memory_active_data, size_t& memory_child_arrays, size_t& memory_param_maps, 
	size_t& param_map_entries, int depth, int& max_depth)
{
	if (cell == nullptr) return;

	total_cells++;
	max_depth = std::max(max_depth, depth);

	// Базовый размер структуры AMR_cell
	memory_base_data += sizeof(AMR_cell);

	// Проверяем, есть ли активные данные
	if (cell->has_active_data()) {
		active_cells++;
		memory_active_data += sizeof(ActiveCellData);
		
		// Анализируем param контейнер
		const auto& param_map = cell->getParam();
		param_map_entries += param_map.size();
		memory_param_maps += estimate_unordered_map_memory(param_map);
	}

	// Анализируем дочерние ячейки
	if (cell->isDivided()) {
		divided_cells++;
		
		// Память на multi_array дочерних ячеек
		const auto& child_cells = cell->cells;
		memory_child_arrays += estimate_multi_array_memory(child_cells);
		
		// Рекурсивно анализируем дочерние ячейки
		const auto& shape = child_cells.shape();
		for (size_t i = 0; i < shape[0]; ++i) {
			for (size_t j = 0; j < shape[1]; ++j) {
				for (size_t k = 0; k < shape[2]; ++k) {
					analyze_cell_memory_recursive(child_cells[i][j][k], total_cells, active_cells,
						divided_cells, leaf_cells, memory_base_data, memory_active_data,
						memory_child_arrays, memory_param_maps, param_map_entries, depth + 1, max_depth);
				}
			}
		}
	} else {
		leaf_cells++;
	}
}

size_t AMR_f::estimate_unordered_map_memory(const unordered_map<string, double>& map)
{
	size_t total_size = sizeof(unordered_map<string, double>);
	
	// Оценка памяти для bucket array
	total_size += map.bucket_count() * sizeof(void*);
	
	// Оценка памяти для элементов
	for (const auto& pair : map) {
		// Размер узла hash table (примерно)
		total_size += sizeof(void*) * 3; // next, hash, bucket pointers
		total_size += sizeof(pair);
		
		// Память строки (с учётом SSO - Small String Optimization)
		if (pair.first.size() > 15) { // Типичный порог SSO
			total_size += pair.first.capacity();
		}
	}
	
	return total_size;
}

size_t AMR_f::estimate_multi_array_memory(const boost::multi_array<AMR_cell*, 3>& arr)
{
	size_t total_size = sizeof(boost::multi_array<AMR_cell*, 3>);
	
	// Память для данных массива
	const auto& shape = arr.shape();
	size_t elements = shape[0] * shape[1] * shape[2];
	total_size += elements * sizeof(AMR_cell*);
	
	// Дополнительные метаданные boost::multi_array
	total_size += sizeof(shape) + 64; // примерная оценка overhead
	
	return total_size;
}

