#include "AMR_cell.h"

AMR_cell::AMR_cell()
{
	this->level = 0;
	this->flags.is_divided = false;
	// active_data создается при инициализации.
	this->ensure_active_data();
}

void AMR_cell::ensure_active_data()
{
	if (!active_data) 
	{
		active_data = std::make_unique<ActiveCellData>();
	}
}

void AMR_cell::clear_active_data()
{
	active_data.reset();
}

void AMR_cell::Get_Moment(AMR_f* AMR, double& m, double& mu, double& mux, double& muu)
{
	if (this->flags.is_divided == false)
	{
		std::array<double, 3> center;
		std::array<double, 3> razmer;
		this->Get_Center(AMR, center, razmer);
		double V = razmer[0] * razmer[1] * razmer[2];
		double u = norm2(center[0], center[1], center[2]);
		m += V * this->getF();
		mu += V * this->getF() * u;
		mux += V * this->getF() * center[0];
		muu += V * this->getF() * kv(u);
		return;
	}
	else
	{
		double SS = 0.0;
		const size_t dim1 = this->cells.shape()[0];
		const size_t dim2 = this->cells.shape()[1];
		const size_t dim3 = this->cells.shape()[2];

		for (size_t i = 0; i < dim1; ++i)
		{
			for (size_t j = 0; j < dim2; ++j)
			{
				for (size_t k = 0; k < dim3; ++k)
				{
					AMR_cell* cell = cells[i][j][k];
					cell->Get_Moment(AMR, m, mu, mux, muu);
				}
			}
		}
	}

	return;
}

void AMR_cell::Get_f(AMR_f* AMR, double& S)
{
	std::array<double, 3> center;
	std::array<double, 3> razmer;
	this->Get_Center(AMR, center, razmer);
	double V = razmer[0] * razmer[1] * razmer[2];

	if (this->flags.is_divided == false)
	{
		S += V * this->getF() * center[0];
		return;
	}
	else
	{
		this->setF(0.0);
		double SS = 0.0;
		const size_t dim1 = this->cells.shape()[0];
		const size_t dim2 = this->cells.shape()[1];
		const size_t dim3 = this->cells.shape()[2];

		for (size_t i = 0; i < dim1; ++i)
		{
			for (size_t j = 0; j < dim2; ++j)
			{
				for (size_t k = 0; k < dim3; ++k)
				{
					AMR_cell* cell = cells[i][j][k];
					double SS = 0.0;
					cell->Get_f(AMR, SS);
					this->setF(this->getF() + SS);
				}
			}
		}
		S += this->getF();
		this->setF(this->getF() / (V * center[0]));
	}

	return;
}

void AMR_cell::Cell_partially_free_space(void)
{
	if (this->flags.is_divided == false)
	{
		return;
	}
	else
	{
		this->clear_active_data();

		const size_t dim1 = this->cells.shape()[0];
		const size_t dim2 = this->cells.shape()[1];
		const size_t dim3 = this->cells.shape()[2];

		for (size_t i = 0; i < dim1; ++i)
		{
			for (size_t j = 0; j < dim2; ++j)
			{
				for (size_t k = 0; k < dim3; ++k)
				{
					AMR_cell* cell = cells[i][j][k];
					cell->Cell_partially_free_space();
				}
			}
		}
	}
}

void AMR_cell::Re_Cell_partially_free_space(void)
{
	if (this->flags.is_divided == false)
	{
		return;
	}
	else
	{
		this->ensure_active_data();

		const size_t dim1 = this->cells.shape()[0];
		const size_t dim2 = this->cells.shape()[1];
		const size_t dim3 = this->cells.shape()[2];

		for (size_t i = 0; i < dim1; ++i)
		{
			for (size_t j = 0; j < dim2; ++j)
			{
				for (size_t k = 0; k < dim3; ++k)
				{
					AMR_cell* cell = cells[i][j][k];
					cell->Re_Cell_partially_free_space();
				}
			}
		}
	}
}

double AMR_cell::Get_SpotokV(void)
{
	if (this->flags.is_divided == false)
	{
		return this->getSpotok();
	}
	else
	{
		double SS = 0.0;
		const size_t dim1 = this->cells.shape()[0];
		const size_t dim2 = this->cells.shape()[1];
		const size_t dim3 = this->cells.shape()[2];

		for (size_t i = 0; i < dim1; ++i) 
		{
			for (size_t j = 0; j < dim2; ++j) 
			{
				for (size_t k = 0; k < dim3; ++k) 
				{
					AMR_cell* cell = cells[i][j][k];
					SS += cell->Get_SpotokV();
				}
			}
		}
		this->setSpotok(SS);
		return SS;
	}

	return 0.0;
}

void AMR_cell::divide(AMR_f* AMR, unsigned short int n1, unsigned short int n2, unsigned short int n3)
{
	this->flags.is_divided = true;

	std::array<double, 3> center;
	std::array<double, 3> razmer;
	this->Get_Center(AMR, center, razmer);
	double x;

	this->cells.resize(boost::extents[n1][n2][n3]);
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k)
			{
				x = (center[0] - razmer[0] / 2.0) + razmer[0] / n1 * i + razmer[0] / n1 / 2.0;
				auto A = new AMR_cell();
				A->nx = i;
				A->ny = j;
				A->nz = k;
				A->parent = this->I_self;
				A->level = this->level + 1;
				this->cells[i][j][k] = A;
				A->I_self = A;

				A->setF(this->getF() * center[0] / x);                        //  Просто сносим значение
			}
		}
	}

	// После разделения можно очистить данные активной ячейки
	//this->clear_active_data();
}

AMR_cell* AMR_cell::find_cell(const double& x, const double& y, const double& z, const double& xL,
	const double& xR, const double& yL, const double& yR, const double& zL, const double& zR)
{
	unsigned int xn = this->cells.shape()[0];
	unsigned int yn = this->cells.shape()[1];
	unsigned int zn = this->cells.shape()[2];

	double dx = (xR - xL) / xn;
	int index1 = static_cast<int>((x - xL) / dx);
	if (index1 == xn) index1 = xn - 1;

	double dy = (yR - yL) / yn;
	int index2 = static_cast<int>((y - yL) / dy);
	if (index2 == yn) index2 = yn - 1;

	double dz = (zR - zL) / zn;
	int index3 = static_cast<int>((z - zL) / dz);
	if (index3 == zn) index3 = zn - 1;

	auto A = this->cells[index1][index2][index3];

	if (A->flags.is_divided == false)
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


void AMR_cell::Culc_gradients(AMR_f* AMR)
{
	// Обнулим градиенты
	this->getParam()["Bx"] = 0.0;
	this->getParam()["By"] = 0.0;
	this->getParam()["Bz"] = 0.0;

	std::array<double, 3> center;
	std::array<double, 3> razmer;
	std::array<double, 3> center1;
	std::array<double, 3> center2;

	this->Get_Center(AMR, center, razmer);

	AMR_cell* S1;
	AMR_cell* S2;
	double B1, B2;

	// !! Градианы только с одной стороны опасны тем, что функция может стать отрицательной, надо отдельно это проверять

	// Bx
	S1 = this->get_sosed(AMR, 0);
	S2 = this->get_sosed(AMR, 1);
	if (S1 != nullptr && S2 == nullptr)
	{
		S1->Get_Center(AMR, center1);
		this->getParam()["Bx"] = (S1->getF() - this->getF()) / (center1[0] - center[0]);
	}
	else if (S1 == nullptr && S2 != nullptr)
	{
		S2->Get_Center(AMR, center2);
		this->getParam()["Bx"] = (this->getF() - S2->getF()) / (center[0] - center2[0]);
	}
	else if (S1 != nullptr && S2 != nullptr)
	{
		S1->Get_Center(AMR, center1);
		S2->Get_Center(AMR, center2);

		B1 = (S1->getF() - this->getF()) / (center1[0] - center[0]);
		B2 = (this->getF() - S2->getF()) / (center[0] - center2[0]);
		this->getParam()["Bx"] = minmod(B1, B2);
	}

	// By
	S1 = this->get_sosed(AMR, 2);
	S2 = this->get_sosed(AMR, 3);
	if (S1 != nullptr && S2 == nullptr)
	{
		S1->Get_Center(AMR, center1);
		this->getParam()["By"] = (S1->getF() - this->getF()) / (center1[1] - center[1]);
	}
	else if (S1 == nullptr && S2 != nullptr)
	{
		S2->Get_Center(AMR, center2);
		this->getParam()["By"] = (this->getF() - S2->getF()) / (center[1] - center2[1]);
	}
	else if (S1 != nullptr && S2 != nullptr)
	{
		S1->Get_Center(AMR, center1);
		S2->Get_Center(AMR, center2);

		B1 = (S1->getF()- this->getF()) / (center1[1] - center[1]);
		B2 = (this->getF() - S2->getF()) / (center[1] - center2[1]);
		this->getParam()["By"] = minmod(B1, B2);
	}

	// Bz
	S1 = this->get_sosed(AMR, 4);
	S2 = this->get_sosed(AMR, 5);
	if (S1 != nullptr && S2 == nullptr)
	{
		S1->Get_Center(AMR, center1);
		this->getParam()["Bz"] = (S1->getF() - this->getF()) / (center1[2] - center[2]);
	}
	else if (S1 == nullptr && S2 != nullptr)
	{
		S2->Get_Center(AMR, center2);
		this->getParam()["Bz"] = (this->getF() - S2->getF()) / (center[2] - center2[2]);
	}
	else if (S1 != nullptr && S2 != nullptr)
	{
		S1->Get_Center(AMR, center1);
		S2->Get_Center(AMR, center2);

		B1 = (S1->getF() - this->getF()) / (center1[2] - center[2]);
		B2 = (this->getF() - S2->getF()) / (center[2] - center2[2]);
		this->getParam()["Bz"] = minmod(B1, B2);
	}

	// Проверка на неотрицательность функции
	// + читаем максимальное значения
	if (true)
	{
		double x, y, z;
		double S = -1.0;
		double SS = 0.0;      // Максимальное значение
		unsigned int kk = 0;
		while (S < 0.0)
		{
			kk++;
			S = 1.0;
			SS = this->getF() * (center[0] + (razmer[0] / 2.0));
			for (short int i = -1; i <= 1; i += 2) 
			{
				for (short int j = -1; j <= 1; j += 2)
				{
					for (short int k = -1; k <= 1; k += 2)
					{
						x = center[0] + i * (razmer[0] / 2.0);
						y = center[1] + j * (razmer[1] / 2.0);
						z = center[2] + k * (razmer[2] / 2.0);
						
						double ff = (this->getF() + this->getParam().at("Bx") * (x - center[0])
							+ this->getParam().at("By") * (y - center[1]) + this->getParam().at("Bz") * (z - center[2])) * x;
						S = min(S, ff);
						SS = max(SS, ff);
					}
				}
			}

			// Проверяем потенциальный максимум функции
			if (this->getParam().at("Bx") < 0.00001)
			{
				x = (this->getParam().at("Bx") * center[0] - this->getF()) / (2.0 * this->getParam().at("Bx"));
				if (x >= center[0] - (razmer[0] / 2.0) && x <= center[0] + (razmer[0] / 2.0))
				{
					y = center[1];
					z = center[2];
					double ff = (this->getF() + this->getParam().at("Bx") * (x - center[0])
						+ this->getParam().at("By") * (y - center[1]) + this->getParam().at("Bz") * (z - center[2])) * x;
					S = min(S, ff);
					SS = max(SS, ff);
				}
			}



			if (S < 0.0)
			{
				this->getParam()["Bx"] /= 2.0;
				this->getParam()["By"] /= 2.0;
				this->getParam()["Bz"] /= 2.0;
			}

			if (kk > 10)
			{
				if (this->getF() < 0.0)
				{
					cout << "Error 839yt478geofrjewfre = " << this->getF() << endl;
				}
				this->getParam()["Bx"] = 0.0;
				this->getParam()["By"] = 0.0;
				this->getParam()["Bz"] = 0.0;
				SS = this->getF();
				break;
			}
		}
		this->getParam()["Max"] = SS;
	}
}



AMR_cell* AMR_cell::get_sosed(AMR_f* AMR, short int nn)
{
	std::array<double, 3> center;
	std::array<double, 3> razmer;
	this->Get_Center(AMR, center, razmer);


	switch (nn) {
	case 0:
		return AMR->find_cell(center[0] + razmer[0] / 2.0 + razmer[0] / 1000.0,
			center[1], center[2]);
		break;
	case 1:
		return AMR->find_cell(center[0] - razmer[0] / 2.0 - razmer[0] / 1000.0,
			center[1], center[2]);
		break;
	case 2:
		return AMR->find_cell(center[0], 
			center[1] + razmer[1] / 2.0 + razmer[1] / 1000.0, center[2]);
		break;
	case 3:
		return AMR->find_cell(center[0], 
			center[1] - razmer[1] / 2.0 - razmer[1] / 1000.0, center[2]);
		break;
	case 4:
		return AMR->find_cell(center[0], center[1], 
			center[2] + razmer[2] / 2.0 + razmer[2] / 1000.0);
		break;
	case 5:
		return AMR->find_cell(center[0], center[1], 
			center[2] - razmer[2] / 2.0 - razmer[2] / 1000.0);
		break;
	default:
		cout << "ERROR 874658767843659837459" << endl;
		exit(-1);
	}

	return nullptr;
}

void AMR_cell::Print_info(void)
{
	if (this->flags.is_divided == false)
	{
		cout << "_________________________________" << endl;
		cout << "level: " << this->level << endl;
		cout << "is devided?: " << this->flags.is_divided << endl;
		cout << "nx: " << this->nx << endl;
		cout << "ny: " << this->ny << endl;
		cout << "nz: " << this->nz << endl;
		cout << "_________________________________" << endl;
	}
	else
	{
		cout << "_________________________________" << endl;
		cout << "level: " << this->level << endl;
		cout << "is devided?: " << this->flags.is_divided << endl;
		cout << "nx: " << this->nx << endl;
		cout << "ny: " << this->ny << endl;
		cout << "nz: " << this->nz << endl;
		cout << "_______Spusk_______" << endl;

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
}

void AMR_cell::Get_random_velosity_in_cell(AMR_f* AMR, const double& ksi,
	const double& Squ, Eigen::Vector3d& Vel, Sensor* Sens)
{
	const size_t dim1 = this->cells.shape()[0];
	const size_t dim2 = this->cells.shape()[1];
	const size_t dim3 = this->cells.shape()[2];
	double SS = 0.0;

	if (this->flags.is_divided == true)
	{
		for (size_t i = 0; i < dim1; ++i)
		{
			for (size_t j = 0; j < dim2; ++j)
			{
				for (size_t k = 0; k < dim3; ++k)
				{
					AMR_cell* cell = cells[i][j][k];
					if (SS + cell->getSpotok() * Squ > ksi)
					{
						return cell->Get_random_velosity_in_cell(AMR, ksi - SS, Squ, Vel, Sens);
					}
					else
					{
						SS += cell->getSpotok() * Squ;
					}
				}
			}
		}
	}
	else
	{
		std::array<double, 3> center;
		std::array<double, 3> razmer;
		this->Get_Center(AMR, center, razmer);

		// Равномерный розыгрышь 
		if (false)
		{
			Vel[0] = (center[0] - razmer[0] / 2.0) +
				Sens->MakeRandom() * (razmer[0]);
			Vel[1] = (center[1] - razmer[1] / 2.0) +
				Sens->MakeRandom() * (razmer[1]);
			Vel[2] = (center[2] - razmer[2] / 2.0) +
				Sens->MakeRandom() * (razmer[2]);
		}
		// Розыгрышь с плотностью: Vx  -  точно нужен для ячеек вблизи нуля!
		else if (center[0] - razmer[0] / 2.0 <= 0.0001)
		{
			double L = center[0] - razmer[0] / 2.0;
			double R = center[0] + razmer[0] / 2.0;
			Vel[0] = sqrt(kv(L) + Sens->MakeRandom() * (kv(R) - kv(L)));

			Vel[1] = (center[1] - razmer[1] / 2.0) +
				Sens->MakeRandom() * (razmer[1]);
			Vel[2] = (center[2] - razmer[2] / 2.0) +
				Sens->MakeRandom() * (razmer[2]);
		}
		// Розыгрышь со вторым порядком методом отказов
		else if (true)
		{
			double x, y, z, ff;
			do
			{
				x = (center[0] - razmer[0] / 2.0) +
					Sens->MakeRandom() * (razmer[0]);
				y = (center[1] - razmer[1] / 2.0) +
					Sens->MakeRandom() * (razmer[1]);
				z = (center[2] - razmer[2] / 2.0) +
					Sens->MakeRandom() * (razmer[2]);
				ff = (this->getF() + this->getParam().at("Bx") * (x - center[0])
					+ this->getParam().at("By") * (y - center[1]) + this->getParam().at("Bz") * (z - center[2])) * x;
			} while (Sens->MakeRandom() * this->getParam().at("Max") > ff);

			Vel[0] = x;
			Vel[1] = y;
			Vel[2] = z;
		}

		return;
	}

	



	cout << "Error 0900061211 " << endl;
	whach(SS);
	whach(ksi);
	whach(this->getSpotok() * Squ);
	whach(this->getSpotok());
	whach(Squ);
	exit(-1);
}

void AMR_cell::Get_index(std::vector<std::array<unsigned int, 3>>& numbers)
{
	//numbers.resize(this->level + 1);
	numbers[this->level] = std::array<unsigned int, 3>{this->nx, this->ny, this->nz};
	if (this->level == 0) return;

	this->parent->Get_index(numbers);
}

void AMR_cell::Get_Center(AMR_f* AMR, std::array<double, 3>& center)
{
	std::vector<std::array<unsigned int, 3>> numbers;
	numbers.resize(this->level + 1);
	this->Get_index(numbers);

	center[2] = center[1] = center[0] = 0.0;

	unsigned int xn = AMR->xn;
	unsigned int yn = AMR->yn;
	unsigned int zn = AMR->zn;

	double xL = AMR->xL;
	double xR = AMR->xR;

	double yL = AMR->yL;
	double yR = AMR->yR;

	double zL = AMR->zL;
	double zR = AMR->zR;

	AMR_cell* cell = nullptr;
	short int stk = 0;

	for (auto& i : numbers)
	{
		stk++;
		double dx = (xR - xL) / xn;
		double dy = (yR - yL) / yn;
		double dz = (zR - zL) / zn;
		xR = xL + (i[0] + 1) * dx;
		xL = xL + i[0] * dx;
		yR = yL + (i[1] + 1) * dy;
		yL = yL + i[1] * dy;
		zR = zL + (i[2] + 1) * dz;
		zL = zL + i[2] * dz;
		

		if (stk == 1)
		{
			cell = AMR->cells[i[0]][i[1]][i[2]];
		}
		else
		{
			/*cout << "A " << endl;
			whach(i[0]);
			whach(i[1]);
			whach(i[2]);
			whach(xL);
			whach(xR);
			whach(dx);
			whach(xn);*/
			cell = cell->cells[i[0]][i[1]][i[2]];
		}

		if (cell->flags.is_divided == true)
		{
			xn = cell->cells.shape()[0];
			yn = cell->cells.shape()[1];
			zn = cell->cells.shape()[2];
		}
	}

	center[0] = (xR + xL) / 2.0;
	center[1] = (yR + yL) / 2.0;
	center[2] = (zR + zL) / 2.0;

	return;
}

void AMR_cell::Get_Center(AMR_f* AMR, std::array<double, 3>& center, std::array<double, 3>& razmer)
{
	std::vector<std::array<unsigned int, 3>> numbers;
	numbers.resize(this->level + 1);
	this->Get_index(numbers);

	center[2] = center[1] = center[0] = 0.0;

	unsigned int xn = AMR->xn;
	unsigned int yn = AMR->yn;
	unsigned int zn = AMR->zn;

	double xL = AMR->xL;
	double xR = AMR->xR;

	double yL = AMR->yL;
	double yR = AMR->yR;

	double zL = AMR->zL;
	double zR = AMR->zR;

	AMR_cell* cell = nullptr;
	short int stk = 0;

	for (auto& i : numbers)
	{
		stk++;
		double dx = (xR - xL) / xn;
		double dy = (yR - yL) / yn;
		double dz = (zR - zL) / zn;
		xR = xL + (i[0] + 1) * dx;
		xL = xL + i[0] * dx;
		yR = yL + (i[1] + 1) * dy;
		yL = yL + i[1] * dy;
		zR = zL + (i[2] + 1) * dz;
		zL = zL + i[2] * dz;


		if (stk == 1)
		{
			cell = AMR->cells[i[0]][i[1]][i[2]];
		}
		else
		{
			cell = cell->cells[i[0]][i[1]][i[2]];
		}

		if (cell->flags.is_divided == true)
		{
			xn = cell->cells.shape()[0];
			yn = cell->cells.shape()[1];
			zn = cell->cells.shape()[2];
		}
	}

	center[0] = (xR + xL) / 2.0;
	center[1] = (yR + yL) / 2.0;
	center[2] = (zR + zL) / 2.0;

	razmer[0] = xR - xL;
	razmer[1] = yR - yL;
	razmer[2] = zR - zL;

	return;
}

void AMR_cell::Get_Centers(AMR_f* AMR, std::vector<std::array<double, 3>>& centers)
{
	if (this->flags.is_divided == false)
	{
		std::array<double, 3> center;
		this->Get_Center(AMR, center);
		centers.push_back(center);
	}
	else
	{
		const size_t dim1 = this->cells.shape()[0];
		const size_t dim2 = this->cells.shape()[1];
		const size_t dim3 = this->cells.shape()[2];

		for (size_t i = 0; i < dim1; ++i) {
			for (size_t j = 0; j < dim2; ++j) {
				for (size_t k = 0; k < dim3; ++k) {
					AMR_cell* cell = cells[i][j][k];
					cell->Get_Centers(AMR, centers);
				}
			}
		}
	}
}

void AMR_cell::Get_all_cells(vector<AMR_cell*>& cells)
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
				if (cell->flags.is_divided == false) {
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

void AMR_cell::Slice_plane(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d, std::vector < std::vector<std::array<double, 3>>>& points)
{
	std::vector< std::array<double, 3> > all_point;
	std::vector< std::array<double, 3> > kyb_point;
	std::array<double, 3> P1;
	std::array<double, 3> P2;
	std::array<double, 3> outIntersection;

	std::array<double, 3> center;
	std::array<double, 3> razmer;
	this->Get_Center(AMR, center, razmer);

	bool aa;

	P1[0] = center[0] - razmer[0] / 2.0;
	P1[1] = center[1] - razmer[1] / 2.0;
	P1[2] = center[2] - razmer[2] / 2.0;
	kyb_point.push_back(P1);

	P1[0] = center[0] + razmer[0] / 2.0;
	P1[1] = center[1] - razmer[1] / 2.0;
	P1[2] = center[2] - razmer[2] / 2.0;
	kyb_point.push_back(P1);

	P1[0] = center[0] + razmer[0] / 2.0;
	P1[1] = center[1] + razmer[1] / 2.0;
	P1[2] = center[2] - razmer[2] / 2.0;
	kyb_point.push_back(P1);

	P1[0] = center[0] - razmer[0] / 2.0;
	P1[1] = center[1] + razmer[1] / 2.0;
	P1[2] = center[2] - razmer[2] / 2.0;
	kyb_point.push_back(P1);

	// ----------------------

	P1[0] = center[0] - razmer[0] / 2.0;
	P1[1] = center[1] - razmer[1] / 2.0;
	P1[2] = center[2] + razmer[2] / 2.0;
	kyb_point.push_back(P1);

	P1[0] = center[0] + razmer[0] / 2.0;
	P1[1] = center[1] - razmer[1] / 2.0;
	P1[2] = center[2] + razmer[2] / 2.0;
	kyb_point.push_back(P1);

	P1[0] = center[0] + razmer[0] / 2.0;
	P1[1] = center[1] + razmer[1] / 2.0;
	P1[2] = center[2] + razmer[2] / 2.0;
	kyb_point.push_back(P1);

	P1[0] = center[0] - razmer[0] / 2.0;
	P1[1] = center[1] + razmer[1] / 2.0;
	P1[2] = center[2] + razmer[2] / 2.0;
	kyb_point.push_back(P1);

	for (size_t i = 0; i < 8; i++)
	{
		for (size_t j = i + 1; j < 8; j++)
		{
			P1 = kyb_point[i];
			P2 = kyb_point[j];

			int kkk = 0;

			if (fabs(P1[0] - P2[0]) < 0.00000001) kkk++;
			if (fabs(P1[1] - P2[1]) < 0.00000001) kkk++;
			if (fabs(P1[2] - P2[2]) < 0.00000001) kkk++;

			if (kkk != 2) continue;

			aa = findIntersection(P1, P2, a, b, c, d, outIntersection);
			if (aa == true) all_point.push_back(outIntersection);
		}
	}

	// std::vector< std::array<double, 3> > all_point;
	// Сейчас тут хранятся все найденные точки, которые надо рассортировать по кругу

	if (all_point.size() < 3)
	{
		return;
	}

	// Находим нормаль к плоскости

	std::array<double, 3> normal;
	normal[0] = a;
	normal[1] = b;
	normal[2] = c;

	double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
	if (length > 0) 
	{
		normal[0] /= length;
		normal[1] /= length;
		normal[2] /= length;
	}

	std::array<double, 3> centroid = { 0.0, 0.0, 0.0 };
	for (const auto& p : all_point) {
		centroid[0] += p[0];
		centroid[1] += p[1];
		centroid[2] += p[2];
	}
	centroid[0] /= all_point.size();
	centroid[1] /= all_point.size();
	centroid[2] /= all_point.size();

	// Выбираем произвольное направление для сортировки (например, ось OX в плоскости)
	std::array<double, 3> reference_dir;
	if (std::abs(normal[0]) > 0.9) 
	{ // Если нормаль близка к OX, выбираем OY
		reference_dir = { 0.0, 1.0, 0.0 };
	}
	else 
	{
		reference_dir = { 1.0, 0.0, 0.0 };
	}

	// Находим вектор в плоскости, перпендикулярный нормали
	std::array<double, 3>  tangent_dir = {
		reference_dir[1] * normal[2] - reference_dir[2] * normal[1],
		reference_dir[2] * normal[0] - reference_dir[0] * normal[2],
		reference_dir[0] * normal[1] - reference_dir[1] * normal[0]
	};

	// Сортируем точки по углу относительно tangent_dir
	std::sort(all_point.begin(), all_point.end(), [centroid, &tangent_dir]
	(std::array<double, 3> aa, std::array<double, 3> bb) 
	{
		std::array<double, 3> vec_a = { aa[0] - centroid[0], aa[1] - centroid[1], aa[2] - centroid[2] };
		std::array<double, 3> vec_b = { bb[0] - centroid[0], bb[1] - centroid[1], bb[2] - centroid[2] };;

		// Угол между vec_a и tangent_dir
		double dot_a = vec_a[0] * tangent_dir[0] + vec_a[1] * tangent_dir[1] + vec_a[2] * tangent_dir[2];
		double cross_a =
			tangent_dir[1] * vec_a[2] - tangent_dir[2] * vec_a[1] -
			tangent_dir[0] * vec_a[2] + tangent_dir[2] * vec_a[0] +
			tangent_dir[0] * vec_a[1] - tangent_dir[1] * vec_a[0];

		double angle_a = std::atan2(cross_a, dot_a);

		// Угол между vec_b и tangent_dir
		double dot_b = vec_b[0] * tangent_dir[0] + vec_b[1] * tangent_dir[1] + vec_b[2] * tangent_dir[2];
		double cross_b =
			tangent_dir[1] * vec_b[2] - tangent_dir[2] * vec_b[1] -
			tangent_dir[0] * vec_b[2] + tangent_dir[2] * vec_b[0] +
			tangent_dir[0] * vec_b[1] - tangent_dir[1] * vec_b[0];

		double angle_b = std::atan2(cross_b, dot_b);

		return angle_a < angle_b;
	});


	points.push_back(all_point);
	return;
}

void AMR_cell::Save_cell(std::ofstream& out)
{
	double h = this->getF();
	if (this->flags.is_divided == true) h = 0.0;
	out.write(reinterpret_cast<const char*>(&h), sizeof(double));

	unsigned short int l = this->level;
	out.write(reinterpret_cast<const char*>(&l), sizeof(unsigned short int));

	size_t dims[3] = { this->cells.shape()[0], this->cells.shape()[1], this->cells.shape()[2] };
	size_t n;
	n = dims[0];
	out.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
	n = dims[1];
	out.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
	n = dims[2];
	out.write(reinterpret_cast<const char*>(&n), sizeof(size_t));

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

void AMR_cell::Read_cell(std::ifstream& in)
{
	double f_value;
	in.read(reinterpret_cast<char*>(&f_value), sizeof(double));
	if (f_value < 0.0) f_value = 0.0;
	this->setF(f_value);
	
	in.read(reinterpret_cast<char*>(&this->level), sizeof(unsigned short int));

	size_t dims[3];
	in.read(reinterpret_cast<char*>(dims), 3 * sizeof(size_t));

	if (dims[0] > 1000 || dims[1] > 1000 || dims[2] > 1000)
	{
		cout << "ERROR  ertygru685467ujtuy3terdg" << endl;
		cout << dims[0] << " " << dims[1] << " " << dims[2] << endl;
		exit(-1);
	}

	// Выделяем память под вложенный массив
	this->cells.resize(boost::extents[dims[0]][dims[1]][dims[2]]);

	if (dims[0] > 0 || dims[1] > 0 || dims[2] > 0) this->flags.is_divided = true;

	// Рекурсивно читаем дочерние ячейки
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
				this->cells[i][j][k]->parent = this->I_self;
				this->cells[i][j][k]->I_self = this->cells[i][j][k];
				this->cells[i][j][k]->Read_cell(in);
			}
		}
	}

}

void AMR_cell::Delete(void)
{
	size_t dims[3] = { this->cells.shape()[0], this->cells.shape()[1], this->cells.shape()[2] };
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
