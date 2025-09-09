#include "Setka.h"
#include "Header.h"
using namespace std;

void Setka::Intitial_read(void)
{
	string str;
	int num, a1, a2, a3, a4;
	double b, b1, b2;

	ifstream fout;
	fout.open("Sergey_100.k");

	if (fout.is_open() == false)
	{
		cout << "ERROR open  " << "Sergey.k" << endl;
		exit(-100);
	}

	for (int i = 0; i < 10; i++)
	{
		fout >> str;
	}

	for (int ii = 0; ii < 1085; ii++) // 502
	{
		for (int i = 0; i < 10; i++)
		{
			fout >> num;
			if (i == 2) a1 = num;
			if (i == 3) a2 = num;
			if (i == 4) a3 = num;
			if (i == 5) a4 = num;
		}

		auto C = new Cell(a1, a2, a3, a4);
		this->all_Cells.push_back(C);
	}

	fout >> str;
	for (int i = 0; i < 6; i++)
	{
		fout >> str;
	}

	for (int ii = 0; ii < 1136; ii++) // 533
	{
		for (int i = 0; i < 6; i++)
		{
			fout >> b;
			if (i == 0) num = int(b);
			if (i == 1) b1 = b;
			if (i == 2) b2 = b;
		}

		auto C = new point(b1, b2, 0.0);
		C->id = num;
		this->all_points.push_back(C);
	}

	cout << "Setka.cpp -> Intitiak_read: Soedinyem" << endl;
	for (auto& I : this->all_Cells)
	{
		for (auto& j : I->int_points)
		{
			auto K = this->Find_point_id(j);
			if (K != nullptr)
			{
				I->points.push_back(K);
				K->cells.push_back(I);
			}
		}
	}
	cout << "Setka.cpp -> Initial_read: END Soedinyem" << endl;

	fout.close();
}

void Setka::regularize(void)
{

	for (int step = 0; step < 8000; step++)
	{
		for (auto& I : this->all_points)
		{
			I->Vx = 0.0;
			I->Vy = 0.0;
			I->Vnum = 0;
		}

		for (auto& I : this->all_Cells)
		{
			int k, kp, km;
			for (k = 0; k < 4; k++)
			{
				km = k - 1;
				kp = k + 1;
				if (kp >= 4) kp = 0;
				if (km < 0) km = 3;
				double x0, y0, x1, y1, x2, y2;
				double u, v, u2, v2;
				double n1, n2;
				double uu, vv, skalar, l;
				x0 = I->points[k]->x;
				y0 = I->points[k]->y;
				x1 = I->points[kp]->x;
				y1 = I->points[kp]->y;
				x2 = I->points[km]->x;
				y2 = I->points[km]->y;
				if (I->points[k]->get_radius() > 0.999) continue;  // 40
				u = x1 - x0;
				v = y1 - y0;
				u2 = x2 - x0;
				v2 = y2 - y0;
				double ll1 = sqrt((u) * (u) + (v) * (v));
				double ll2 = sqrt((u2) * (u2) + (v2) * (v2));
				double stepen = 5.0;

				I->points[k]->Vx += (u * pow(ll1, stepen) + u2 * pow(ll2, stepen)) / 10.0 / pow((ll1 + ll2), stepen); //10
				I->points[k]->Vy += (v * pow(ll1, stepen) + v2 * pow(ll2, stepen)) / 10.0 / pow((ll1 + ll2), stepen); // 10
				I->points[k]->Vnum += 1;

				/*n1 = sqrt(u * u + v * v);
				n2 = sqrt(u2 * u2 + v2 * v2);
				l = (n1 + n2) / 2.0;
				if(l < 3) continue;
				u = u / n1;
				v = v / n1;
				u2 = u2 / n2;
				v2 = v2 / n2;
				uu = (u + u2) / 2.0;
				vv = (v + v2) / 2.0;
				n1 = sqrt(uu * uu + vv * vv);
				uu = uu / n1;
				vv = vv / n1;
				skalar = acos(u * u2 + v * v2);
				I->points[k]->Vx += uu * (pi / 2 - skalar) / 20.0;
				I->points[k]->Vy += vv * (pi / 2 - skalar) / 20.0;
				I->points[k]->Vnum += 1;*/
			}
		}

		for (auto& I : this->all_points)
		{
			if (I->Vnum > 0)
			{
				I->x += I->Vx / I->Vnum;
				I->y += I->Vy / I->Vnum;
			}
		}

		for (auto& I : this->all_points)
		{
			I->Vx = 0.0;
			I->Vy = 0.0;
			I->Vnum = 0;
		}

		for (auto& I : this->all_Cells)
		{
			int k, kp, km;
			for (k = 0; k < 4; k++)
			{
				km = k - 1;
				kp = k + 1;
				if (kp >= 4) kp = 0;
				if (km < 0) km = 3;
				double x0, y0, x1, y1, x2, y2;
				double u, v, u2, v2;
				double n1, n2;
				double uu, vv, skalar, l;
				x0 = I->points[k]->x;
				y0 = I->points[k]->y;
				x1 = I->points[kp]->x;
				y1 = I->points[kp]->y;
				x2 = I->points[km]->x;
				y2 = I->points[km]->y;
				if (I->points[k]->get_radius() > 0.999) continue; // 40.0
				u = x1 - x0;
				v = y1 - y0;
				u2 = x2 - x0;
				v2 = y2 - y0;

				n1 = sqrt(u * u + v * v);
				n2 = sqrt(u2 * u2 + v2 * v2);
				l = (n1 + n2) / 2.0;
				if(l < 3) continue;
				u = u / n1;
				v = v / n1;
				u2 = u2 / n2;
				v2 = v2 / n2;
				uu = (u + u2) / 2.0;
				vv = (v + v2) / 2.0;
				n1 = sqrt(uu * uu + vv * vv);
				uu = uu / n1;
				vv = vv / n1;
				skalar = acos(u * u2 + v * v2);
				I->points[k]->Vx += uu * pow((pi / 2 - skalar) * 3, 3) / 200.0; // 200
				I->points[k]->Vy += vv * pow((pi / 2 - skalar) * 3, 3) / 200.0; // 200
				I->points[k]->Vnum += 1;
			}
		}

		for (auto& I : this->all_points)
		{
			if (I->Vnum > 0)
			{
				I->x += I->Vx / I->Vnum;
				I->y += I->Vy / I->Vnum;
			}
		}


	}
}

void Setka::regularize2(void)
{
	for (int step = 0; step < 8000; step++) {
		double step_size = 0.00001;  // Уменьшаем шаг со временем

		for (auto& point : this->all_points) {
			point->Vx = point->Vy = 0.0;
			point->Vnum = 0;
		}

		for (auto& cell : this->all_Cells) {
			for (int k = 0; k < 4; k++) {
				if (cell->points[k]->get_radius() > 0.999) continue;

				int kp = (k + 1) % 4;
				int km = (k + 3) % 4;

				point* p0 = cell->points[k];
				point* p1 = cell->points[kp];
				point* p2 = cell->points[km];

				// Выравнивание длин
				double dx1 = p1->x - p0->x;
				double dy1 = p1->y - p0->y;
				double dx2 = p2->x - p0->x;
				double dy2 = p2->y - p0->y;

				double len1 = sqrt(dx1 * dx1 + dy1 * dy1);
				double len2 = sqrt(dx2 * dx2 + dy2 * dy2);

				p0->Vx += (dx1 / len1 + dx2 / len2) * step_size;
				p0->Vy += (dy1 / len1 + dy2 / len2) * step_size;

				// Выравнивание углов
				double dot = (dx1 * dx2 + dy1 * dy2) / (len1 * len2);
				dot = std::max(-1.0, std::min(1.0, dot));  // Ограничиваем для acos
				double angle = acos(dot);

				double target_angle = pi / 2.0;
				double angle_diff = target_angle - angle;

				// Биссектриса угла
				double bisect_x = (dx1 / len1 + dx2 / len2);
				double bisect_y = (dy1 / len1 + dy2 / len2);
				double bisect_len = sqrt(bisect_x * bisect_x + bisect_y * bisect_y);

				if (bisect_len > 1e-10) {
					//p0->Vx += bisect_x / bisect_len * angle_diff * step_size * 0.5;
					//p0->Vy += bisect_y / bisect_len * angle_diff * step_size * 0.5;
				}

				p0->Vnum += 1;
			}
		}

		// Обновление координат
		for (auto& point : this->all_points) {
			if (point->Vnum > 0 && point->get_radius() <= 0.999) {
				point->x += point->Vx / point->Vnum;
				point->y += point->Vy / point->Vnum;
			}
		}
	}
}


void Setka::Intitial_build(void)
{
	int N0, N1;
	double n, R;

	N0 = this->all_points.size();
	N1 = this->all_Cells.size();
	R = 0.0;

	for (auto& I : this->all_points)
	{
		n = sqrt(I->x * I->x + I->y * I->y + I->z * I->z);
		R = max(R, n);
	}

	for (auto& I : this->all_points)
	{
		n = sqrt(I->x * I->x + I->y * I->y + I->z * I->z);
		I->x /= (R * (1.0 + (R - n)/R/9));
		I->y /= (R * (1.0 + (R - n)/R/9));
		I->z /= (R * (1.0 + (R - n)/R/9));
		/*I->y /= R;
		I->z /= R;
		I->x /= R;*/
	}

	for (auto& I : this->all_points)
	{
		I->z = I->x;
		I->x = 1.0/tan(phi);
	}

	for (auto& I : this->all_points)
	{
		n = sqrt(I->x * I->x + I->y * I->y + I->z * I->z);
		I->x /= n;
		I->y /= n;
		I->z /= n;
	}

	// Теперь нужно создать такую же сетку с другой стороны оси Х

	for (int i = 0; i < N0; i++)
	{
		auto& I = this->all_points[i];
		auto C = new point(-I->x, I->y, I->z);
		C->id = -I->id;
		this->all_points.push_back(C);
	}

	for (int i = 0; i < N1; i++)
	{
		auto& I = this->all_Cells[i];
		auto C = new Cell(-I->int_points[0], -I->int_points[1], -I->int_points[2], -I->int_points[3]);
		this->all_Cells.push_back(C);
	}

	cout << "Setka.cpp -> Intitial_build: Soedinyem" << endl;
	for (int i = N1; i < this->all_Cells.size(); i++)
	{
		auto& I = this->all_Cells[i];
		for(auto& j : I->int_points)
		{
			auto K = this->Find_point_id(j);
			if (K != nullptr)
			{
				I->points.push_back(K);
				K->cells.push_back(I);
			}
		}
	}
	cout << "Setka.cpp -> Intitial_build: END Soedinyem" << endl;
	
}

void Setka::Intitial_build2(void)
{
	double dphi = 0.0;
	double dtheta = 0.0;
	double x, y, z, d, dmin;
	int N_do = 0;
	int i_min;

	N_do = all_points.size();

	PP = new point * *[Nphi];
	for (int i = 0; i < Nphi; i++)
	{
		PP[i] = new point *[Ntheta];
	}

	for (int i = 0; i < Nphi; i++)
	{
		dphi = (i + 0.5) * 2.0 * pi / Nphi;
		for (int j = 0; j < Ntheta; j++)
		{
			dtheta = pi - phi - (j + 1.0) * (pi - 2.0 * phi) / (Ntheta + 1.0);
			y = sin(dtheta) * cos(dphi);
			z = sin(dtheta) * sin(dphi);
			x = cos(dtheta);
			PP[i][j] = new point(x, y, z);
		}
	}

	for (int i = 0; i < Nphi - 1; i++)
	{
		for (int j = 0; j < Ntheta - 1; j++)
		{
			auto C = new Cell(PP[i][j], PP[i + 1][j], PP[i + 1][j + 1], PP[i][j + 1]);
			this->all_Cells.push_back(C);
			PP[i][j]->cells.push_back(C);
			PP[i + 1][j]->cells.push_back(C);
			PP[i + 1][j + 1]->cells.push_back(C);
			PP[i][j + 1]->cells.push_back(C);
		}
	}

	for (int j = 0; j < Ntheta - 1; j++)
	{
		auto C = new Cell(PP[Nphi - 1][j], PP[0][j], PP[0][j + 1], PP[Nphi - 1][j + 1]);
		this->all_Cells.push_back(C);
		PP[Nphi - 1][j]->cells.push_back(C);
		PP[0][j]->cells.push_back(C);
		PP[0][j + 1]->cells.push_back(C);
		PP[Nphi - 1][j + 1]->cells.push_back(C);
	}

	for (int i = 0; i < Nphi; i++)
	{
		for (int j = 0; j < Ntheta; j++)
		{
			this->all_points.push_back(PP[i][j]);
		}
	}

	// Соединяем круг с цилиндром (первая часть)
	for (int j = 0; j < Nphi; j++)
	{
		int jj = j + 1;
		if (jj == Nphi) jj = 0;
		int i_min2 = 0;

		dmin = 100000000.0;
		i_min = 0;

		for (int i = 0; i < N_do; i++)
		{
			auto K = this->all_points[i];
			d = Pdistance(K, PP[j][0]);
			if (d < dmin)
			{
				dmin = d;
				i_min = i;
			}
		}

		dmin = 100000000.0;
		for (int i = 0; i < N_do; i++)
		{
			auto K = this->all_points[i];
			d = Pdistance(K, PP[jj][0]);
			if (d < dmin)
			{
				dmin = d;
				i_min2 = i;
			}
		}


		auto C = new Cell(PP[j][0], this->all_points[i_min], this->all_points[i_min2], PP[jj][0]);
		this->all_Cells.push_back(C);
		PP[j][0]->cells.push_back(C);
		this->all_points[i_min]->cells.push_back(C);
		PP[jj][0]->cells.push_back(C);
		this->all_points[i_min2]->cells.push_back(C);
	}

	// Соединяем круг с цилиндром (вторая часть)
	for (int j = 0; j < Nphi; j++)
	{
		int jj = j + 1;
		if (jj == Nphi) jj = 0;
		int i_min2 = 0;

		dmin = 100000000.0;
		i_min = 0;

		for (int i = 0; i < N_do; i++)
		{
			auto K = this->all_points[i];
			d = Pdistance(K, PP[j][Ntheta - 1]);
			if (d < dmin)
			{
				dmin = d;
				i_min = i;
			}
		}

		dmin = 100000000.0;
		for (int i = 0; i < N_do; i++)
		{
			auto K = this->all_points[i];
			d = Pdistance(K, PP[jj][Ntheta - 1]);
			if (d < dmin)
			{
				dmin = d;
				i_min2 = i;
			}
		}


		auto C = new Cell(PP[j][Ntheta - 1], this->all_points[i_min], this->all_points[i_min2], PP[jj][Ntheta - 1]);
		this->all_Cells.push_back(C);
		PP[j][Ntheta - 1]->cells.push_back(C);
		this->all_points[i_min]->cells.push_back(C);
		PP[jj][Ntheta - 1]->cells.push_back(C);
		this->all_points[i_min2]->cells.push_back(C);
	}

}

void Setka::Intitial_build3(void)
{
	// Соединяем круг с цилиндром

}

void Setka::Find_cell_soseds(void)
{
	//  будем склеивать ячейки через узлы
	for (auto& i : this->all_points)
	{
		for (auto& j : i->cells)
		{
			if (j->soseds.size() == 4) continue;
			for (auto& jj : i->cells)
			{
				if (jj == j) continue;

				// Проверям, может быть ячейки уже являются найденными соседями
				if (find(j->soseds.begin(), j->soseds.end(), jj) != j->soseds.end()) continue;

				if (Cell_is_soseds(j, jj))
				{
					j->soseds.push_back(jj);
					jj->soseds.push_back(j);
				}
			}
		}
	}
}

bool Setka::Cell_is_soseds(Cell* A, Cell* B)
{
	int k = 0;
	for (auto& jj : B->points)
	{
		if (find(A->points.begin(), A->points.end(), jj) != A->points.end()) k++;
	}

	if (k == 2) return true;
	return false;
}

point* Setka::Find_point_id(const int& num)
{
	for (auto& I : this->all_points)
	{
		if (I->id == num) return I;
	}
	cout << "Setka.cpp -> Find_point_id: Net takogo yzla  NID:0000000001" << endl;
	return nullptr;
}

point* Setka::Get_point(const double& x, const double& y, const double& z, const double& r)
{
	for (auto& j : this->all_points)
	{
		double n1 = x - j->x;
		double n2 = y - j->y;
		double n3 = z - j->z;
		double rr = sqrt(n1 * n1 + n2 * n2 + n3 * n3);
		if (rr < r)
		{
			return j;
		}
	}
	return nullptr;
}

void Setka::Print_Setka_TecPlot(void)
{
	// Рисует сетку, основанную на ЛИНИЯХ, а не поверхностях
	ofstream fout;
	int N;

	N = this->all_Cells.size();

	fout.open("Setka_krug.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=  " << 4 * N;
	fout << " , E= " << 4 * N;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;

	for (auto& I : this->all_Cells)
	{
		for (auto& j : I->points)
		{
			fout << j->x << " " << j->y << " " << j->z << endl;
		}
	}

	for (int j = 0; j < N; j++)
	{
		fout << j * 4 + 1 << " " << j * 4 + 2 << endl;
		fout << j * 4 + 2 << " " << j * 4 + 3 << endl;
		fout << j * 4 + 3 << " " << j * 4 + 4 << endl;
		fout << j * 4 + 4 << " " << j * 4 + 1 << endl;
	}
}

void Setka::Print_Setka_TecPlot_surface(void)
{
	ofstream fout;
	int N;

	N = this->all_Cells.size();

	fout.open("Setka_surface_krug.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=  " << 4 * N;
	fout << " , E= " << N;
	fout << " , F=FEPOINT, ET=quadrilateral  " << endl;

	for (auto& I : this->all_Cells)
	{
		for (auto& j : I->points)
		{
			fout << j->x << " " << j->y << " " << j->z << endl;
		}
	}

	for (int j = 0; j < N; j++)
	{
		fout << j * 4 + 1 << " " << j * 4 + 2 << " " << j * 4 + 3 << " " << j * 4 + 4 << endl;
	}
}


void Setka::Print_cell_soseds()
{
	int ll = 0;
	int nn = 0;
	for (auto& i : this->all_Cells)
	{
		ll += i->soseds.size();
		i->Set_center();
	}
	nn = 2 * ll;

	ofstream fout;
	fout.open("Setka_all_soseds.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << nn;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;

	double x0, y0, x1, y1;

	for (auto& i : this->all_Cells)
	{
		x0 = i->xc;
		y0 = i->yc;

		for (auto& j : i->soseds)
		{
			x1 = j->xc;
			y1 = j->yc;
			fout << x0 << " " << y0 << endl;
			fout << (x0 + x1) / 2.0 << " " << (y0 + y1) / 2.0 << endl;
		}
	}

	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}


	fout.close();
}

void Setka::All_numerate()
{
	int k = 1;
	for (auto& i : this->all_Cells)
	{
		i->number = k;
		k++;
	}

	k = 1;
	for (auto& i : this->all_points)
	{
		i->number = k;
		k++;
	}
}


bool comp(point* a, point* b) {
	return polar_angle(a->x, a->y) < polar_angle(b->x, b->y);
}

void Setka::Save_for_3D(string name)
{
	this->All_numerate();

	vector<point*> obolochka;

	// Нужно сначала вывести все точки на границе круга, ещё и в правильном порядке

	double Rmax = 0.0;
	for (auto& i : this->all_points)
	{
		Rmax = max(Rmax, sqrt(kv(i->x) + kv(i->y)));
	}

	for (auto& i : this->all_points)
	{
		i->x = i->x / Rmax;
		i->y = i->y / Rmax;
		i->z = i->z / Rmax;
	}


	double r;
	for (auto& i : this->all_points)
	{
		r = sqrt(kv(i->x) + kv(i->y));
		if (r < Rmax - 0.01) continue;
		obolochka.push_back(i);
	}

	sort(obolochka.begin(), obolochka.end(), comp);

	/*for (auto& i : obolochka)
	{
		cout << i->x << " " << i->y << " " << polar_angle(i->x, i->y) << endl;
		system("pause");
	}*/



	ofstream out(name + "_krug_setka.bin", ios::binary | ios::out);
	int n;
	double x, y;

	// Записываем номер версии файла вывода (для того, чтобы можно было отслеживать версии)
	n = 1;
	out.write((char*)&n, sizeof n);


	n = this->all_points.size();
	out.write((char*)&n, sizeof n);
	
	for (auto& i : obolochka)
	{
		x = i->x;
		y = i->y;
		//out.write((char*)&x, sizeof x);
		//out.write((char*)&y, sizeof y);
	}

	for (auto& i : this->all_points)
	{
		//if (find(obolochka.begin(), obolochka.end(), i) != obolochka.end()) continue;
		x = i->x;
		y = i->y;
		out.write((char*)&x, sizeof x);
		out.write((char*)&y, sizeof y);
	}

	// Записываем количество ячеек и все ячейки (номера узлов, из которых состоит ячейка)
	n = this->all_Cells.size();
	out.write((char*)&n, sizeof n);
	for (auto& i : this->all_Cells)
	{
		n = i->points.size();
		out.write((char*)&n, sizeof n);
		for (auto& j : i->points)
		{
			n = j->number;
			out.write((char*)&n, sizeof n);
		}
	}

	// Теперь записываем всех соседей ячейки
	for (auto& i : this->all_Cells)
	{
		n = i->soseds.size();
		out.write((char*)&n, sizeof n);
		for (auto& j : i->soseds)
		{
			n = j->number;
			out.write((char*)&n, sizeof n);
		}
	}


	int kkk = 1000;
	// Для будующих добавок какой-либо информации 
	// Если добавки нет, пишется ноль, если есть, пишется 1 и сама добавка
	// Это позволяет впоследствии считывать файл, даже если чего-то нет в середине вывода
	// (например, где-то будут пикапы, где-то нет)

	// Если есть добавка пишем так:
	// {---------------------------------
	//kkk--;
	//n = 1;
	//out.write((char*)&n, sizeof n);
	// Тело добавки
	// }---------------------------------

	n = 0;
	for (int i = 0; i < kkk; i++)
	{
		out.write((char*)&n, sizeof n);
	}


	out.close();

}
