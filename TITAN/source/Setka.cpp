#include "Setka.h"
#include <algorithm>

using namespace std;

Setka::Setka()
{
	this->Surf1 = nullptr;
	this->geo = new Geo_param();
	Luch::geo = this->geo;
	this->geo->Nphi = 60;     // &INIT&
	// ! ��� ����� ������� �� ����� ������������ ����� (������� ��� �������� �� ���� ��� �������������)

	this->All_name_luch["A_Luch"] = &this->A_Luch;
	this->All_name_luch["B_Luch"] = &this->B_Luch;
	this->All_name_luch["C_Luch"] = &this->C_Luch;
	this->All_name_luch["D_Luch"] = &this->D_Luch;
	this->All_name_luch["E_Luch"] = &this->E_Luch;
	this->All_name_luch["H_Luch"] = &this->H_Luch;
	this->All_name_luch["G_Luch"] = &this->G_Luch;
	this->name_luch.reserve(7);
	this->name_luch.push_back("A_Luch");
	this->name_luch.push_back("B_Luch");
	this->name_luch.push_back("C_Luch");
	this->name_luch.push_back("D_Luch");
	this->name_luch.push_back("E_Luch");
	this->name_luch.push_back("H_Luch");
	this->name_luch.push_back("G_Luch");
	this->New_initial();

}

void Setka::Renumerate(void)
{
	int kkk = 1;
	for (auto& i : this->All_Yzel)
	{
		i->number = kkk;
		kkk++;
	}

	kkk = 1;
	for (auto& i : this->All_Cell)
	{
		i->number = kkk;
		kkk++;
	}
}

void Setka::Set_luch_parametr()
{
	// ��������� ����� ����������� ��������� 
	// (��������, ���� the � phi ��� ���������� ����� A, B, C �����) - ��� ������� �������
	for (auto& i : this->A2_Luch)
	{
		auto A = i->Yzels[0];
		double x = A->coord[0][0];
		double y = A->coord[0][1];
		double z = A->coord[0][2];
		double r = std::sqrt(x * x + y * y + z * z);
		i->parameters["the"] = std::acos(x / r);
		i->parameters["phi"] = polar_angle(y, z);
	}

	for (auto& j : this->A_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->B_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->D_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->E_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->H_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& j : this->G_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}

	for (auto& i : this->C2_Luch)
	{
		auto A = i->Yzels[0];
		double x = A->coord[0][0];
		double y = A->coord[0][1];
		double z = A->coord[0][2];
		double r = std::sqrt(x * x + y * y + z * z);
		i->parameters["the"] = std::acos(x / r);
		i->parameters["phi"] = polar_angle(y, z);
	}

	for (auto& j : this->C_Luch)
	{
		for (auto& i : j)
		{
			auto A = i->Yzels[0];
			double x = A->coord[0][0];
			double y = A->coord[0][1];
			double z = A->coord[0][2];
			double r = std::sqrt(x * x + y * y + z * z);
			i->parameters["the"] = std::acos(x / r);
			i->parameters["phi"] = polar_angle(y, z);
		}
	}
}

void Setka::Read_old_surface(string name)
{
	cout << "Start Read_old_surface: " << name << endl;
	if (this->Surf1 != nullptr)
	{
		this->Surf1->~Surfaces();
		delete this->Surf1;
	}

	this->Surf1 = new Surfaces();
	this->Surf1->Read_old(name);
	cout << "End Read_old_surface: " << name << endl;
}

void Setka::Move_to_surf(Surfaces* Surf)
{
	double x, y, z, phi, the, r, rr;
	for (int st = 0; st < 1; st++)
	{
		for (auto& i : this->A_Luch)
		{
			for (auto& j : i)
			{
				x = j->Yzels_opor[1]->coord[0][0];
				y = j->Yzels_opor[1]->coord[0][1];
				z = j->Yzels_opor[1]->coord[0][2];
				r = sqrt(kvv(x, y, z));
				the = polar_angle(x, sqrt(kv(y) + kv(z)));
				phi = polar_angle(y, z);

				if (phi == 6.1784655520599268)
				{
					//cout << 'rr' << endl;
				}

				rr = Surf->Get_TS(phi, the);

				j->Yzels_opor[1]->coord[0][0] *= rr / r;
				j->Yzels_opor[1]->coord[0][1] *= rr / r;
				j->Yzels_opor[1]->coord[0][2] *= rr / r;

				if (r < 0.0001 || rr < 0.0001 || std::isnan(rr) || std::fpclassify(rr) == FP_SUBNORMAL)
				{
					cout << "0989898653   errjr" << endl;
				}
				/*if (r > rr)
				{
					j->Yzels_opor[1]->coord[0][0] *= 0.99;
					j->Yzels_opor[1]->coord[0][1] *= 0.99;
					j->Yzels_opor[1]->coord[0][2] *= 0.99;
				}
				else
				{
					j->Yzels_opor[1]->coord[0][0] *= 1.01;
					j->Yzels_opor[1]->coord[0][1] *= 1.01;
					j->Yzels_opor[1]->coord[0][2] *= 1.01;
				}*/
			}
		}
	}

	for (auto& i : this->All_Luch)
	{
		i->dvigenie(0);
	}

}

void Setka::New_initial()
{
	cout << "---START New_initial---" << endl;

	int n, i, j, n2, n3;
	double x, y;
	vector<Yzel*> Yz_;

	// ��������� ���� 2� �����
	ifstream ffin("SDK1_2D_Setka.bin", ios::binary | ios::in);
	if (!ffin)
	{
		cout << "Net takogo fajla (fajl 2D setki)" << endl;
		exit(-1);
	}

	ffin.read((char*)&n, sizeof n);  // ������ �����
	cout << "Versiya fajla 2D setki: " << n << endl;

	// ��������� �������������� ���������
	if (true)
	{
		ffin.read((char*)&x, sizeof x);  this->geo->tetta0 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->tetta1 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->tetta2 = x;

		ffin.read((char*)&n, sizeof n);  this->geo->N1 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->N2 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->N3 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->N4 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->N5 = n;

		ffin.read((char*)&n, sizeof n);  this->geo->M0 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M1 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M11 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M2 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M3 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->M4 = n;
		ffin.read((char*)&n, sizeof n);  this->geo->MF = n;

		ffin.read((char*)&x, sizeof x);  this->geo->R0 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R1 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R2 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R3 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R4 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->R5 = x;

		ffin.read((char*)&x, sizeof x);  this->geo->L6 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->L7 = x;

		ffin.read((char*)&x, sizeof x);  this->geo->dd1 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd2 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd3 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd4 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd5 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd6 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd7 = x;
		ffin.read((char*)&x, sizeof x);  this->geo->dd8 = x;
	}

	// ��������� ����
	int sdvig_yz = 0;
	if (true)
	{
		ffin.read((char*)&n, sizeof n); // ����� �����
		for (int i = 0; i < n; i++)
		{
			ffin.read((char*)&x, sizeof x);
			ffin.read((char*)&y, sizeof y);
			auto yzel = new Yzel(x, y, 0.0);  // ������ ����
			this->All_Yzel.push_back(yzel);
			Yz_.push_back(yzel);
			sdvig_yz++;
		}

		Yzel_2D.push_back(Yz_);  // � ���� ������ ���������� ����������� ���������. ������ �������� ������������, �.�. ����� ������ Yz_, ��� �� �������� �� Yzel_2D
		Yz_.clear();
	}

	// ��������� ����
	if (true)
	{
		vector<Luch*> Luch_;
		//for (auto& [name, vec_ptr] : this->All_name_luch)
		for(auto& name : this->name_luch)
		{
			auto vec_ptr = this->All_name_luch[name];
			cout << "Zopolnyaem  " << name << endl;
			ffin.read((char*)&n, sizeof n); // ����� �����
			for (i = 0; i < n; i++)
			{
				auto Luch_0 = new Luch();  // ������ ���
				All_Luch.push_back(Luch_0);
				Luch_0->type = name;
				Luch_.push_back(Luch_0);
				ffin.read((char*)&n2, sizeof n2); // ����� �����
				for (j = 0; j < n2; j++)
				{
					ffin.read((char*)&n3, sizeof n3); // ����� ����
					if (Yzel_2D[0].size() < static_cast<size_t>(n3))
					{
						cout << "86740576536 PROBLEM   " << Yzel_2D[0].size() << " " << n3;
					}
					Luch_0->Yzels.push_back(Yzel_2D[0][n3 - 1]);
				}

				ffin.read((char*)&n2, sizeof n2); // ����� ������� �����
				//cout << name << "  Opor yzel = " << n2 << endl;
				for (j = 0; j < n2; j++)
				{
					ffin.read((char*)&n3, sizeof n3); // ����� ����
					Luch_0->Yzels_opor.push_back(Yzel_2D[0][n3 - 1]);
				}
			}
			vec_ptr->push_back(Luch_);
			cout << "Dobavili luchey:  " << Luch_.size() << " " << n << endl;
			Luch_.clear();
		}
	}

	// ������ ���� ������� ����� � 3�, ���������� ����� ����� � ����� � ���

	// ����������� ����� �����
	if (true)
	{
		for (int iii1 = 1; iii1 < this->geo->Nphi; iii1++)
		{
			Yz_.clear();
			for (auto& i : Yzel_2D[0])
			{
				auto yzel = new Yzel(i->coord[0][0], i->coord[0][1], i->coord[0][2]);  // ������ ����
				Yz_.push_back(yzel);
				this->All_Yzel.push_back(yzel);
			}
			this->Yzel_2D.push_back(Yz_);

			// ������� ���������� (������������ ���������)
			double tet = iii1 * (2.0 * const_pi) / this->geo->Nphi;
			for (auto& i : Yz_)
			{
				double newy = cos(tet) * i->coord[0][1] - sin(tet) * i->coord[0][2];
				double newz = sin(tet) * i->coord[0][1] + cos(tet) * i->coord[0][2];
				i->coord[0][1] = newy;
				i->coord[0][2] = newz;
			}

			Yz_.clear();
		}
	}

	// ����������� ����� �����
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}

		vector<Luch*> Luch_;
		for (auto& [name, vec_ptr] : this->All_name_luch)
		{
			vec_ptr->resize(this->geo->Nphi);
			cout << "Yvelichivaim luchi = " << name << endl;
			for (int iii1 = 1; iii1 < this->geo->Nphi; iii1++)
			{
				//cout << "cds1 = " << (*vec_ptr).size() << endl;
				//cout << "cds2 = " << (*vec_ptr)[0].size() << endl;
				//cout << "cds3 = " << (*vec_ptr)[1].size() << endl;
				//cout << "cds2 = " << (*vec_ptr)[0][0]->Yzels.size() << endl;

				(*vec_ptr)[iii1].resize((*vec_ptr)[0].size());
				int kk = 0;
				for (auto& ii : (*vec_ptr)[0])
				{
					auto Luch_0 = new Luch();  // ������ ���
					All_Luch.push_back(Luch_0);
					Luch_0->type = name;

					(*vec_ptr)[iii1][kk] = Luch_0;
					kk++;

					for (auto& j : ii->Yzels)
					{
						Luch_0->Yzels.push_back(this->All_Yzel[j->number - 1 + iii1 * (sdvig_yz)]);
					}
					for (auto& j : ii->Yzels_opor)
					{
						Luch_0->Yzels_opor.push_back(this->All_Yzel[j->number - 1 + iii1 * (sdvig_yz)]);
					}
				}
			}
		}
	}

	// ��������� ������
	if (true)
	{
		this->Cell_2D.resize(this->geo->Nphi);

		ffin.read((char*)&n, sizeof n); // ����� �����
		for (int i = 0; i < this->geo->Nphi; i++)
		{
			this->Cell_2D[i].resize(n);
		}
		this->All_Cell.reserve(this->geo->Nphi* n * 2);

		for (int i = 0; i < n; i++) // ��������� ������ ������
		{
			ffin.read((char*)&n2, sizeof n2);
			auto cell_ = new Cell();  // ������ ������
			this->All_Cell.push_back(cell_);
			this->Cell_2D[0][i] = cell_;

			for (int j = 0; j < n2; j++)
			{
				ffin.read((char*)&n3, sizeof n3);
				cell_->yzels.push_back(Yzel_2D[0][n3 - 1]);
			}
		}

		for (int iii1 = 1; iii1 < this->geo->Nphi; iii1++)
		{
			for (int i = 0; i < n; i++) // ��������� ������ ������
			{
				n2 = this->Cell_2D[0][i]->yzels.size();
				auto cell_ = new Cell();  // ������ ������
				this->All_Cell.push_back(cell_);
				this->Cell_2D[iii1][i] = cell_;

				for (int j = 0; j < n2; j++)
				{
					n3 = this->Cell_2D[0][i]->yzels[j]->number;
					cell_->yzels.push_back(Yzel_2D[iii1][n3 - 1]);
				}
			}
		}
	}

	// ������ ���� � ������ ������ �������� ���� � �������� ���� �����
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}

		for (int i = 0; i < this->Cell_2D.size(); i++)
		{
			j = i + 1;
			if (j > this->Cell_2D.size() - 1) j = 0;
			for (int k = 0; k < this->Cell_2D[0].size(); k++)
			{
				for (int kk = 0; kk < 4; kk++)
				{
					int nn = this->Cell_2D[j][k]->yzels[kk]->number;
					this->Cell_2D[i][k]->yzels.push_back(All_Yzel[nn - 1]);
				}
			}
		}
	}

	// ��������� ���� �����
	ifstream fin("SDK1_krug_setka.bin", ios::binary | ios::in);
	if (!fin)
	{
		cout << "Net takogo fajla (fajl setki v krugu)" << endl;
		exit(-1);
	}

	fin.read((char*)&n, sizeof n);  // ������ �����
	cout << "Versiya fajla setki v krugu: " << n << endl;

	// ��������� ���� (��������� �� � ������ ���� �������� �����)
	if (true) 
	{
		this->Krug_Yzel.resize(this->A_Luch[0][0]->Yzels.size());
		this->Krug_Yzel_2.resize(this->C_Luch[0][0]->Yzels.size());

		fin.read((char*)&n, sizeof n); // ����� �����

		for (int i = 0; i < this->Krug_Yzel.size(); i++)
		{
			this->Krug_Yzel[i].resize(n);
		}

		for (int i = 0; i < this->Krug_Yzel_2.size(); i++)
		{
			this->Krug_Yzel_2[i].resize(n);
		}
		
		for (int i = 0; i < n; i++)
		{
			fin.read((char*)&x, sizeof x);
			fin.read((char*)&y, sizeof y);

			for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
			{
				auto yzel = new Yzel(x, y, 0.0);
				//this->All_Yzel.push_back(yzel);
				yzel->number = -1;
				this->Krug_Yzel[ii][i] = yzel;
			}

			for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
			{
				auto yzel = new Yzel(x, y, 0.0);
				yzel->number = -1;
				//this->All_Yzel.push_back(yzel);
				this->Krug_Yzel_2[ii][i] = yzel;
			}
		}

		double Rmax = 1.0;
		double thetamin = 4.0;
		// ������� ������ �����

		//Rmax = sqrt(kv(this->Krug_Yzel[0][0]->coord[0][0]) + kv(this->Krug_Yzel[0][0]->coord[0][1]));
		//cout << "R max: " << Rmax << endl;

		Yzel* AAA = nullptr;
		for (int i = 0; i < this->Krug_Yzel[0].size(); i++)
		{
			if (this->Krug_Yzel[0][i]->func_R(0) > 1.0 - 0.00001)
			{
				AAA = this->Krug_Yzel[0][i];
				break;
			}
		}

		thetamin = polar_angle(AAA->coord[0][0], AAA->coord[0][1]);
		cout << "Ugol vrashcheniya setki: " << thetamin * 180.0 / const_pi << endl;
		cout << "this->geo.tetta0: " << this->geo->tetta0 << endl;
		cout << "this->geo.tetta1: " << this->geo->tetta1 << endl;

		// ������� ���������� ����� 1 �� �����
		for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel[ii])
			{
				i->coord[0][0] /= Rmax;
				i->coord[0][1] /= Rmax;
				i->coord[0][2] = i->coord[0][1];
				i->coord[0][1] = i->coord[0][0];
				i->coord[0][0] = 1.0 / tan(this->geo->tetta0);// 1.0 / tan(this->geo.tetta0 - (const_pi / 2.0 - this->geo.tetta0) / (this->geo.N1 - 1));
				double r = sqrt(kv(i->coord[0][0]) + kv(i->coord[0][1]) + kv(i->coord[0][2]));
				i->coord[0][0] *= this->geo->R0 / r * (ii + 1);
				i->coord[0][1] *= this->geo->R0 / r * (ii + 1);
				i->coord[0][2] *= this->geo->R0 / r * (ii + 1);
			}
		}

		// ������� ���������� ��� ���������� � ������� ������ �����
		for (int ii = 0; ii < this->Krug_Yzel.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel[ii])
			{
				double newy = cos(-thetamin) * i->coord[0][1] - sin(-thetamin) * i->coord[0][2];
				double newz = sin(-thetamin) * i->coord[0][1] + cos(-thetamin) * i->coord[0][2];
				i->coord[0][1] = newy;
				i->coord[0][2] = newz;
			}
		}

		// ������� ���������� ����� 2 �� �����
		for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel_2[ii])
			{
				i->coord[0][0] /= Rmax;
				i->coord[0][1] /= Rmax;
				i->coord[0][2] = i->coord[0][1];
				i->coord[0][1] = i->coord[0][0];
				i->coord[0][0] = 1.0 / tan(this->geo->tetta1);// 1.0 / tan(this->geo.tetta0 - (const_pi / 2.0 - this->geo.tetta0) / (this->geo.N1 - 1));
				double r = sqrt(kv(i->coord[0][0]) + kv(i->coord[0][1]) + kv(i->coord[0][2]));
				i->coord[0][0] *= this->geo->R0 / r * (ii + 1);
				i->coord[0][1] *= this->geo->R0 / r * (ii + 1);
				i->coord[0][2] *= this->geo->R0 / r * (ii + 1);

			}
		}

		// ������� ���������� ��� ���������� � ������� ������ �����
		for (int ii = 0; ii < this->Krug_Yzel_2.size(); ii++)
		{
			for (auto& i : this->Krug_Yzel_2[ii])
			{
				double newy = cos(-thetamin) * i->coord[0][1] - sin(-thetamin) * i->coord[0][2];
				double newz = sin(-thetamin) * i->coord[0][1] + cos(-thetamin) * i->coord[0][2];
				i->coord[0][1] = newy;
				i->coord[0][2] = newz;
			}
		}

	}

	// ���������������� ��� ����
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}
	}

	// ���� ������� ������ ���� (������ ���������)
	// ��������� ���� ���� ���������
	// ��� ����� ���� �� �� ������, ���� �� ����� � ����� � ����� ������ � ���������� ������� (������� ���������)
	cout << "Yzlov do dobavleniya  " << this->All_Yzel.size() << endl;
	for (int i = 0; i < this->A_Luch.size(); i++)
	{
		auto A = this->A_Luch[i][0]->Yzels[0];
		for (int j = 0; j < this->Krug_Yzel[0].size(); j++)
		{
			auto B = this->Krug_Yzel[0][j];


			double d = sqrt(kv(A->coord[0][0] - B->coord[0][0]) + kv(A->coord[0][1] - B->coord[0][1]) +
				kv(A->coord[0][2] - B->coord[0][2]));
			if (d < this->geo->R0 / 1000.0)
			{
				//cout << "Nashol " << j << endl;

				for (int k = 0; k < this->Krug_Yzel.size(); k++)
				{
					delete this->Krug_Yzel[k][j];
					this->Krug_Yzel[k][j] = this->A_Luch[i][0]->Yzels[k];
				}

				break;
			}
		}
		
	}

	for (int i = 0; i < this->C_Luch.size(); i++)
	{
		auto A = this->C_Luch[i].back()->Yzels[0];
		for (int j = 0; j < this->Krug_Yzel_2[0].size(); j++)
		{
			auto B = this->Krug_Yzel_2[0][j];


			double d = sqrt(kv(A->coord[0][0] - B->coord[0][0]) + kv(A->coord[0][1] - B->coord[0][1]) +
				kv(A->coord[0][2] - B->coord[0][2]));
			if (d < this->geo->R0 / 1000.0)
			{
				//cout << "Nashol " << j << endl;

				for (int k = 0; k < this->Krug_Yzel_2.size(); k++)
				{
					delete this->Krug_Yzel_2[k][j];
					this->Krug_Yzel_2[k][j] = this->C_Luch[i].back()->Yzels[k];
				}

				break;
			}
		}
	}

	// ��������� ���� � ����� ������ �����
	for (int i = 0; i < this->Krug_Yzel.size(); i++)
	{
		for (int j = 0; j < this->Krug_Yzel[0].size(); j++)
		{
			if (this->Krug_Yzel[i][j]->number < 0)
			{
				this->All_Yzel.push_back(this->Krug_Yzel[i][j]);
			}
		}
	}

	cout << "Yzlov posle dobavleniya  " << this->All_Yzel.size() << endl;

	// ���������������� ��� ����
	if (true)
	{
		int kkk = 1;
		for (auto& i : this->All_Yzel)
		{
			i->number = kkk;
			kkk++;
		}
	}

	// ���� ������� ���� �� ����� �� ������ (A2 � C2 �����)  // ---------------------------------------------------------------------------------------------
	if (true)
	{
		this->A2_Luch.reserve(this->Krug_Yzel[0].size() - 60);
		this->C2_Luch.reserve(this->Krug_Yzel_2[0].size() - 60); // ����� ������� ����� �� �������, ��� ���
			// ��� ��� ������ � ���� A

		for (int j = 0; j < this->Krug_Yzel[0].size(); j++) 
		{
			int nf = this->Krug_Yzel[0][j]->number;

			bool hasN = std::any_of(  // �������� �� ������ �������, ��������������� �������
				this->A_Luch.begin(),
				this->A_Luch.end(),
				[nf](auto x) { return x[0]->Yzels[0]->number == nf; }
			);

			if (hasN) continue;

			auto Luch_0 = new Luch();  // ������ ���
			All_Luch.push_back(Luch_0);
			Luch_0->type = "A2_Luch";
			this->A2_Luch.push_back(Luch_0);
			for (int i = 0; i < this->Krug_Yzel.size(); i++)
			{
				Luch_0->Yzels.push_back(this->Krug_Yzel[i][j]);
			}


			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 + this->geo->M1 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M2 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M2 + 1 + this->geo->M3 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M2 + 1 + this->geo->M3 + 1 + this->geo->M4 + 1 - 1][j]);
		}


		for (int j = 0; j < this->Krug_Yzel_2[0].size(); j++)
		{
			int nf = this->Krug_Yzel_2[0][j]->number;

			bool hasN = std::any_of(  // �������� �� ������ �������, ��������������� �������
				this->C_Luch.begin(),
				this->C_Luch.end(),
				[nf](auto x) { return x.back()->Yzels[0]->number == nf; }
			);

			if (hasN) continue;

			auto Luch_0 = new Luch();  // ������ ���
			All_Luch.push_back(Luch_0);
			Luch_0->type = "C2_Luch";
			this->C2_Luch.push_back(Luch_0);
			for (int i = 0; i < this->Krug_Yzel_2.size(); i++)
			{
				Luch_0->Yzels.push_back(this->Krug_Yzel_2[i][j]);
			}


			Luch_0->Yzels_opor.push_back(this->Krug_Yzel_2[this->geo->M0 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel_2[this->geo->M0 + 1 + this->geo->M1 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel_2[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M11 + this->geo->N4 + 1 - 1][j]);
			Luch_0->Yzels_opor.push_back(this->Krug_Yzel_2[this->geo->M0 + 1 + this->geo->M1 + 1 +
				this->geo->M11 + this->geo->N4 + this->geo->N5 - 1][j]);
		}
		
	}


	// ������ ���� ������� ���� � ������������ � ��������� ��������
	if (true) 
	{
		this->Set_luch_parametr();
		// ���������� ������� ����� � �2-����� � ���������� �������
		double x, y, z, r;
		for (auto& i : this->A2_Luch)
		{
			auto A = i->Yzels_opor[0];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2]; 
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R1;
			A->coord[0][1] = y / r * this->geo->R1;
			A->coord[0][2] = z / r * this->geo->R1;

			A = i->Yzels_opor[1];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R2;
			A->coord[0][1] = y / r * this->geo->R2;
			A->coord[0][2] = z / r * this->geo->R2;

			A = i->Yzels_opor[2];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R3;
			A->coord[0][1] = y / r * this->geo->R3;
			A->coord[0][2] = z / r * this->geo->R3;

			A = i->Yzels_opor[3];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R4;
			A->coord[0][1] = y / r * this->geo->R4;
			A->coord[0][2] = z / r * this->geo->R4;

			A = i->Yzels_opor[4];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R5;
			A->coord[0][1] = y / r * this->geo->R5;
			A->coord[0][2] = z / r * this->geo->R5;
		}

		for (auto& i : this->C2_Luch)
		{
			auto A = i->Yzels_opor[0];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R1;
			A->coord[0][1] = y / r * this->geo->R1;
			A->coord[0][2] = z / r * this->geo->R1;

			A = i->Yzels_opor[1];
			x = A->coord[0][0];
			y = A->coord[0][1];
			z = A->coord[0][2];
			r = sqrt(kvv(x, y, z));
			A->coord[0][0] = x / r * this->geo->R2;
			A->coord[0][1] = y / r * this->geo->R2;
			A->coord[0][2] = z / r * this->geo->R2;
		}

		// ���� �� ������������� ������ �����, ����� ������� �� � ���������� �������
		std::sort(this->All_Luch.begin(), this->All_Luch.end(), [](Luch* a, Luch* b) {
			if (a->type == "A_Luch" && b->type != "A_Luch")
			{
				return true;
			}

			if (a->type == "B_Luch" && b->type != "A_Luch" && b->type != "B_Luch")
			{
				return true;
			}

			if (a->type == "C_Luch" && b->type != "A_Luch" && 
				b->type != "B_Luch" && b->type != "C_Luch")
			{
				return true;
			}
			if (a->type == "D_Luch" && b->type != "A_Luch" &&
				b->type != "B_Luch" && b->type != "C_Luch" &&
				b->type != "D_Luch")
			{
				return true;
			}

			if (a->type == "E_Luch" && b->type != "A_Luch" &&
				b->type != "B_Luch" && b->type != "C_Luch" && 
				b->type != "D_Luch" && b->type != "E_Luch")
			{
			return true;
			}

			if (a->type == "H_Luch" && b->type != "A_Luch" &&
				b->type != "B_Luch" && b->type != "C_Luch" &&
				b->type != "D_Luch" && b->type != "E_Luch" &&
				b->type != "H_Luch")
			{
			return true;
			}

			return false;
			});

		for (auto& i : this->All_Luch)
		{
			i->dvigenie(0);
		}
		//for (auto& nn : this->name_luch) // ���-�� ���������� ���� � ���������� �������
		//{
		//	auto A = *(All_name_luch[nn]);
		//	for (auto& j : A)
		//	{
		//		for (auto& i : j)
		//		{
		//			i->dvigenie(0);
		//		}
		//	}
		//}
	}

	// ��������� ������ (��� ����� � ������)
	if (true)
	{
		this->Cell_layer_head.resize(this->Krug_Yzel.size() - 1);
		this->Cell_layer_tail.resize(this->Krug_Yzel_2.size() - 1);

		int nn, ny;
		fin.read((char*)&n, sizeof n); // ����� �����

		for (int i = 0; i < this->Cell_layer_head.size(); i++)
		{
			this->Cell_layer_head[i].resize(n);
		}

		for (int i = 0; i < this->Cell_layer_tail.size(); i++)
		{
			this->Cell_layer_tail[i].resize(n);
		}

		for (int i = 0; i < this->Cell_layer_head.size(); i++)
		{
			for (int j = 0; j < this->Cell_layer_head[0].size(); j++)
			{
				auto cell = new Cell();
				this->Cell_layer_head[i][j] = cell;
				this->All_Cell.push_back(cell);
			}
		}

		for (int i = 0; i < this->Cell_layer_tail.size(); i++)
		{
			for (int j = 0; j < this->Cell_layer_tail[0].size(); j++)
			{
				auto cell = new Cell();
				this->Cell_layer_tail[i][j] = cell;
				this->All_Cell.push_back(cell);
			}
		}


		for (int i = 0; i < n; i++)
		{
			fin.read((char*)&nn, sizeof nn); // ����� ����� � ������
			for (int j = 0; j < nn; j++)
			{
				fin.read((char*)&ny, sizeof ny); // ����� ����

				for (int k = 0; k < this->Cell_layer_head.size(); k++)
				{
					this->Cell_layer_head[k][i]->yzels.push_back(Krug_Yzel[k][ny - 1]);
					this->Cell_layer_head[k][i]->yzels.push_back(Krug_Yzel[k + 1][ny - 1]);
				}

				for (int k = 0; k < this->Cell_layer_tail.size(); k++)
				{
					this->Cell_layer_tail[k][i]->yzels.push_back(Krug_Yzel_2[k][ny - 1]);
					this->Cell_layer_tail[k][i]->yzels.push_back(Krug_Yzel_2[k + 1][ny - 1]);
				}
			}
			
		}

		// ������ ������� ����, ����� ��� ��� � ���������� �������
		if (true)
		{
			for (int i = 0; i < this->Cell_layer_head.size(); i++)
			{
				for (int j = 0; j < this->Cell_layer_head[0].size(); j++)
				{
					auto a1 = this->Cell_layer_head[i][j]->yzels[0];
					auto a2 = this->Cell_layer_head[i][j]->yzels[1];
					auto a3 = this->Cell_layer_head[i][j]->yzels[2];
					auto a4 = this->Cell_layer_head[i][j]->yzels[3];
					auto a5 = this->Cell_layer_head[i][j]->yzels[4];
					auto a6 = this->Cell_layer_head[i][j]->yzels[5];
					auto a7 = this->Cell_layer_head[i][j]->yzels[6];
					auto a8 = this->Cell_layer_head[i][j]->yzels[7];

					this->Cell_layer_head[i][j]->yzels[1] = a3;
					this->Cell_layer_head[i][j]->yzels[2] = a5;
					this->Cell_layer_head[i][j]->yzels[3] = a7;
					this->Cell_layer_head[i][j]->yzels[4] = a2;
					this->Cell_layer_head[i][j]->yzels[5] = a4;
					this->Cell_layer_head[i][j]->yzels[6] = a6;
					this->Cell_layer_head[i][j]->yzels[7] = a8;
				}
			}
		}

		
	}


	// ���� ��������� ����������� �������� �����, ������� ��������� ���������� R1 � �.�. ������� �������������� ����������
	// ����� �������� ��������� ������������ ������������ � ��������� �� ����� (������ ����� �� ������ ���������)
	
	// ������� �����, ����, ������� ��
	// ������� ����������� ������ �� �����

	fin.close();
	ffin.close();
	cout << "---END New_initial---" << endl;
}

void Setka::Winslow_method(void)
{
	int N, M;
	double R = 1.0;
	double** x;
	double** y;
	double** z;
	double** x2;
	double** y2;
	double xx, yy;
	double du, dv;
	double dxu, dxv, dyu, dyv, ddxuv, ddyuv, A, B, C, err;

	N = 16;
	M = 16;
	

	// ������ �������
	if (true)
	{
		x = new double* [N];
		for (int i = 0; i < N; i++) x[i] = new double[M];

		y = new double* [N];
		for (int i = 0; i < N; i++) y[i] = new double[M];

		x2 = new double* [N];
		for (int i = 0; i < N; i++) x2[i] = new double[M];

		y2 = new double* [N];
		for (int i = 0; i < N; i++) y2[i] = new double[M];

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				x[i][j] = 0.0;
				x2[i][j] = 0.0;
				y[i][j] = 0.0;
				y2[i][j] = 0.0;
			}
		}
	}

	// ������������� ��������
	if (true)
	{
		for (int i = 0; i < N; i++)
		{
			xx = R * cos(-const_pi/ 2.0 - const_pi/ 4.0 + (i)*const_pi/ 2.0 / (N - 1.0));
			yy = R * sin(-const_pi/ 2.0 - const_pi/ 4.0 + (i)*const_pi/ 2.0 / (N - 1.0));
			x[i][0] = xx;
			x2[i][0] = xx;
			y[i][0] = yy;
			y2[i][0] = yy;
		}

		for (int i = 0; i < N; i++)
		{
			xx = R * cos(const_pi/ 2.0 + const_pi/ 4.0 - (i)*const_pi/ 2.0 / (N - 1.0));
			yy = R * sin(const_pi/ 2.0 + const_pi/ 4.0 - (i)*const_pi/ 2.0 / (N - 1.0));
			x[i][M - 1] = xx;
			x2[i][M - 1] = xx;
			y[i][M - 1] = yy;
			y2[i][M - 1] = yy;
		}

		for (int j = 0; j < M; j++)
		{
			xx = R * cos(-const_pi/ 4.0 + (j)*const_pi/ 2.0 / (M - 1.0));
			yy = R * sin(-const_pi/ 4.0 + (j)*const_pi/ 2.0 / (M - 1.0));
			x[N - 1][j] = xx;
			x2[N - 1][j] = xx;
			y[N - 1][j] = yy;
			y2[N - 1][j] = yy;
		}

		for (int j = 0; j < M; j++)
		{
			xx = R * cos(const_pi+ const_pi/ 4.0 - (j)*const_pi/ 2.0 / (M - 1.0));
			yy = R * sin(const_pi+ const_pi/ 4.0 - (j)*const_pi/ 2.0 / (M - 1.0));
			x[0][j] = xx;
			x2[0][j] = xx;
			y[0][j] = yy;
			y2[0][j] = yy;
		}
	}

	du = 1.0 / N;
	dv = 1.0 / M;

	do
	{
		err = 0.0;
		for (int i = 1; i < N - 1; i++)
		{
			for (int j = 1; j < M - 1; j++)
			{
				dxu = (x[i + 1][j] - x[i - 1][j]) / (2.0 * du);
				dxv = (x[i][j + 1] - x[i][j - 1]) / (2.0 * dv);
				dyu = (y[i + 1][j] - y[i - 1][j]) / (2.0 * du);
				dyv = (y[i][j + 1] - y[i][j - 1]) / (2.0 * dv);
				ddxuv = (x[i + 1][j + 1] + x[i - 1][j - 1] - x[i + 1][j - 1] - x[i - 1][j + 1]) / (4.0 * du * dv);
				ddyuv = (y[i + 1][j + 1] + y[i - 1][j - 1] - y[i + 1][j - 1] - y[i - 1][j + 1]) / (4.0 * du * dv);
				A = dxv * dxv + dyv * dyv;
				B = 2.0 * (dxu * dxv + dyu * dyv) * ddxuv;
				C = dxu * dxu + dyu * dyu;
				/*if (i == 1 && j == 1)
				{
					cout << const_pi<< " " << dxu << " " << dxv << " " << A << " a " << B << " " << C << endl;
					cout << x[i + 1][j] << " " << x[i - 1][j] << " " << du << endl;
				}*/
				if (fabs(2.0 * A / (du * du) + 2.0 * C / (dv * dv)) > 0.000001)
				{
					x2[i][j] = (A * (x[i + 1][j] + x[i - 1][j]) / (du * du) - B + C * (x[i][j + 1] + x[i][j - 1]) / (dv * dv)) / (2.0 * A / (du * du) + 2.0 * C / (dv * dv));
					err = max(err, fabs(x2[i][j] - x[i][j]));
				}
				B = 2.0 * (dxu * dxv + dyu * dyv) * ddyuv;
				if (fabs(2.0 * A / (du * du) + 2.0 * C / (dv * dv)) > 0.000001)
				{
					y2[i][j] = (A * (y[i + 1][j] + y[i - 1][j]) / (du * du) - B + C * (y[i][j + 1] + y[i][j - 1]) / (dv * dv)) / (2.0 * A / (du * du) + 2.0 * C / (dv * dv));
					err = max(err, fabs(y2[i][j] - y[i][j]));
				}
			}
		}
		z = x2;
		x2 = x;
		x = z;

		z = y2;
		y2 = y;
		y = z;

		//cout << err << " " << x[1][1] << " " << x2[1][1] << endl;
		//return;
	} while (err > 0.000001);

	for (int i = 0; i < int(N/4); i++)
	{
		for (int j = 0; j < int(M/4); j++)
		{
			x[i][j] = (x[i][j] - x[N - 1 - i][M - 1 - j] + x[i][M - 1 - j] - x[N - 1 - i][j]) / 4.0;
			x[i][M - 1 - j] = x[i][j];
			x[N - 1 - i][j] = -x[i][j];
			x[N - 1 - i][M - 1 - j] = -x[i][j];

			y[i][j] = (y[i][j] - y[N - 1 - i][M - 1 - j] - y[i][M - 1 - j] + y[N - 1 - i][j]) / 4.0;
			y[i][M - 1 - j] = -y[i][j];
			y[N - 1 - i][j] = y[i][j];
			y[N - 1 - i][M - 1 - j] = -y[i][j];
		}
	}

	ofstream fout;
	fout.open("winslow.txt");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			//if( i%10 == 0 && j%10 == 0)  
			fout << x[i][j] << " " << y[i][j] << endl;
			/*else if( i == N - 1 && j % 10 == 0)  fout << x[i][j] << " " << y[i][j] << endl;
			else if(i % 10 == 0 && j == M - 1)  fout << x[i][j] << " " << y[i][j] << endl;
			else if(i == N - 1 && j == M - 1)  fout << x[i][j] << " " << y[i][j] << endl;*/
		}
	}
	fout.close();

}

void Setka::Tecplot_print_all_yzel_in_3D(string name)
{
	// name - ��� ��� �����
	int k = 1;
	ofstream fout;
	string name_f = "Tecplot_all_yzel_in_3D_" + name + ".txt";
	

	fout.open(name_f);



	fout << "TITLE = HP  VARIABLES = X, Y, Z"  << endl;

	for (auto& ii : this->Yzel_2D)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();


	name_f = "Tecplot_all_kryg_yzel_in_3D_" + name + ".txt";
	k = 1;
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->Krug_Yzel)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();


	name_f = "Tecplot_all_kryg_yzel_2_in_3D_" + name + ".txt";
	k = 1;
	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->Krug_Yzel_2)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_krug_yzel_in_3D(int num)
{
	// num - ����� ���� ������� - �������� ��� ���������
	int k = 1;
	vector <vector<Yzel*>> VVVV;
	ofstream fout;
	string name_f = "Tecplot_krug_yzel_in_3D_.txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	VVVV = Krug_Yzel;
	if (num == 2) VVVV = Krug_Yzel_2;

	for (auto& ii : VVVV)
	{
		fout << "ZONE T = HP, SOLUTIONTIME = " << k << endl;
		for (auto& i : ii)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
		k++;
		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_opor_yzel_in_luchs_3D(string name)
{
	ofstream fout;
	string name_f = "Tecplot_opor_yzel_in_luchs_3D_" + name + ".txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	for (auto& ii : this->All_Luch)
	{
		if (ii->type != name) continue;
		for (auto& i : ii->Yzels_opor)
		//for (auto& i : ii->Yzels)
		{
			fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
		}
	}

	fout.close();
}
 
void Setka::Tecplot_print_all_lush_in_3D(string name)
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_all_lush_in_3D_" + name + ".txt";

	auto A = *All_name_luch[name];
	int k = 0;
	for (auto& j : A[0])
	{
		k += (j->Yzels.size() - 1);
	}

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int kkk = 1;

	
	for(auto& i : A)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << 2 * k << ", E = " << k << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& j : i)
		{
			for (int m = 0; m < j->Yzels.size() - 1; m++)
			{
				fout << j->Yzels[m]->coord[0][0] << " " << j->Yzels[m]->coord[0][1] << " " << j->Yzels[m]->coord[0][2] << endl;
				fout << j->Yzels[m + 1]->coord[0][0] << " " << j->Yzels[m + 1]->coord[0][1] << " " << j->Yzels[m + 1]->coord[0][2] << endl;
			}
		}


		for (int m = 0; m < k; m++)
		{
			fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_lush_in_2D()
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_all_lush_in_2D.txt";

	int k = 0;

	for (auto& nn : this->name_luch)
	{
		auto A = *(this->All_name_luch[nn]);
		for (auto& j : A[0])
		{
			k += (j->Yzels.size() - 1);
		}
	}

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y" << endl;

	int kkk = 1;


	for (int i = 0; i < this->geo->Nphi; i++)
	{
		//auto i = A[1];
		fout << "ZONE T=HP, N = " << 2 * k << ", E = " << k << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& nn : this->name_luch)
		{
			auto A = *(this->All_name_luch[nn]);
			for (auto& j : A[i])
			{
				for (int m = 0; m < j->Yzels.size() - 1; m++)
				{
					fout << j->Yzels[m]->coord[0][0] << " " << sqrt(kv(j->Yzels[m]->coord[0][1]) + kv(j->Yzels[m]->coord[0][2])) << endl;
					fout << j->Yzels[m + 1]->coord[0][0] << " " << sqrt(kv(j->Yzels[m + 1]->coord[0][1]) + kv(j->Yzels[m + 1]->coord[0][2])) << endl;
				}
			}
		}


		for (int m = 0; m < k; m++)
		{
			fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}

void Setka::Tecplot_print_all_cell_in_3D()
{
	// name - ��� ��� �����
	ofstream fout;
	string name_f = "Tecplot_all_cell_in_3D_.txt";


	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	int kkk = 1;


	for (auto& i : this->Cell_layer_head[0])
	{
		/*if (kkk > 200)
		{
			break;
		}*/

		fout << "ZONE T=HP, N = " << i->yzels.size() << ", E = " << i->yzels.size() - 1 << ", F=FEPOINT, ET=LINESEG, SOLUTIONTIME = " << kkk << endl;
		kkk++;

		for (auto& j : i->yzels)
		{
			fout << j->coord[0][0] << " " << j->coord[0][1] << " " << j->coord[0][2] << endl;
		}

		for (int m = 0; m < i->yzels.size() - 1; m++)
		{
			fout << m + 1 << " " << m + 2 << endl;
		}

		fout << endl;
	}

	fout.close();
}
