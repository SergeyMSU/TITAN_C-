#include "Help_funk.h"
#include "Yzel.h"

double polar_angle(const double& x, const double& y)
{
	if (fabs(x) + fabs(y) < 0.000001)
	{
		return 0.0;
	}

	if (x < 0)
	{
		return atan(y / x) + 1.0 * const_pi;
	}
	else if (x > 0 && y >= 0)
	{
		return atan(y / x);
	}
	else if (x > 0 && y < 0)
	{
		return atan(y / x) + 2.0 * const_pi;
	}
	else if (y > 0 && x >= 0 && x <= 0)
	{
		return const_pi / 2.0;
	}
	else if (y < 0 && x >= 0 && x <= 0)
	{
		return  3.0 * const_pi / 2.0;
	}
	return 0.0;
}


double Yzel_distance(Yzel* A, Yzel* B, int time)
{
	return sqrt(kv(A->coord[time][0] - B->coord[time][0]) +
		kv(A->coord[time][1] - B->coord[time][1]) + kv(A->coord[time][2] - B->coord[time][2]));
}

double Yzel_distance_x(Yzel* A, Yzel* B, int time)
{
	return sqrt(kv(A->coord[time][1] - B->coord[time][1]) + kv(A->coord[time][2] - B->coord[time][2]));
}

bool areCellsEqual(const Gran& cell1, const Gran& cell2) {
	if (cell1.yzels.size() != cell2.yzels.size()) {
		return false;
	}

	// Создаем мультимножество для каждого вектора
	std::unordered_multiset<Yzel*> set1(cell1.yzels.begin(), cell1.yzels.end());
	std::unordered_multiset<Yzel*> set2(cell2.yzels.begin(), cell2.yzels.end());

	return set1 == set2;
}

bool areCellsEqual(const Gran* cell1, const Gran* cell2) 
{
	if (cell1->yzels.size() != cell2->yzels.size()) 
	{
		return false;
	}

	// Создаем мультимножество для каждого вектора
	std::unordered_multiset<Yzel*> set1(cell1->yzels.begin(), cell1->yzels.end());
	std::unordered_multiset<Yzel*> set2(cell2->yzels.begin(), cell2->yzels.end());

	return set1 == set2;
}

bool areCellsEqual_my(const Gran* cell1, const Gran* cell2)
{
	if (cell1->yzels.size() != cell2->yzels.size())
	{
		return false;
	}

	for (int i = 0; i < cell1->yzels.size(); i++)
	{
		if (std::find(cell2->yzels.begin(), cell2->yzels.end(), cell1->yzels[i]) == cell2->yzels.end())
		{
			return false;
		}
	}

	return true;
}


double triangleArea3D(
	const double & x1, const double& y1, const double& z1,
	const double& x2, const double& y2, const double& z2,
	const double& x3, const double& y3, const double& z3)
{
	// Вычисляем векторы сторон
	double v1x = x2 - x1;
	double v1y = y2 - y1;
	double v1z = z2 - z1;

	double v2x = x3 - x1;
	double v2y = y3 - y1;
	double v2z = z3 - z1;

	// Вычисляем векторное произведение
	double cross_x = v1y * v2z - v1z * v2y;
	double cross_y = v1z * v2x - v1x * v2z;
	double cross_z = v1x * v2y - v1y * v2x;

	// Вычисляем длину вектора произведения
	double cross_magnitude = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);

	// Площадь равна половине длины векторного произведения
	return 0.5 * cross_magnitude;
}


void crossProductFast(
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz,
	double& rx, double& ry, double& rz)
{
	rx = ay * bz - az * by;
	ry = az * bx - ax * bz;
	rz = ax * by - ay * bx;
}

double scalarProductFast(
	const double& ax, const double& ay, const double& az,
	const double& bx, const double& by, const double& bz)
{
	return ax * bx + ay * by + az * bz;
}

double tetrahedronVolume(
	const double& x1, const double& y1, const double& z1,  // Вершина A
	const double& x2, const double& y2, const double& z2,  // Вершина B
	const double& x3, const double& y3, const double& z3,  // Вершина C
	const double& x4, const double& y4, const double& z4)  // Вершина D
{
	// Вычисляем векторы рёбер из вершины A
	double abx = x2 - x1, aby = y2 - y1, abz = z2 - z1;  // Вектор AB
	double acx = x3 - x1, acy = y3 - y1, acz = z3 - z1;  // Вектор AC
	double adx = x4 - x1, ady = y4 - y1, adz = z4 - z1;  // Вектор AD

	// Вычисляем векторное произведение AB × AC
	double cross_x = aby * acz - abz * acy;
	double cross_y = abz * acx - abx * acz;
	double cross_z = abx * acy - aby * acx;

	// Вычисляем скалярное произведение (AB × AC) · AD
	double dot_product = cross_x * adx + cross_y * ady + cross_z * adz;

	// Объём равен 1/6 модуля этого произведения
	return std::abs(dot_product) / 6.0;
}