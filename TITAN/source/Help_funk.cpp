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