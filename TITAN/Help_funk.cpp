#include "Help_funk.h"

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

