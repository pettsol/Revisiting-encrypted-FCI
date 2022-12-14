#include <iostream>

#include "acc_efci.h"
#include "Rabbit/rabbit.h"

int main()
{
	double a = -0.3;
	double b = 0.4;

	uint64_t gamma = 10000000000;

	uint64_t c, d;

	rho(&c, a, gamma);
	rho(&d, b, gamma);

	uint64_t s1 = 0xffff15502fffffff;
	uint64_t s2 = 0xaabf35356fffffff;
	uint64_t s3 = s1 + s2;

	double e, f;

	inv_rho(&e, c, gamma);
	inv_rho(&f, d, gamma);

	std::cout << "e = " << e << std::endl;
	std::cout << "f = " << f << std::endl;

	double g;
	uint64_t h = (c + s1) + (d + s2);
	h = h - s3;

	inv_rho(&g, h, gamma);

	std::cout << "g = " << g << std::endl;
}
