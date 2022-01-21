#ifndef DRAND48_H
#define DRAND48_H


#include<stdlib.h>
static unsigned long long seed_pre = 0;
#define const_a 0x5DEECE66DLL
#define const_c 0xB
#define DRAND48_MAX 0xFFFFFFFFFFFFLL
#define LRAND_MAX 0xFFFFFFFFL

void srand48(long int seedval)
{
	seed_pre = ((long long int)seedval)<<16 | 0x00000000330ELL;
}

double drand48()
{
	seed_pre = (const_a * seed_pre + const_c) % (DRAND48_MAX+1);
	unsigned long int seed = seed_pre >> 16;
	return seed / (LRAND_MAX + 1.);
}

unsigned long int lrand48()
{
	seed_pre = (const_a * seed_pre + const_c) % (DRAND48_MAX+1);
	return seed_pre >> 16;
}
/*#define PRECISION 2.82e14
double drand48()
{
	double x = 0;
	double denom = RAND_MAX + 1.;
	double need;
	for (need = PRECISION; need > 1; need /= (RAND_MAX + 1.))
	{
		x += rand() / denom;
		denom *= (RAND_MAX + 1.);
	}
	return x;
}*/

#endif