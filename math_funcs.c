#include "./include/math_funcs.h"

double gen_norm_dist_rn(const double mean, const double sd) {
    static bool haveSpare = false;
    static double spare;

    if(haveSpare)
    {
        haveSpare = false;
        return spare;
    }
    haveSpare = true;

    double u1, u2;
    do
    {
        u1 = rand()/((double)RAND_MAX);
        u2 = rand()/((double)RAND_MAX);
    }while(u1 < DBL_EPSILON);

    double mag = sd * sqrt(-2.0 * log(u1));
    spare = mag * sin(2*M_PI * u2) + mean;

    return mag * cos(2*M_PI * u2) + mean;
}

double gen_unif_dist_rn(const double a, const double b) {
    double u1;
    do
    {
        u1 = rand()/((double)RAND_MAX);
    }while(u1 < DBL_EPSILON);
		return (b-a)*u1+a;
}

double max(const double x, const double y) {
	return (x > y) ? x : y;
}

double min(const double x, const double y) {
	return (x < y) ? x : y;
}
