#ifndef MATH_FUNCS_H
#define MATH_FUNCS_H

#define _USE_MATH_DEFINES
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

double gen_norm_dist_rn(const double mean, const double sd);
double gen_unif_dist_rn(const double a, const double b);

double max(const double a, const double b);
double min(const double a, const double b);

#endif
