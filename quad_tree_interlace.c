#include <math.h>
#include <stdio.h>
#include <time.h>

#include "./include/math_funcs.h"
//#include "./include/ibmxvals.h"

// IBM parameters
const double alpha = 11.85566;
const double nu = log(0.234);//0.9345938;
const double psi = 1.f;//4.13415;
const double mu = 0.04588;

const int M = 300;
const double h = 1.f/252.f;
const int n = 1000;
const double T = 42.f/252.f;
const double r = 0.0343;
const int N = 100;
const double evals[] = {60.f, 70.f, 75.f, 80.f, 85.f, 90.f, 95.f};
double E = 95.f;
double p = 0.135;
const double dt = T/(double)N;
const double X_0 = log(83.7);

double xvals[] = {
	76.65,
    75.48,
    72.01,
    74.03,
    74.21,
    74.61,
    75.43,
    77.05,
    75.91,
    76.38,
    76.51,
    76.47,
    77.08,
    75.50,
    75.26,
    74.98,
    73.30,
    73.28,
    72.62,
    73.16,
    74.34,
    74.29,
    76.36,
    77.16,
    76.41,
    76.51,
    75.81,
    76.00,
    77.14,
    77.10,
    75.55,
    76.84,
    77.35,
    75.79,
    75.00,
    75.04,
    74.80,
    74.93,
    74.77,
    75.05,
    74.89,
    76.30,
    77.05,
    76.39,
    76.55,
    76.41,
    77.23,
    75.41,
    74.01,
    73.88,
    75.30,
    74.73,
    74.20,
    74.67,
    74.79,
    75.81,
    77.38,
    79.30,
    78.96,
    80.04,
    81.45,
    82.42,
    82.38,
    81.81};
const int x_len = 64;
//const int x_len = 1329;
//const int x_len = 42;
/*double xvals[] = {76.36,77.16,76.41,76.51,75.81,76.0,77.14,77.1,75.55,76.84,77.35,75.79,75.0,75.04,74.8,
				74.93,74.77,75.05,74.89,76.3,77.05,76.39,76.55,76.41,77.23,75.41,74.01,73.88,75.3,74.73,
				74.2,74.67,74.79,75.81,77.38,79.3,78.96,80.04,81.45,82.42,82.38,81.81};*/

double calc_phi(const double x) {
	if(x > -1.f && x < 1.f)
		return 1-fabs(x);
	else
		return 0.f;
}

void calc_cdf(double* cdf, const double* X, const double xn, const int n) {
	double C = 0.f;
	for(int i=0; i<n; i++) {
		C += calc_phi(X[i]-xn);
	}
	cdf[0] = calc_phi(X[0]-xn)/C;
	for(int i=1; i<n; i++) {
		cdf[i] = calc_phi(X[i]-xn)/C + cdf[i-1];
	}
}

double calc_sigma(const double x) {
	return exp(-fabs(x));
}

double payoff_func(const double S, const double E, const double r, const double T) {
	return max(S-E, 0)*exp(-r*T);
}

double gen_y_val(const double* cdf, const double* y, const int n) {
	double val = gen_unif_dist_rn(0.f, 1.f);
	int i=0;
	while(i < n) {
		if(cdf[i] == 0.f) i++;
		else break;
	}
	for(int i=0; i<n; i++) {
		if(cdf[i] >= val) {
			return y[i];
		}
	}
	return y[n-1];
}

void find_qs(const double x, const int node_j, const double mul_pt, double* q) {
	const double d1 = x - node_j*mul_pt;
	const double d2 = x - (node_j-1)*mul_pt;
	if(-d1 < d2) {
		const double b = d1/mul_pt;
		const double bsq = pow(b, 2);
		q[3] = p;
		q[0] = 0.5*(1.f+b+bsq)-p;
		q[1] = 3.f*p-bsq;
		q[2] = 0.5*(1.f-b+bsq)-3*p;
	}
	else {
		const double b = d2/mul_pt;
		const double bsq = pow(b, 2);
		q[0] = p;
		q[1] = 0.5*(1.f+b+bsq)-3*p;
		q[2] = 3.f*p-bsq;
		q[3] = 0.5*(1.f-b+bsq)-p;
	}
}

double gen_tree(const double* Y_bar, double** xlist, double** xlist2, 
								double** qlist, double** qlist2, const int list_size) {
	double* curr_x = *xlist2;
	double* curr_q = *qlist2;
	double* next_x = *xlist;
	double* next_q = *qlist;
	static int size_mul = 1;

	double sig = calc_sigma(Y_bar[0]);
	double add_pt = (r-pow(sig, 2.f)/2.f)*dt;
	double mul_pt = sig*sqrt(dt);
	double q[4];

	//Calculate the successors of the first
	int node_j = (int)ceil(X_0/mul_pt);

	find_qs(X_0, node_j, mul_pt, q);

	for(int i=0; i<4; i++) {
		next_x[i] = (node_j+1.f-i)*mul_pt+add_pt;
		next_q[i] = q[i];
	}
	
	int top;
	int j_upp;
	int j_downn;
	int nodes_used = 4;

	for(int i=1; i<N; i++) {
		sig = calc_sigma(Y_bar[i]);
		add_pt = (r-pow(sig, 2.f)/2.f)*dt;
		mul_pt = sig*sqrt(dt);

		top = nodes_used; //nodes_used is over the highest point by 1
		j_upp = ceil(next_x[0]/mul_pt)+2; //use next_x for now
		j_downn = ceil(next_x[top-1]/mul_pt)-2;

		// resize the array if it won't fit
		if(j_upp-j_downn+1 > list_size*size_mul) {
			printf("\nRESIZING...\n");
			while(j_upp-j_downn+1 > list_size*size_mul)
				size_mul++;
			*xlist = realloc(*xlist, sizeof(double)*list_size*size_mul);
			*qlist = realloc(*qlist, sizeof(double)*list_size*size_mul);
			*xlist2 = realloc(*xlist2, sizeof(double)*list_size*size_mul);
			*qlist2 = realloc(*qlist2, sizeof(double)*list_size*size_mul);
			if(*xlist == 0 || *qlist == 0 || *xlist2 == 0 || *qlist2 == 0) {
				printf("ERROR: memory reallocation failed");
				exit(-1);
			}
		}
		if(i % 2 == 0) {
			curr_x = *xlist2;
			curr_q = *qlist2;
			next_x = *xlist;
			next_q = *qlist;
		}
		else {
			curr_x = *xlist;
			curr_q = *qlist;
			next_x = *xlist2;
			next_q = *qlist2;
		} 

		nodes_used = 4;
		node_j = (int)ceil(curr_x[0]/mul_pt);
		next_x[0] = (node_j+1.f)*mul_pt+add_pt;
		next_x[1] = node_j*mul_pt+add_pt;
		next_x[2] = (node_j-1.f)*mul_pt+add_pt;
		next_x[3] = (node_j-2.f)*mul_pt+add_pt;

		find_qs(curr_x[0], node_j, mul_pt, q);
		next_q[0] = curr_q[0]*q[0];
		next_q[1] = curr_q[0]*q[1];
		next_q[2] = curr_q[0]*q[2];
		next_q[3] = curr_q[0]*q[3];
		
		// find the new probabilities
		for(int j=1; j<top; j++) {
			int last_j = node_j;	
			node_j = (int)ceil(curr_x[j]/mul_pt);
			int j_diff = last_j-node_j;
			j_diff = j_diff > 4 ? 4 : j_diff;
			int pos = nodes_used-(4-j_diff);

			find_qs(curr_x[j], node_j, mul_pt, q);

			for(int k=0; k<4; k++) {
				if(3-j_diff < k) {
					next_x[pos] = (node_j+1.f-k)*mul_pt+add_pt;
					next_q[pos] = curr_q[j]*q[k];
				}
				else
					next_q[pos] += curr_q[j]*q[k];
				pos++;
			}
			nodes_used += j_diff;
		}
	}
	
	// Use the last calculated values
	curr_x = next_x;
	curr_q = next_q;
	double expected_val = 0;
	for(int j=0; j<nodes_used; j++) {
		expected_val += payoff_func(exp(curr_x[j]), E, r, T)*curr_q[j];
	}
	return expected_val;
}

int main() {
	clock_t start_prob, end_prob;
	start_prob = clock();
	time_t tim;
	srand((unsigned) time(&tim));

	for(int i=0; i<x_len; i++) {
		xvals[i] = log(xvals[i]);
	}
	
	double x_res[n];
	double y_res[n];
	double y_res2[n];
	double cdf[n];
	// Generate probability distribution
	printf("Calculating cdf...\n");

	// Step 0 - Mutation
	double X;
	double Y;
	for(int k=0; k<n; k++) {
		X = xvals[0];
		Y = nu;
		for(int j=0; j<M; j++) {
			const double sig = calc_sigma(Y);
			Y = Y + h/M * alpha * (nu-Y) + sqrt(h/M) * psi * gen_norm_dist_rn(0.f, 1.f) ;
			X = X + h/M * (mu-pow(sig,2)/2) + sqrt(h/M) * sig * gen_norm_dist_rn(0.f, 1.f);
		}
		x_res[k] = X;
		y_res[k] = Y;
	}

	// Step 0 - Selection
	for(int k=0; k<n; k++) {
		calc_cdf(&cdf[0], &x_res[0], xvals[1], n);
	}
	double* y_curr;
	double* y_past;
	// Steps 1 to K-1
	for(int i=1; i<x_len-1; i++) {
		if(i % 2 == 0) {
			y_past = &y_res2[0];
			y_curr = &y_res[0];
		} else {
			y_past = &y_res[0];
			y_curr = &y_res2[0];
		}
		// Mutation
		for(int k=0; k<n; k++) {
			X = xvals[i];
			Y = gen_y_val(&cdf[0], y_past, n);
			for(int j=0; j<M; j++) {
				const double sig = calc_sigma(Y);
				Y = Y + h/M * alpha * (nu-Y) + sqrt(h/M) * psi * gen_norm_dist_rn(0.f, 1.f);
				X = X + h/M * (mu-pow(sig,2)/2) + sqrt(h/M) * sig * gen_norm_dist_rn(0.f, 1.f);
			}
			x_res[k] = X;
			y_curr[k] = Y;
		}

		// Selection
		for(int k=0; k<n; k++) {
			calc_cdf(&cdf[0], &x_res[0], xvals[i+1], n);
		}
	}

	if(x_len % 2 == 0) {
		y_curr = &y_res[0];
	}
	else {
		y_curr = &y_res2[0];
	}
	end_prob = clock();
	double avg_sigma = cdf[0];
	for(int i=1; i<n; i++) {
		avg_sigma += (cdf[i]-cdf[i-1])*calc_sigma(y_curr[i]);
	}
	printf("Average Sigma = %f\n", avg_sigma);
	printf("Generation time = %f\n", (double)(end_prob-start_prob)/CLOCKS_PER_SEC);
	
	// Generate the tree
	// allocate 1,000 doubles initially (8kb)
	clock_t start_mc, end_mc;
	start_mc = clock();
	const int list_size = 1000;
	
	double* xlist = malloc(sizeof(double)*list_size);
	double* qlist = malloc(sizeof(double)*list_size);
	double* xlist2 = malloc(sizeof(double)*list_size);
	double* qlist2 = malloc(sizeof(double)*list_size);

	if(xlist == 0 || qlist == 0 || xlist2 == 0 || qlist2 == 0) {
		printf("ERROR: memory allocation failed");
		exit(-1);
	}
	double expected_val[7];
	double Y_bar[N];
	int nruns = 100;
	for(int i=0; i<nruns; i++) {
		// Generate the values of Y_bar before making the tree
		for(int j=0; j<N; j++) {
			double rn = gen_unif_dist_rn(0.f, 1.f);
			for(int k=0; k<n; k++) {
				while(cdf[j] == 0)	{
					k++;
					if(k == n) {
						printf("Problem exists with generating probability distribution");
						exit(-1);
					}
				}
				if(rn <= cdf[k]) {
					Y_bar[j] = y_curr[k];
					break;
				}
			}
		}
		printf("\rRun %d/%d", i+1, nruns);
		for(int numer = 0; numer<7; numer++) {
			//p = (16.f+numer)/192.f;
			E = evals[numer];
			expected_val[numer] += gen_tree(Y_bar, &xlist, &xlist2, &qlist, &qlist2, list_size);
		}
	}
	end_mc = clock();
	printf("\n");
	for(int i = 0; i<7; i++)
		printf("E=%f, Expected Value: %f\n", evals[i], expected_val[i]/(double)nruns);
	printf("MC runtime: %f\n", (double)(end_mc-start_mc)/CLOCKS_PER_SEC);

	free(xlist);
	free(qlist);
	free(xlist2);
	free(qlist2);
}
