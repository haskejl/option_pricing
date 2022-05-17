#include <math.h>
#include <stdio.h>
#include <time.h>

#include "./include/math_funcs.h"
#include "./include/tree.h"

// IBM parameters
const double alpha = 11.85566;
const double nu = 0.9345938;
const double psi = 4.13415;
const double mu = 0.04588;

const int M = 300;
const double h = 1.f/252.f;
const int n = 1000;
const double T = 42.f/252.f;
const double r = 0.0343;
// NOTE: when changing this value, the sizes of the arrays in tree.h must also be changed
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

double calc_phi(const double x) {
	if(x > -0.1 && x < 0.1)
		return 10.f*(1-fabs(x*10.f));
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

void find_qs(const struct Tree_Node* node, const int node_j, 
						 const double mul_pt, double* q) {
	const double d1 = node->x - node_j*mul_pt;
	const double d2 = node->x - (node_j-1)*mul_pt;
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

double gen_tree(const double* Y_bar, struct Tree_Node* nodes,
								const int list_size) {
	static int size_mul = 1;
	struct Tree tree;
	struct Tree_Node* curr_node = nodes;
	int nodes_used = 1;
	
	//Setup the tree
	tree.pHead = nodes;
	tree.pTops[0] = nodes;
	tree.pTops[1] = nodes+1;
	tree.pBottoms[0] = nodes;
	tree.pBottoms[1] = nodes+4;
	tree.n_nodes[0] = 1;
	tree.n_nodes[1] = 4;
	curr_node->x = X_0;

	double sig = calc_sigma(Y_bar[0]);
	double add_pt = (r-pow(sig, 2.f)/2.f)*dt;
	double mul_pt = sig*sqrt(dt);
	double q[4];

	//Calculate the successors of the first
	int node_j = (int)ceil(X_0/mul_pt);

	find_qs(curr_node, node_j, mul_pt, q);
	
	// Setup the node and set children's x value
	struct Tree_Node* curr_child = tree.pTops[1];
	for(int i=0; i<4; i++) {
		curr_child->x = (node_j+1-i)*mul_pt+add_pt;
		curr_node->q[i] = q[i];
		curr_node->pNext[i] = curr_child;
		curr_child++;
	}
	nodes_used += 4;
	
	int j_upp;
	int j_downn;

	for(int i=1; i<N; i++) {
		// Initial calculations for time moment
		sig = calc_sigma(Y_bar[i]);
		add_pt = (r-pow(sig, 2.f)/2.f)*dt;
		mul_pt = sig*sqrt(dt);

		j_upp = ceil(tree.pTops[i]->x/mul_pt)+2; 
		j_downn = ceil(tree.pBottoms[i]->x/mul_pt)-2;

		// Quit if tree won't fit to prevent seg fault
		if(nodes_used+j_upp-j_downn > list_size*size_mul) {
			printf("ERROR: Amount of memory too small\n");
			exit(-1);
		}
		// move the current node to the top node of of the time moment
		curr_node = tree.pTops[i];
		tree.pTops[i+1] = tree.pBottoms[i]+1;
		curr_child = tree.pTops[i+1];
		tree.n_nodes[i+1] = 4;

		node_j = ceil(curr_node->x/mul_pt);
		find_qs(curr_node, node_j, mul_pt, q);

		for(int j=0; j<4; j++) {
			curr_child->x = (node_j+1-j)*mul_pt+add_pt;
			curr_node->q[j] = q[j];
			curr_node->pNext[j] = curr_child;
			curr_child++;
		}
		nodes_used += 4;
		
		
		//handle each node in the time moment
		int last_j;
		while(curr_node != tree.pBottoms[i]) {
			//next node here in order to have the operation performed on the bot node
			curr_node++;
			last_j = node_j;
			node_j = ceil(curr_node->x/mul_pt);
			find_qs(curr_node, node_j, mul_pt, q);
			//Interlace values
			int j_diff = last_j - node_j;
			if(j_diff < 0) {
				printf("\nInvalid node value, current j cannot be higher than the last!, %d, %d\n", last_j, node_j);
				printf("%f, %f\n", (curr_node-1)->x, curr_node->x);
				exit(-1);
			}

			// match up the same child nodes
			j_diff = j_diff > 4 ? 4 : j_diff;
			curr_child -= (4-j_diff);

			for(int j=0; j<4; j++) {
				//only need to calculate this if it isn't an overlapping node
				if(3-j_diff < j)
					curr_child->x = (node_j+1-j)*mul_pt+add_pt;
				curr_node ->q[j] = q[j];
				curr_node ->pNext[j] = curr_child;
				curr_child++;
			}
			tree.n_nodes[i+1] += j_diff;
			nodes_used += j_diff;
		}
		tree.pBottoms[i+1] = tree.pTops[i+1]+(tree.n_nodes[i+1]-1);
	}
	curr_node = tree.pTops[N];
	while(curr_node != tree.pBottoms[N]){
		curr_node->v = payoff_func(exp(curr_node->x), E, r, T);
		curr_node++;
	}
	//bottom node
	curr_node->v = payoff_func(exp(curr_node->x), E, r, T);
	for(int i=N-1; i>0; i--) {
		curr_node = tree.pTops[i];
		while(curr_node != tree.pBottoms[i]){
			curr_node->v = 0.f; // initialize data
			for(int j=0; j<4; j++){
				curr_node->v += curr_node->pNext[j]->v*curr_node->q[j];
			} 
			curr_node++;
		}
		//bottom node
		curr_node->v = 0.f; // initialize data
		for(int j=0; j<4; j++)
			curr_node->v += curr_node->pNext[j]->v*curr_node->q[j];
	}
	curr_node = tree.pTops[0];
	curr_node->v = 0.f; // initialize data
	for(int j=0; j<4; j++)
		curr_node->v += curr_node->pNext[j]->v*curr_node->q[j];
	return curr_node->v;
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
	// allocate 1M 80byte objects initially (80MB)
	clock_t start_mc, end_mc;
	start_mc = clock();
	// need enough to run this with N=100
	const int num_nodes = 10000000;

	struct Tree_Node* nodes = malloc(sizeof(struct Tree_Node)*num_nodes);

	if(nodes == 0) {
		printf("ERROR: memory allocation failed");
		exit(-1);
	}
	double expected_val[7];
	double Y_bar[N];
	int nruns = 100;
	int n_st = 0;
	while(cdf[n_st] == 0) {
		n_st++;
		if(n_st == n) {
			printf("Problem exists with generating probability distribution");
			exit(-1);
		}
	}

	for(int i=0; i<nruns; i++) {
		// Generate the values of Y_bar before making the tree
		for(int j=0; j<N; j++) {
			double rn = gen_unif_dist_rn(0.f, 1.f);
			for(int k=n_st; k<n; k++) {
				if(rn <= cdf[k]) {
					Y_bar[j] = y_curr[k];
					break;
				}
			}
		}
		printf("\rRun %d/%d", i+1, nruns);
		for(int j=0; j<7; j++) {
			E = evals[j];
			expected_val[j] += gen_tree(Y_bar, nodes, num_nodes);
		}
	}
	end_mc = clock();
	for(int j=0; j<7; j++) {
		printf("E=%f\t%f\n", evals[j], expected_val[j]/(double)nruns);
	}
	printf("MC runtime: %f\n", (double)(end_mc-start_mc)/CLOCKS_PER_SEC);

	free(nodes);
}
