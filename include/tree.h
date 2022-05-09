struct Tree_Node {
	double x;
	double v;
	double q[4];
	struct Tree_Node* pNext[4];
};

struct Tree {
	struct Tree_Node* pHead;
	struct Tree_Node* pTops[2001];
	struct Tree_Node* pBottoms[2001];
	int n_nodes[2001];
};
