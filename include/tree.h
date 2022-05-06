struct Tree_Node {
	double x;
	double v;
	double q[4];
	struct Tree_Node* pNext[4];
};

struct Tree {
	struct Tree_Node* pHead;
	struct Tree_Node* pTops[101];
	struct Tree_Node* pBottoms[101];
	int n_nodes[101];
};
