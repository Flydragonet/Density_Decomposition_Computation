#include <bits/stdc++.h>
#include <unordered_map>
#include <unordered_set>
using namespace std;
const int INF = 2000000000;

struct Timer {
	double start_time, end_time;
	void start() { start_time = clock(); }
	void end() { end_time = clock(); }
	double time() { return (end_time - start_time) / CLOCKS_PER_SEC; }
};
struct Old_Queue
{
	int* node; int head, tail;
	Old_Queue() {}
	Old_Queue(int size) { node = (int*)malloc(size * sizeof(int)); head = tail = 0; }
	~Old_Queue() { free(node); }
	bool empty() { return head == tail; }
	int pop() { return node[head++]; }
	void push(int x) { node[tail++] = x; }
};
struct Old_Set
{
	int* nodes;
	bool* in;
	int size, E_number;
	bool have_found;
	Old_Set() {}
	Old_Set(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); E_number = 0; }
	void reset(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); E_number = 0; }
	void free_memory() { free(nodes), free(in); }
	void push(int x) { nodes[size++] = x; in[x] = true; }
};
template <class T>
struct Set {
	T* nodes; bool* in; int size = -1;
	Set() {}
	Set(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(T)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void alloc(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void insert(T x) { nodes[size++] = x; in[x] = true; }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false; size = 0; }
	~Set() { free(nodes), free(in); }
};
template <class T>
struct Map {
	int* nodes; bool* in; int size = -1; T* value;
	Map() {}
	Map(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void alloc(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void freememory() { free(nodes), free(in), free(value); }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false, value[nodes[i]] = 0; size = 0; }
	T& operator[](int x) { if (!in[x]) nodes[size++] = x, in[x] = true; return value[x]; }
	~Map() { free(nodes), free(in), free(value); }
};
template <class T>
struct Queue {
	T* nodes; int head, tail;
	Queue() {}
	Queue(int size) { nodes = (T*)malloc(size * sizeof(T)); head = tail = 0; }
	void alloc(int sz) { head = tail = 0; nodes = (int*)malloc(sz * sizeof(int)); }
	~Queue() { free(nodes); }
	bool empty() { return head == tail; }
	int pop() { return nodes[head++]; }
	void push(T x) { nodes[tail++] = x; }
	void clear() { head = tail = 0; }
};
struct Stack
{
	int* node; int head;
	Stack() {}
	void init(int size) { node = (int*)malloc(size * sizeof(int)); head = 0; }
	~Stack() { free(node); }
	bool empty() { return head == 0; }
	int pop() { return node[--head]; }
	void push(int x) { node[head++] = x; }
	void clear() { head = 0; }
};
template <class T>
struct List_Node {
	List_Node<T>* pre, * next;
	T value;
	List_Node() : pre(nullptr), next(nullptr), value(T()) {}
	List_Node(const T& val) : pre(nullptr), next(nullptr), value(val) {}
};
template <class T>
struct List {
	List_Node<T>* head_node;
	int size;
	List() { head_node = new List_Node<T>, head_node->pre = head_node->next = nullptr; size = 0; }
	~List() { clear(); }
	void clear() {
		List_Node<T>* current = head_node->next;
		while (current != nullptr) {
			List_Node<T>* temp = current;
			current = current->next;
			delete temp;
		}
		head_node->next = nullptr, size = 0;
	}
	void insert(const T& value) {
		List_Node<T>* new_node = new List_Node<T>(value);
		new_node->next = head_node->next, new_node->pre = head_node;
		head_node->next = new_node;
		if (new_node->next != nullptr) new_node->next->pre = new_node;
		size++;
	}
	void insert(List_Node<T>* new_node) {
		new_node->next = head_node->next, new_node->pre = &head_node;
		head_node->next = new_node;
		if (new_node->next != nullptr) new_node->next->pre = new_node;
		size++;
	}
	void delete_node(List_Node<T>* node) {
		if (node == nullptr) return;
		if (node->pre != nullptr) node->pre->next = node->next;
		if (node->next != nullptr) node->next->pre = node->pre;
		delete node;
		size--;
	}
};

struct Edge { int t1, t2, to; };
struct Graph
{
	Graph() {}
	void read_edge(char* graph_address);
	void init_orientation();
	void decomposition();
	void output_distribution(char* graph_name, FILE* distribution);

	Edge* e;
	int n, m;
	int* d, * r, * c;
	int** adj;
	int* undeg;
	int test_value, pseudoarboricity;

	int get_max_d();
	void dfs_layer(int X_rank, int Y_rank); //Give R[X_rank] and R[Y_rank], find all R between range (Y_rank, X_rank)
	void ReTest(int Z_rank); //find R[Z_rank]
	Old_Set* R;
	int* dist;
	int* cur;
	int* p;
	int* high_indegree, high_num;
	bool DinicBFS(int X_rank, int Y_rank);
	bool DinicDFS(int);
	void update_r(int X_rank, int Y_rank, int Z_rank);

	void check_correctness(); void dfs_correctness(int x, int fa);
	void display_every_r();

	void read_update_edge_from_file(char* update_edge_address);
	void generate_update_edge(int UPDATE_EDGE_NUMBER);

	int answer_k;
	Queue<int> Q; Set<int> vis; Set<int> D; Map<int> parent;
	int* edge_label;
	void output_dense_subgraph();
	int * father;
	int find_root(int x) {
		if (father[x] == x) return x; return father[x] = find_root(father[x]);
	}
	void merge(int x, int y) {
		int root_x = find_root(x), root_y = find_root(y);
		father[root_x] = root_y;
	}

	bool* deleted;
	Set<int> idn_change_node, pend;
	void process_insertion();
	void process_deletion();
	int dfs_clock, scc_cnt, r0;
	int scc_dfs(int now);
	unordered_set<int> can_get_r0;
	unordered_map<int, int> dfsn;
	unordered_map<int, int> sccno;
	vector<vector<int> > scc;
	Stack S;

	void construct_bucket();
	unordered_set<int>* bucket;
	unordered_set<int>* index_adj;
	unordered_set<int> index_e;

	~Graph()
	{
		free(e); free(d); free(r); free(adj); free(undeg); free(cur); free(dist); free(p); free(high_indegree); free(c);
	}
};
Graph G;

enum Algorithm {
	BASIC, IMPROVE
};
Algorithm algorithm_used;
inline void check(bool flag, const char* message) {
	if (!flag) {
		printf("!!!!! CHECK ERROR !!!!!\n");
		printf("Error message: %s\n", message);
		assert(0);
	}
}

vector<int> update_edge; int now_update_edge;
void Graph::read_update_edge_from_file(char* update_edge_address) {
	FILE* file = fopen(update_edge_address, "r");
	check(file != NULL, "Can not open file update_edge_address\n");
	int t1, t2;
	while (fscanf(file, "%d %d", &t1, &t2) == 2) {
		for (int j = 0; j < undeg[t1]; j++) {
			Edge& ne = e[adj[t1][j]];
			if ((ne.t1 == t1 && ne.t2 == t2) || (ne.t2 == t1 && ne.t1 == t2)) {
				update_edge.push_back(adj[t1][j]);
				break;
			}
		}
	}
}
void Graph::generate_update_edge(int UPDATE_EDGE_NUMBER) {
	set<int> unique;
	srand(0);
	while (unique.size() < UPDATE_EDGE_NUMBER)
		unique.insert(rand() % m);
	for (auto i : unique)
		update_edge.push_back(i);
}

void Graph::read_edge(char* graph_address)
{
	const int ARRAY_SIZE_INCREMENT = 500000;
	FILE* file;
	file = fopen(graph_address, "r");

	int Emax = ARRAY_SIZE_INCREMENT;
	n = m = 0;
	e = (Edge*)malloc(Emax * sizeof(Edge));

	char line[200];
	while (fgets(line, 200, file))
	{
		int i = 0;
		e[m].t1 = e[m].t2 = 0;
		while (line[i] < '0' || line[i]>'9') i++;
		while (line[i] >= '0' && line[i] <= '9') e[m].t1 = e[m].t1 * 10 + line[i] - '0', i++;
		while (line[i] < '0' || line[i]>'9') i++;
		while (line[i] >= '0' && line[i] <= '9') e[m].t2 = e[m].t2 * 10 + line[i] - '0', i++;
		n = max(n, max(e[m].t1, e[m].t2));
		m++;
		if (m == Emax)
		{
			Emax += ARRAY_SIZE_INCREMENT;
			e = (Edge*)realloc(e, Emax * sizeof(Edge));
		}
	}

	fclose(file);
	n++;
	Q.alloc(n); vis.alloc(n); D.alloc(n); parent.alloc(n); idn_change_node.alloc(n); pend.alloc(n);
	S.init(n);

	e = (Edge*)realloc(e, m * sizeof(Edge));

	deleted = (bool*)malloc(m * sizeof(bool)); memset(deleted, 0, m * sizeof(bool));
	d = (int*)malloc(n * sizeof(int)); memset(d, 0, n * sizeof(int));
	r = (int*)malloc(n * sizeof(int)); memset(r, 0, n * sizeof(int));
	c = (int*)malloc(n * sizeof(int)); memset(c, -1, n * sizeof(int));
	dist = (int*)malloc(n * sizeof(int));
	cur = (int*)malloc(n * sizeof(int));
	p = (int*)malloc(n * sizeof(int));
	high_indegree = (int*)malloc(n * sizeof(int)); high_num = 0;
	adj = (int**)malloc(n * sizeof(int*));
	undeg = (int*)malloc(n * sizeof(int));
	int* now_edge_num = (int*)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) undeg[i] = now_edge_num[i] = 0;
	for (int i = 0; i < m; i++) { undeg[e[i].t1]++, undeg[e[i].t2]++; }
	for (int i = 0; i < n; i++) adj[i] = (int*)malloc(undeg[i] * sizeof(int));
	for (int i = 0; i < m; i++)
	{
		int t1 = e[i].t1, t2 = e[i].t2;
		adj[t1][now_edge_num[t1]++] = i;
		adj[t2][now_edge_num[t2]++] = i;
	}
	free(now_edge_num);

	printf("|V|\t\t\t%ld\n|E|\t\t\t%ld\n", n, m);
	return;
}
bool cmp(int a, int b) { return G.d[a] < G.d[b]; }
void Graph::init_orientation()
{
	int* p2node = (int*)malloc(n * sizeof(int)), * node2p = (int*)malloc(n * sizeof(int)), * d_start = (int*)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++)
		d[i] = undeg[i], p2node[i] = i;
	sort(p2node, p2node + n, &cmp);
	int nowr = -1, pointer = 0;
	while (pointer < n)
	{
		node2p[p2node[pointer]] = pointer;
		if (d[p2node[pointer]] > nowr)
			d_start[++nowr] = pointer;
		else
			pointer++;
	}
	nowr = 0;
	for (pointer = 0; pointer < n; pointer++)
	{
		int now = p2node[pointer];
		if (c[now] != -1) continue;
		nowr = max(nowr, d[now]);
		c[now] = nowr;
		for (int j = 0; j < undeg[now]; j++)
		{
			Edge& ne = e[adj[now][j]];
			int tar = ne.t1 == now ? ne.t2 : ne.t1;
			if (c[tar] != -1) continue;
			ne.to = now;
			int lp = d_start[d[now]], rp = node2p[now], ln = p2node[lp], rn = now;
			node2p[ln] = rp;
			node2p[rn] = lp;
			p2node[lp] = rn;
			p2node[rp] = ln;
			d_start[d[now]]++;
			d[now]--;

			lp = d_start[d[tar]], rp = node2p[tar], ln = p2node[lp], rn = tar;
			node2p[ln] = rp;
			node2p[rn] = lp;
			p2node[lp] = rn;
			p2node[rp] = ln;
			d_start[d[tar]]++;
			d[tar]--;
		}
	}
	memset(d, 0, n * sizeof(int));
	for (int i = 0; i < m; i++)
		d[e[i].to]++;
}
int Graph::get_max_d()
{
	int maxd = 0;
	for (int i = 0; i < n; i++)
		maxd = maxd > d[i] ? maxd : d[i];
	return maxd;
}
void Graph::decomposition()
{
	int max_d = get_max_d();
	//printf("Degeneracy: %ld\n", max_d);
	R = (Old_Set*)malloc((max_d + 2) * sizeof(Old_Set));
	for (int i = 0; i <= max_d + 1; i++)
	{
		R[i].reset(n);
		R[i].have_found = false;
	}
	//R[0] is V
	for (int i = 0; i < n; i++)
		R[0].push(i);
	R[0].E_number = m;
	R[0].have_found = true;
	//R[max_d+1] is empty set
	R[max_d + 1].have_found = true;
	//begin recursion
	dfs_layer(max_d + 1, 0);
}
void Graph::dfs_layer(int X_rank, int Y_rank)
{
	Old_Set& X = R[X_rank], & Y = R[Y_rank];
	//get the number of edges
	int XY_E = R[Y_rank].E_number - R[X_rank].E_number;
	//define variables of the binary search
	int l = Y_rank, r = X_rank, l_layer = Y_rank, r_layer = X_rank;
	//first find the R[i], such that i is the maximum number satisfying |XY_E - E(R[i])| < |XY_E|/2
	while (r > l)
	{
		int mid = (r + l + 1) / 2;
		ReTest(mid);
		if (R[Y_rank].E_number - R[mid].E_number < XY_E / 2)
			l = mid, l_layer = mid;
		else
			r = mid - 1, r_layer = mid;
	}
	int k = l;
	//begin deeper recursion
	//cout << "X_rank, k, Y_rank: " << X_rank << ", " << k << ", " << Y_rank << "\n";
	if (k - Y_rank >= 2 && R[Y_rank].size != R[k].size)
		dfs_layer(k, Y_rank);
	k++;
	if (X_rank - k >= 2 && R[X_rank].size != R[k].size)
		dfs_layer(X_rank, k);
	return;
}
void Graph::ReTest(int Z_rank)
{
	//pruning: if R[tem] have been found, then return directly
	if (R[Z_rank].have_found)
		return;
	//find R[Z_rank] within the range R[Y_rank] - R[X_rank], now find the integer X_rank and Y_rank
	int X_rank, Y_rank;
	X_rank = Z_rank + 1;
	while (!R[X_rank].have_found)
		X_rank++;
	Y_rank = Z_rank - 1;
	while (!R[Y_rank].have_found)
		Y_rank--;
	Old_Set& X = R[X_rank], & Y = R[Y_rank];
	test_value = Z_rank - 1;
	while (DinicBFS(X_rank, Y_rank))
	{
		for (int p = 0; p < Y.size; p++)
			cur[Y.nodes[p]] = 0;
		for (int i = 0; i < high_num; i++)
		{
			p[high_indegree[i]] = INF + 1;
			DinicDFS(high_indegree[i]);
		}
	}
	update_r(X_rank, Y_rank, Z_rank);
	R[Z_rank].have_found = true;
	return;
}
bool Graph::DinicBFS(int X_rank, int Y_rank)
{
	int t_dist = INF;
	high_num = 0;
	Old_Set& X = R[X_rank], & Y = R[Y_rank];
	Old_Queue Q(Y.size);
	for (int p = 0; p < Y.size; p++)
	{
		int i = Y.nodes[p];
		dist[i] = INF;
		if (X.in[i]) continue;
		if (d[i] > test_value)
		{
			dist[i] = 0;
			Q.push(i);
			high_indegree[high_num++] = i;
		}
	}
	while (!Q.empty())
	{
		int x = Q.pop();
		for (int i = 0; i < undeg[x]; i++)
		{
			Edge& ne = e[adj[x][i]];
			if (ne.to != x) continue;
			int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (dist[from] != INF || X.in[from]) continue;
			dist[from] = dist[x] + 1;
			if (d[from] < test_value) t_dist = dist[from] + 1;
			Q.push(from);
		}
	}
	return t_dist != INF;
}
bool Graph::DinicDFS(int x)
{
	if (d[x] < test_value)
	{
		d[e[p[x]].to]--;
		d[x]++;
		e[p[x]].to = x;
		return true;
	}
	for (int& i = cur[x]; i < undeg[x]; i++)
	{
		Edge& ne = e[adj[x][i]];
		if (ne.to != x) continue;
		int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
		if (dist[from] != dist[x] + 1) continue;
		p[from] = adj[x][i];
		if (DinicDFS(from))
		{
			if (p[x] >= INF)
			{
				if (d[x] == test_value) return true;
				continue;
			}
			d[e[p[x]].to]--;
			d[x]++;
			e[p[x]].to = x;
			return true;
		}
	}
	return false;
}
void Graph::update_r(int X_rank, int Y_rank, int Z_rank)
{
	Old_Set& X = R[X_rank], & Y = R[Y_rank], & Z = R[Z_rank];
	bool* vis = (bool*)malloc(n * sizeof(bool));
	Old_Queue Q(Y.size);
	for (int p = 0; p < Y.size; p++)
	{
		int i = Y.nodes[p];
		if (X.in[i]) continue;
		vis[i] = false;
		if (d[i] > test_value)
		{
			r[i] = max(r[i], test_value + 1);
			vis[i] = true;
			Q.push(i);
		}
	}
	while (!Q.empty())
	{
		int x = Q.pop();
		Z.push(x);
		for (int i = 0; i < undeg[x]; i++)
		{
			Edge& ne = e[adj[x][i]];
			if (ne.to != x) continue;
			int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (vis[from] || X.in[from]) continue;
			vis[from] = true;
			r[from] = max(r[from], test_value + 1);
			Q.push(from);
		}
	}
	for (int p = 0; p < X.size; p++)
		Z.push(X.nodes[p]);
	Z.E_number = 0;
	for (int p = 0; p < Z.size; p++)
		Z.E_number += d[Z.nodes[p]];
	return;
}
void Graph::display_every_r()
{
	for (int i = 0; i < m; i++)
	{
		int from = e[i].to == e[i].t1 ? e[i].t2 : e[i].t1;
		printf("%ld %ld\n", from, e[i].to);
	}
	for (int i = 0; i < n; i++)
	{
		printf("r[%ld] = %ld, d[%ld] = %ld\n", i, r[i], i, d[i]);
	}
}
void Graph::output_distribution(char* graph_name, FILE* distribution)
{
	fprintf(distribution, "%s,", graph_name);
	int max_d = get_max_d();
	max_d++;
	int* R = (int*)malloc(max_d * sizeof(int)); memset(R, 0, max_d * sizeof(int));
	for (int i = 0; i < n; i++)
		R[r[i]]++;
	for (int i = 0; i < max_d; i++)
		fprintf(distribution, "%ld,", R[i]);
	fprintf(distribution, "\n");
}
bool index_cmp(int a, int b) { return G.r[a] > G.r[b]; }
void Graph::check_correctness()
{
	int now_d = get_max_d();
	while (now_d >= 0)
	{
		memset(dist, 0, n * sizeof(int));
		Old_Queue Q(n);
		for (int i = 0; i < n; i++)
		{
			if (d[i] >= now_d)
			{
				Q.push(i);
				dist[i] = 1;
			}
		}
		while (!Q.empty())
		{
			int x = Q.pop();
			check(r[x] >= now_d, "r[x] error");
			check(d[x] >= now_d - 1, "d[x] error");
			for (int i = 0; i < undeg[x]; i++)
			{
				Edge& ne = e[adj[x][i]];
				if (deleted[adj[x][i]]) continue;
				if (ne.to != x) continue;
				int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (dist[from]) continue;
				dist[from] = 1;
				Q.push(from);
			}
		}
		now_d--;
	}

	for (int i = 0; i < m; i++)
		if (!deleted[i])
			check(edge_label[i] == min(r[e[i].t1], r[e[i].t2]), "edge_label error");

	vis.clear();
	for (int x = 0; x < n; x++)
		if (!vis.in[x])
			vis.insert(x), dfs_correctness(x, -1);

	int sum1 = 0, sum2 = 0, index_m2 = 0;
	for (auto edge_id : index_e) sum1 += edge_label[edge_id], check(!deleted[edge_id], "mst edge deleted error");
	int* sorted_node = (int*)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) sorted_node[i] = i, father[i] = i;
	sort(sorted_node, sorted_node + n, &index_cmp);
	for (int i = 0; i < n; i++) {
		int x = sorted_node[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (deleted[adj[x][j]]) continue;
			if (edge_label[adj[x][j]] != r[x]) continue;
			int y = ne.t1 == x ? ne.t2 : ne.t1;
			int root_x = find_root(x), root_y = find_root(y);
			if (root_x != root_y) {
				father[root_x] = root_y;
				sum2 += edge_label[adj[x][j]];
				index_m2++;
			}
		}
	}
	free(sorted_node);
	check(sum1 == sum2, "sum error");
}
void Graph::dfs_correctness(int x, int fa) {
	for (auto edge_id : index_adj[x]) {
		Edge& ne = e[edge_id];
		int y = ne.t1 == x ? ne.t2 : ne.t1;
		if (fa == y) continue;
		check(!vis.in[y], "mst contain cycle");
		vis.insert(y);
		dfs_correctness(y, x);
	}
}
void Graph::output_dense_subgraph() {
	printf("- %-20s: %d\n", "Number of D nodes", D.size);
	if (false) {
		sort(D.nodes, D.nodes + D.size);
		for (int i = 0; i < D.size; i++) {
			int x = D.nodes[i];
			printf("%d ", x);
		}
		printf("\n");
	}
	return;
}
void Graph::process_deletion() {
	check(!deleted[now_update_edge], "the deleting edge already deleted");
	int u, v;
	v = e[now_update_edge].to, u = e[now_update_edge].t1 == v ? e[now_update_edge].t2 : e[now_update_edge].t1;
	r0 = r[v];
	if (d[v] == r[v] - 1) {
		Q.clear(), parent.clear(); Q.push(v), parent[v] = -1;
		int w = -1;
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (deleted[adj[x][j]]) continue;
				if (ne.to == x) continue;
				if (parent.in[ne.to]) continue;
				if (r[ne.to] != r0) continue;
				Q.push(ne.to), parent[ne.to] = adj[x][j];
				if (d[ne.to] == r0) {
					w = ne.to; break;
				}
			}
			if (w != -1) break;
		}
		check(w != -1, "w equal -1 error");
		Q.clear(), pend.clear(); Q.push(w), pend.insert(w);
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (deleted[adj[x][j]]) continue;
				if (ne.to != x) continue;
				int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (r[from] != r0) continue;
				if (d[from] != r0 - 1) continue;
				if (pend.in[from]) continue;
				pend.insert(from), Q.push(from);
			}
		}
		while (w != v) {
			int edge_id = parent[w];
			d[e[edge_id].to]--;
			e[edge_id].to = e[edge_id].to == e[edge_id].t1 ? e[edge_id].t2 : e[edge_id].t1;
			d[e[edge_id].to]++;
			w = e[edge_id].to;
		}
	}
	else if (d[v] == r[v]) {
		Q.clear(), pend.clear(); Q.push(v), pend.insert(v);
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (deleted[adj[x][j]]) continue;
				if (ne.to != x) continue;
				int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (r[from] != r0) continue;
				if (d[from] != r0 - 1) continue;
				if (pend.in[from]) continue;
				pend.insert(from), Q.push(from);
			}
		}
	}
	deleted[now_update_edge] = true, d[v]--;
	int original_index_m = index_e.size();
	check(bucket[edge_label[now_update_edge]].erase(now_update_edge), "bucket erase error");
	if (index_e.erase(now_update_edge)) {
		// printf("delete edge %d (%d, %d), label %d\n", now_update_edge, e[now_update_edge].t1, e[now_update_edge].t2, edge_label[now_update_edge]);
		check(index_adj[e[now_update_edge].t1].erase(now_update_edge), "index_adj erase error1");
		check(index_adj[e[now_update_edge].t2].erase(now_update_edge), "index_adj erase error2");
	}

	idn_change_node.clear();
	dfs_clock = scc_cnt = 0;
	can_get_r0.clear();
	dfsn.clear();
	scc.clear();
	sccno.clear();
	S.clear();
	for (int i = 0; i < pend.size; i++)
	{
		int now_pend = pend.nodes[i];
		if (sccno.find(now_pend) == sccno.end())
			scc_dfs(now_pend);
		if (can_get_r0.find(now_pend) == can_get_r0.end())
			r[now_pend]--, idn_change_node.insert(now_pend);
	}

	for (int i = 0; i < idn_change_node.size; i++) {
		int x = idn_change_node.nodes[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (deleted[adj[x][j]]) continue;
			check(edge_label[adj[x][j]] <= r0, "edge_label > r0 error");
			if (edge_label[adj[x][j]] == r0) {
				check(bucket[r0].erase(adj[x][j]), "bucket[r0] erase error");
				edge_label[adj[x][j]]--;
				bucket[r0 - 1].insert(adj[x][j]);
				if (index_e.erase(adj[x][j])) {
					// printf("delete edge %d (%d, %d), label %d\n", adj[x][j], e[adj[x][j]].t1, e[adj[x][j]].t2, edge_label[adj[x][j]]);
					check(index_adj[ne.t1].erase(adj[x][j]), "index_adj[ne.t1] erase error");
					check(index_adj[ne.t2].erase(adj[x][j]), "index_adj[ne.t2] erase error");
				}
			}
		}
	}

	for (int x = 0; x < n; x++) father[x] = x;
	for (auto ne_index : index_e) {
		int root1 = find_root(e[ne_index].t1), root2 = find_root(e[ne_index].t2);
		check(root1 != root2, "root error");
		merge(e[ne_index].t1, e[ne_index].t2);
	}

	bool break_loop = false;
	for (int i = r0; i >= 0; i--) {
		for (auto ne_index : bucket[i]) {
			if (index_e.count(ne_index)) continue;
			int root1 = find_root(e[ne_index].t1), root2 = find_root(e[ne_index].t2);
			if (root1 != root2) {
				merge(e[ne_index].t1, e[ne_index].t2);
				index_e.insert(ne_index);
				// printf("insert edge %d (%d, %d), label %d\n", ne_index, e[ne_index].t1, e[ne_index].t2, edge_label[ne_index]);
				index_adj[e[ne_index].t1].insert(ne_index);
				index_adj[e[ne_index].t2].insert(ne_index);
				if (index_e.size() == original_index_m) {
					break_loop = true;
					break;
				}
			}
		}
		if (break_loop)
			break;
	}
}
int Graph::scc_dfs(int now)
{
	if (d[now] == r0)
		can_get_r0.insert(now);
	int lownow = dfsn[now] = ++dfs_clock;
	S.push(now);
	for (int i = 0; i < undeg[now]; i++)
	{
		Edge& ne = e[adj[now][i]];
		if (ne.to == now || deleted[adj[now][i]]) continue;
		if (r[ne.to] != r0) continue;
		if (dfsn.find(ne.to) == dfsn.end())
			lownow = min(lownow, scc_dfs(ne.to));
		else if (sccno.find(ne.to) == sccno.end())
			lownow = min(lownow, dfsn[ne.to]);
		if (can_get_r0.find(ne.to) != can_get_r0.end())
			can_get_r0.insert(now);
	}
	if (lownow == dfsn[now])
	{
		bool this_scc_can_reach_r0 = false;
		scc.resize(scc_cnt + 1);
		while (true)
		{
			int x = S.pop();
			if (can_get_r0.find(x) != can_get_r0.end())
				this_scc_can_reach_r0 = true;
			scc[scc_cnt].push_back(x);
			sccno[x] = scc_cnt;
			if (x == now)
				break;
		}
		if (this_scc_can_reach_r0)
			for (auto i : scc[scc_cnt])
				can_get_r0.insert(i);
		scc_cnt++;
	}
	return lownow;
}
void Graph::process_insertion() {
	check(deleted[now_update_edge], "the inserting edge already inserted");
	int u, v;
	if (r[e[now_update_edge].t1] > r[e[now_update_edge].t2]) u = e[now_update_edge].t1, v = e[now_update_edge].t2;
	else if (d[e[now_update_edge].t1] > r[e[now_update_edge].t2]) u = e[now_update_edge].t1, v = e[now_update_edge].t2;
	else u = e[now_update_edge].t2, v = e[now_update_edge].t1;
	r0 = r[v]; idn_change_node.clear();
	if (d[v] == r0 - 1) {
		e[now_update_edge].to = v, d[v]++, deleted[now_update_edge] = false;
	}
	else {
		e[now_update_edge].to = v, d[v]++, deleted[now_update_edge] = false;
		int w = -1;
		Q.clear(), parent.clear(); Q.push(v), parent[v] = -1;
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (deleted[adj[x][j]]) continue;
				if (ne.to != x) continue;
				int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (r[from] != r0) continue;
				if (parent.in[from]) continue;
				Q.push(from), parent[from] = adj[x][j];
				if (d[from] == r0 - 1) {
					w = from; break;
				}
			}
		}
		if (w != -1) {
			while (w != v) {
				int edge_id = parent[w];
				w = e[edge_id].to;
				d[e[edge_id].to]--;
				e[edge_id].to = e[edge_id].to == e[edge_id].t1 ? e[edge_id].t2 : e[edge_id].t1;
				d[e[edge_id].to]++;
			}
		}
		else {
			r[v]++, idn_change_node.insert(v);
			Q.clear(); Q.push(v);
			while (!Q.empty()) {
				int x = Q.pop();
				for (int j = 0; j < undeg[x]; j++) {
					Edge& ne = e[adj[x][j]];
					if (deleted[adj[x][j]]) continue;
					if (ne.to != x) continue;
					int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
					if (r[from] != r0) continue;
					r[from]++, Q.push(from), idn_change_node.insert(from);
				}
			}
		}
	}

	edge_label[now_update_edge] = min(r[e[now_update_edge].t1], r[e[now_update_edge].t2]);
	bucket[edge_label[now_update_edge]].insert(now_update_edge);
	List<int> label_change_edge; label_change_edge.insert(now_update_edge);
	for (int i = 0; i < idn_change_node.size; i++) {
		int x = idn_change_node.nodes[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (deleted[adj[x][j]]) continue;
			if (min(r[ne.t1], r[ne.t2]) > edge_label[adj[x][j]]) {
				check(min(r[ne.t1], r[ne.t2]) == r0 + 1, "r0 + 1 error");
				edge_label[adj[x][j]]++;
				if (!index_e.count(adj[x][j])) label_change_edge.insert(adj[x][j]);
			}
		}
	}

	for (auto node = label_change_edge.head_node; node != nullptr; node = node->next) {
		int edge_id = node->value;
		int t1 = e[edge_id].t1, t2 = e[edge_id].t2;
		Q.clear(), parent.clear(); Q.push(t1), parent[t1] = -1;
		bool break_loop = false;
		while (!Q.empty()) {
			int x = Q.pop();
			for (auto ne_id : index_adj[x]) {
				Edge& ne = e[ne_id];
				int y = ne.t1 == x ? ne.t2 : ne.t1;
				if (parent.in[y]) continue;
				parent[y] = ne_id, Q.push(y);
				if (y == t2) {
					break_loop = true;
					break;
				}
			}
			if (break_loop) break;
		}
		int min_label = INF, min_label_edge = -1;
		if (break_loop) {
			int x = t2;
			while (x != t1) {
				if (edge_label[parent[x]] < min_label)
					min_label = edge_label[parent[x]], min_label_edge = parent[x];
				x = e[parent[x]].t1 == x ? e[parent[x]].t2 : e[parent[x]].t1;
			}
		}
		if (min_label == INF) {
			index_e.insert(edge_id);
			index_adj[t1].insert(edge_id);
			index_adj[t2].insert(edge_id);
		}
		else if (edge_label[edge_id] > min_label) {
			index_e.erase(min_label_edge);
			index_adj[e[min_label_edge].t1].erase(min_label_edge);
			index_adj[e[min_label_edge].t2].erase(min_label_edge);
			index_e.insert(edge_id);
			index_adj[t1].insert(edge_id);
			index_adj[t2].insert(edge_id);
		}
	}
	
}
void Graph::construct_bucket() {
	edge_label = (int*)malloc(m * sizeof(int));
	for (int i = 0; i < m; i++)
		edge_label[i] = min(r[e[i].t1], r[e[i].t2]);

	bucket = new unordered_set<int>[pseudoarboricity + 1];
	for (int i = 0; i < m; i++)
		bucket[edge_label[i]].insert(i);

	int* sorted_node = (int*)malloc(n * sizeof(int));
	father = (int*)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) sorted_node[i] = i, father[i] = i;
	sort(sorted_node, sorted_node + n, &index_cmp);
	for (int i = 0; i < n; i++) {
		int x = sorted_node[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (edge_label[adj[x][j]] != r[x]) continue;
			int y = ne.t1 == x ? ne.t2 : ne.t1;
			int root_x = find_root(x), root_y = find_root(y);
			if (root_x != root_y) {
				father[root_x] = root_y;
				index_e.insert(adj[x][j]);
			}
		}
	}
	printf("- %-20s: %d\n", "Index edge number", index_e.size());
	free(sorted_node);

	index_adj = new unordered_set<int>[n];
	for (auto edge_id : index_e) {
		// printf("%d %d\n", e[edge_id].t1, e[edge_id].t2);
		index_adj[e[edge_id].t1].insert(edge_id);
		index_adj[e[edge_id].t2].insert(edge_id);
	}
}

int main(int argc, char** argv)
{
	if (argc != 3 && argc != 2) {
	argument_error:
		printf("Usage: ./main <dataset_address> [<update_edge_address>]\n");
		return 0;
	}
	Timer timer;
	char dataset_address[1024], update_edge_address[1024]; strcpy(dataset_address, argv[1]);

	printf("----------Now processing graph: %s----------\n", argv[1]);

	double runtime;
	timer.start(); G.read_edge(argv[1]); timer.end(); runtime = timer.time();
	printf("- %-20s: %lf seconds\n", "Reading edges runtime", runtime);

	timer.start(); G.init_orientation(); timer.end(); double appro_runtime = timer.time();

	timer.start(); G.decomposition(); timer.end(); runtime = timer.time();
	printf("- %-20s: %lf seconds\n", "Initialize orientation", appro_runtime);
	printf("- %-20s: %lf seconds\n", "Get the decomposition", runtime);

	G.pseudoarboricity = G.get_max_d();
	printf("- %-20s: %d\n", "Pseudoarboricity", G.pseudoarboricity);

	if (argc == 3) {
		strcpy(update_edge_address, argv[2]);
		G.read_update_edge_from_file(update_edge_address);
	}
	else {
		G.generate_update_edge(1000);
	}

	timer.start();
	G.construct_bucket();
	timer.end(); double index_runtime = timer.time();

	printf("- %-20s: %lf\n", "Construct index and bucket time", index_runtime);

	double totaldeltime=0.0;
	double totalinstime=0.0;
	for (int i = 0; i < update_edge.size(); i++) {
		now_update_edge = update_edge[i];
		// printf("- %-20s: %d %d\n", "Now delete edge", G.e[now_update_edge].t1, G.e[now_update_edge].t2);
		timer.start();
		G.process_deletion();
		timer.end(); double delete_runtime = timer.time();
		totaldeltime+=delete_runtime;
		printf("- %-20s: %d time %lf\n", "Now delete edge", i, totaldeltime);

		timer.start();
		G.process_insertion();
		timer.end(); double insert_runtime = timer.time();
		totalinstime+=insert_runtime;
		printf("- %-20s: %d time %lf\n", "Now insert edge", i, totalinstime);
	}
	
	// G.check_correctness();

	// double totalinstime=0.0;
	// for (int i = 0; i < update_edge.size(); i++) {
	// 	timer.start();
	// 	now_update_edge = update_edge[i];
	// 	// printf("- %-20s: %d %d\n", "Now insert edge", G.e[now_update_edge].t1, G.e[now_update_edge].t2);
	// 	G.process_insertion();
	// 	timer.end(); double insert_runtime = timer.time();
	// 	totalinstime+=insert_runtime;
	// 	printf("- %-20s: %d time %lf\n", "Now insert edge", i, totalinstime);
	// }
	
	// G.check_correctness();

	// printf("- %-20s: %lf\n", "Total delete time", delete_runtime);
	// printf("- %-20s: %lf\n", "Total insert time", insert_runtime);
	printf("- %-20s: %d\n", "Number of update edges", update_edge.size());

	return 0;
}
