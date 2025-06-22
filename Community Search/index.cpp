#include <bits/stdc++.h>
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

	bool check_correctness();
	void display_every_r();

	int answer_k;
	Queue<int> Q; Set<int> vis; Set<int> D; Map<int> parent;
	int* edge_label;
	void construct_basic_index();
	void basic_query();
	void construct_improve_index();
	void improve_query();
	void output_dense_subgraph();
	int index_m, * index_e, ** index_adj, * index_undeg, * father;
	int find_root(int x) {
		if (father[x] == x) return x; return father[x] = find_root(father[x]);
	}
	void merge(int x, int y) {
		int root_x = find_root(x), root_y = find_root(y);
		father[root_x] = root_y;
	}

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

vector<set<int> > Query; int now_query;
void read_query_from_file(char* query_address) {
	FILE* file = fopen(query_address, "r");
	check(file != NULL, "Can not open file query_address\n");
	int MAX_LINE_LENGTH = 131072, N;
	char line[131072];
	fgets(line, MAX_LINE_LENGTH, file);
	sscanf(line, "%d", &N);
	Query.resize(N);
	for (int i = 0; i < N; i++) {
		if (fgets(line, MAX_LINE_LENGTH, file) != NULL) {
			int num;
			char* ptr = line;

			while (sscanf(ptr, "%d", &num) == 1) {
				Query[i].insert(num);
				while (*ptr != ' ' && *ptr != '\n' && *ptr != '\0') {
					ptr++;
				}
				while (*ptr == ' ') {
					ptr++;
				}
			}
		}
		else {
			check(0, "query_file number of lines error\n");
			break;
		}
	}
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
	Q.alloc(n); vis.alloc(n); D.alloc(n); parent.alloc(n);

	e = (Edge*)realloc(e, m * sizeof(Edge));

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
bool Graph::check_correctness()
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
			if (r[x] < now_d || d[x] < now_d - 1)
				return false;
			for (int i = 0; i < undeg[x]; i++)
			{
				Edge& ne = e[adj[x][i]];
				if (ne.to != x) continue;
				int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (dist[from]) continue;
				dist[from] = 1;
				Q.push(from);
			}
		}
		now_d--;
	}
	return true;
}
void Graph::construct_basic_index() {
	edge_label = (int*)malloc(m * sizeof(int));
	for (int i = 0; i < m; i++)
		edge_label[i] = min(r[e[i].t1], r[e[i].t2]);
}
void Graph::basic_query() {
	vector<queue<int> > queue; queue.resize(pseudoarboricity + 1);
	vis.clear(); set<int> visited_query_node;
	int now_k = INF; for (auto x : Query[now_query]) now_k = min(now_k, r[x]);
	int root_node = *(Query[now_query].begin());
	queue[now_k].push(root_node), vis.insert(root_node), visited_query_node.insert(root_node);
	while (true) {
		while (!queue[now_k].empty()) {
			int x = queue[now_k].front(); queue[now_k].pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				int y = ne.t1 == x ? ne.t2 : ne.t1;
				if (vis.in[y]) continue;
				queue[min(edge_label[adj[x][j]], now_k)].push(y), vis.insert(y);
				if (Query[now_query].count(y)) visited_query_node.insert(y);
			}
		}
		if (visited_query_node.size() == Query[now_query].size()) break;
		now_k--; if (now_k < 0) break;
	}
	answer_k = now_k;
	if (answer_k >= 0) {
		// printf("- %-20s: %d\n", "Layer number", answer_k);
		D.clear();
		for (int x = 0; x < n; x++) if (r[x] >= answer_k) D.insert(x);
		// output_dense_subgraph();
	}
	else {
		// printf("- %-20s\n", "Query nodes are not connected");
	}
}
bool index_cmp(int a, int b) { return G.r[a] > G.r[b]; }
void Graph::construct_improve_index() {
	edge_label = (int*)malloc(m * sizeof(int));
	for (int i = 0; i < m; i++)
		edge_label[i] = min(r[e[i].t1], r[e[i].t2]);

	int* sorted_node = (int*)malloc(n * sizeof(int));
	father = (int*)malloc(n * sizeof(int)), index_e = (int*)malloc((n - 1) * sizeof(int)), index_undeg = (int*)malloc(n * sizeof(int)), memset(index_undeg, 0, n * sizeof(int)), index_adj = (int**)malloc(n * sizeof(int*));
	index_m = 0;
	for (int i = 0; i < n; i++) sorted_node[i] = i, father[i] = i;
	sort(sorted_node, sorted_node + n, &index_cmp);
	int sum = 0;
	for (int i = 0; i < n; i++) {
		int x = sorted_node[i];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = ne.t1 == x ? ne.t2 : ne.t1;
			if (r[y] < r[x]) continue;
			int root_x = find_root(x), root_y = find_root(y);
			if (root_x != root_y) {
				father[root_x] = root_y;
				index_e[index_m++] = adj[x][j];
				sum += min(r[x], r[y]);
			}
		}
	}
	printf("%d\n", sum);
	printf("- %-20s: %d\n", "Index edge number", index_m);
	for (int i = 0; i < index_m; i++) index_undeg[e[index_e[i]].t1]++, index_undeg[e[index_e[i]].t2]++;
	for (int x = 0; x < n; x++) index_adj[x] = (int*)malloc(index_undeg[x] * sizeof(int));
	int* now_edge_num = (int*)malloc(n * sizeof(int)); memset(now_edge_num, 0, n * sizeof(int));
	for (int i = 0; i < index_m; i++) {
		int t1 = e[index_e[i]].t1, t2 = e[index_e[i]].t2;
		index_adj[t1][now_edge_num[t1]++] = index_e[i];
		index_adj[t2][now_edge_num[t2]++] = index_e[i];
	}
	return;
}
void Graph::improve_query() {
	parent.clear(); Q.clear();
	int root_node = *(Query[now_query].begin()); set<int> visited_query_node;
	Q.push(root_node), parent[root_node] = -1, visited_query_node.insert(root_node);
	while (!Q.empty()) {
		int x = Q.pop();
		for (int i = 0; i < index_undeg[x]; i++) {
			Edge& ne = e[index_adj[x][i]];
			int y = ne.t1 == x ? ne.t2 : ne.t1;
			if (parent.in[y]) continue;
			parent[y] = index_adj[x][i], Q.push(y);
			if (Query[now_query].count(y)) {
				visited_query_node.insert(y);
				if (visited_query_node.size() == Query[now_query].size())
					goto find_tree_finished;
			}
		}
	}
find_tree_finished:
	if (visited_query_node.size() != Query[now_query].size()) {
		answer_k = -1;
		printf("- %-20s\n", "Query nodes are not connected");
		return;
	}
	Q.clear(), vis.clear();
	int now_k = INF;
	for (auto x : Query[now_query]) Q.push(x), vis.insert(x);
	while (!Q.empty()) {
		int x = Q.pop();
		if (parent[x] == -1) continue;
		now_k = min(now_k, edge_label[parent[x]]);
		int parent_node = e[parent[x]].t1 == x ? e[parent[x]].t2 : e[parent[x]].t1;
		if (vis.in[parent_node]) continue;
		Q.push(parent_node), vis.insert(parent_node);
	}
	answer_k = now_k;
	if (answer_k >= 0) {
		// printf("- %-20s: %d\n", "Layer number", answer_k);
		D.clear();
		for (int x = 0; x < n; x++) if (r[x] >= answer_k) D.insert(x);
		// output_dense_subgraph();
	}
	else {
		// printf("- %-20s\n", "Query nodes are not connected");
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

int main(int argc, char** argv)
{
	if (argc != 4) {
	argument_error:
		printf("Usage: ./main <dataset_address> <query_address> <algorithm>\n");
		printf("algorithm:\n");
		printf("-b: basic_search\n");
		printf("-i: improve_search\n");
		return 0;
	}
	Timer timer;
	char dataset_address[1024], query_address[1024]; strcpy(dataset_address, argv[1]), strcpy(query_address, argv[2]);
	if (strcmp(argv[3], "-b") == 0) algorithm_used = BASIC;
	else if (strcmp(argv[3], "-i") == 0) algorithm_used = IMPROVE;
	else goto argument_error;

	printf("----------Now processing graph: %s----------\n", argv[1]);

	double runtime;
	timer.start(); G.read_edge(argv[1]); timer.end(); runtime = timer.time();
	printf("- %-20s: %lf seconds\n", "Reading edges runtime", runtime);

	timer.start(); G.init_orientation(); timer.end(); double appro_runtime = timer.time();

	timer.start(); G.decomposition(); timer.end(); runtime = timer.time();
	// printf("- %-20s: %lf seconds\n", "Initialize orientation", appro_runtime);
	// printf("- %-20s: %lf seconds\n", "Get the decomposition", runtime);

	double indexbuildtime=0.0;
	indexbuildtime+=appro_runtime;
	indexbuildtime+=runtime;

	G.pseudoarboricity = G.get_max_d();
	printf("- %-20s: %d\n", "Pseudoarboricity", G.pseudoarboricity);

	timer.start();
	if (algorithm_used == BASIC) G.construct_basic_index();
	else if (algorithm_used == IMPROVE) G.construct_improve_index();
	timer.end(); runtime = timer.time();
	indexbuildtime+=runtime;
	printf("- %-20s: %lf\n", "Index construct time", indexbuildtime);

	timer.start();
	read_query_from_file(query_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d\n", "Number of queries", Query.size());

	double total_runtime=0.0;
	for (now_query = 0; now_query < Query.size(); now_query++) {
		printf("----- Query #%d -----\n", now_query);
		timer.start();
		if (algorithm_used == BASIC) G.basic_query();
		else if (algorithm_used == IMPROVE) G.improve_query();
		timer.end(); runtime = timer.time();
		printf("- %-20s: %lf\n", "Get D time", runtime);
		total_runtime+=runtime;
		printf("- %-20s: %lf\n", "the total time", total_runtime);
		// G.check_correctness();
	}
	// timer.end(); runtime = timer.time();
	// printf("- %-20s: %lf\n", "Query time", runtime);


	return 0;
}
