#include <bits/stdc++.h>
using namespace std;
const int INF = 2000000000;

struct Timer {
	double start_time, end_time;
	void start() { start_time = clock(); }
	void end() { end_time = clock(); }
	double time() { return (end_time - start_time) / CLOCKS_PER_SEC; }
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

inline int read_number(FILE* in) {
	int x = 0; char ch = 0; while (ch < '0' || ch > '9') ch = fgetc(in); while (ch >= '0' && ch <= '9') { x = x * 10 + (ch - '0'); ch = fgetc(in); } return x;
}
inline void check(bool flag, const char* message) {
	if (!flag) {
		printf("!!!!! CHECK ERROR !!!!!\n");
		printf("Error message: %s\n", message);
		assert(0);
	}
}

enum Algorithm {
	ENUMERATE, BINARY
};
Algorithm algorithm_used;

struct Edge { int u, v, to; };
struct Graph {
	int M, N;
	Edge* e;
	int* undeg, * indeg;
	int** adj;
	void read_graph_from_dataset(char* dataset_address);

	inline bool in_S(int x) { return indeg[x] < test_value; }
	inline bool in_T(int x) { return indeg[x] > test_value; }

	void initialize_orientation();
	int appro_p;

	int answer_k;
	void solve_enumerate();
	void solve_binary();

	int test_value;
	Set<int> D;
	void test();
	Set<int> T_node; Queue<int> Q;
	Map<int> parent; Set<int> vis; Map<int> dist, cur;

	bool DinicBFS();
	bool DinicDFS(int x);

	void output_dense_subgraph();

	void check_correctness();
};
Graph G;

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

void Graph::read_graph_from_dataset(char* dataset_address) {
	FILE* file = fopen(dataset_address, "r");
	check(file != NULL, "Can not open file dataset_address\n");
	M = N = 0; int t1, t2;
	while (fscanf(file, "%d %d\n", &t1, &t2) == 2) {
		M++;
		N = max(N, max(t1, t2));
	}
	N++;
	fclose(file);

	FILE* in = fopen(dataset_address, "r");
	check(in != NULL, "Can not open file dataset_address\n");

	e = (Edge*)malloc(M * sizeof(Edge));
	undeg = (int*)malloc(N * sizeof(int)); indeg = (int*)malloc(N * sizeof(int));
	memset(undeg, 0, N * sizeof(int)); memset(indeg, 0, N * sizeof(int));
	adj = (int**)malloc(N * sizeof(int*));

	T_node.alloc(N); Q.alloc(N); parent.alloc(N); vis.alloc(N); D.alloc(N); dist.alloc(N); cur.alloc(N);

	for (int i = 0; i < M; i++) {
		e[i].u = read_number(in), e[i].v = read_number(in);
		undeg[e[i].u]++, undeg[e[i].v]++;
	}
	for (int i = 0; i < N; i++) {
		adj[i] = (int*)malloc(undeg[i] * sizeof(Edge));
	}
	memset(undeg, 0, N * sizeof(int));
	for (int i = 0; i < M; i++) {
		int u = e[i].u, v = e[i].v;
		adj[u][undeg[u]++] = i;
		adj[v][undeg[v]++] = i;
	}
}
bool cmp(int a, int b) { return G.undeg[a] < G.undeg[b]; }
void Graph::initialize_orientation() {
	int* p2node = (int*)malloc(N * sizeof(int)), * node2p = (int*)malloc(N * sizeof(int)), * d_start = (int*)malloc(N * sizeof(int)), * c = (int*)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++)
		indeg[i] = undeg[i], p2node[i] = i, c[i] = -1;
	sort(p2node, p2node + N, &cmp);
	int nowr = -1, pointer = 0;
	while (pointer < N)
	{
		node2p[p2node[pointer]] = pointer;
		if (indeg[p2node[pointer]] > nowr)
			d_start[++nowr] = pointer;
		else
			pointer++;
	}
	nowr = 0;
	for (pointer = 0; pointer < N; pointer++)
	{
		int now = p2node[pointer];
		if (c[now] != -1) continue;
		nowr = max(nowr, indeg[now]);
		c[now] = nowr;
		appro_p = max(appro_p, c[now]);
		for (int j = 0; j < undeg[now]; j++)
		{
			Edge& ne = e[adj[now][j]];
			int tar = ne.u == now ? ne.v : ne.u;
			if (c[tar] != -1) continue;
			ne.to = now;
			int lp = d_start[indeg[now]], rp = node2p[now], ln = p2node[lp], rn = now;
			node2p[ln] = rp;
			node2p[rn] = lp;
			p2node[lp] = rn;
			p2node[rp] = ln;
			d_start[indeg[now]]++;
			indeg[now]--;

			lp = d_start[indeg[tar]], rp = node2p[tar], ln = p2node[lp], rn = tar;
			node2p[ln] = rp;
			node2p[rn] = lp;
			p2node[lp] = rn;
			p2node[rp] = ln;
			d_start[indeg[tar]]++;
			indeg[tar]--;
		}
	}
	memset(indeg, 0, N * sizeof(int));
	for (int i = 0; i < M; i++)
		indeg[e[i].to]++;
}
void Graph::solve_enumerate() {
	int k;
	for (k = appro_p; k >= 0; k--) {
		test_value = k - 1, test();
		bool all_in = true;
		for (auto x : Query[now_query]) if (!D.in[x]) all_in = false;
		if (!all_in) continue;
		Q.clear(), vis.clear(); set<int> visited_query_node;
		int root_node = *(Query[now_query].begin());
		Q.push(root_node), vis.insert(root_node), visited_query_node.insert(root_node);
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				int y = ne.u == x ? ne.v : ne.u;
				if (!D.in[y]) continue;
				if (vis.in[y]) continue;
				Q.push(y), vis.insert(y);
				if (Query[now_query].count(y)) visited_query_node.insert(y);
			}
		}
		if (visited_query_node.size() == Query[now_query].size()) break;
	}
	answer_k = k;
	if (answer_k >= 0) {
		printf("- %-20s: %d\n", "Layer number", answer_k);
		output_dense_subgraph();
	}
	else {
		printf("- %-20s\n", "Query nodes are not connected");
	}
}
void Graph::solve_binary() {
	int k_l = -1, k_u = appro_p;
	while (k_l < k_u) {
		int k_m = (k_l + k_u + 1) / 2;
		test_value = k_m - 1, test();
		bool all_in = true;
		for (auto x : Query[now_query]) if (!D.in[x]) all_in = false;
		if (!all_in) {
			k_u = k_m - 1;
			continue;
		}
		Q.clear(), vis.clear(); set<int> visited_query_node;
		int root_node = *(Query[now_query].begin());
		Q.push(root_node), vis.insert(root_node), visited_query_node.insert(root_node);
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				int y = ne.u == x ? ne.v : ne.u;
				if (!D.in[y]) continue;
				if (vis.in[y]) continue;
				Q.push(y), vis.insert(y);
				if (Query[now_query].count(y)) visited_query_node.insert(y);
			}
		}
		if (visited_query_node.size() == Query[now_query].size()) k_l = k_m;
		else k_u = k_m - 1;
	}
	answer_k = k_l;
	if (answer_k >= 0) {
		test_value = answer_k - 1, test();
		printf("- %-20s: %d\n", "Layer number", answer_k);
		output_dense_subgraph();
	}
	else {
		printf("- %-20s\n", "Query nodes are not connected");
	}
}
void Graph::test() {
	while (DinicBFS()) {
		for (int i = 0; i < T_node.size; i++) {
			int x = T_node.nodes[i];
			if (in_T(x))
				parent[x] = -2, cur[x] = 0, DinicDFS(x);
		}
	}
	D.clear(), Q.clear();
	for (int x = 0; x < N; x++) {
		if (in_T(x)) {
			Q.push(x), D.insert(x);
		}
	}
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = ne.to == ne.u ? ne.v : ne.u;
			if (D.in[from]) continue;
			Q.push(from), D.insert(from);
		}
	}
	return;
}
bool Graph::DinicBFS() {
	int dist_t = INF;

	Q.clear(), dist.clear(), parent.clear(), cur.clear(), T_node.clear();
	for (int x = 0; x < N; x++) {
		if (in_T(x))
			dist[x] = 1, Q.push(x), T_node.insert(x);
	}

	bool break_loop = false;
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = ne.to == ne.u ? ne.v : ne.u;
			if (in_S(from)) {
				dist_t = dist[x] + 2; break_loop = true; break;
			}
			if (dist.in[from]) continue;
			dist[from] = dist[x] + 1;
			Q.push(from);
		}
		if (break_loop) break;
	}
	return dist_t != INF;
}
bool Graph::DinicDFS(int x) {
	if (in_S(x)) {
		indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
		return true;
	}
	for (int& j = cur[x]; j < undeg[x]; j++) {
		Edge& ne = e[adj[x][j]];
		if (ne.to != x) continue;
		int from = ne.to == ne.u ? ne.v : ne.u;
		if ((dist[from] != dist[x] + 1) && !in_S(from)) continue;
		parent[from] = adj[x][j];
		if (DinicDFS(from)) {
			if (parent[x] == -2) {
				if (indeg[x] == test_value) return true;
				continue;
			}
			indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
			return true;
		}
	}
	return false;
}
void Graph::check_correctness() {
	if (answer_k != -1) {
		test_value = answer_k - 1;
		test();
		for (auto x : Query[now_query]) check(D.in[x], "x is not in D");
		Q.clear(), vis.clear(); set<int> visited_query_node;
		int root_node = *(Query[now_query].begin());
		Q.push(root_node), vis.insert(root_node), visited_query_node.insert(root_node);
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				int y = ne.u == x ? ne.v : ne.u;
				if (!D.in[y]) continue;
				if (vis.in[y]) continue;
				Q.push(y), vis.insert(y);
				if (Query[now_query].count(y)) visited_query_node.insert(y);
			}
		}
		check(visited_query_node.size() == Query[now_query].size(), "correctness error 1");
	}

	test_value = answer_k;
	test();
	bool all_in = true;
	for (auto x : Query[now_query]) if (!D.in[x]) all_in = false;
	if (!all_in) return;
	Q.clear(), vis.clear(); set<int> visited_query_node;
	int root_node = *(Query[now_query].begin());
	Q.push(root_node), vis.insert(root_node), visited_query_node.insert(root_node);
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = ne.u == x ? ne.v : ne.u;
			if (!D.in[y]) continue;
			if (vis.in[y]) continue;
			Q.push(y), vis.insert(y);
			if (Query[now_query].count(y)) visited_query_node.insert(y);
		}
	}
	check(visited_query_node.size() != Query[now_query].size(), "correctness error 2");

	return;
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

int main(int argc, char** argv) {
	if (argc != 4) {
	argument_error:
		printf("Usage: ./main <dataset_address> <query_address> <algorithm>\n");
		printf("algorithm:\n");
		printf("-e: online (enumerate)\n");
		printf("-b: online++ (binary)\n");
		return 0;
	}
	Timer timer;
	char dataset_address[1024], query_address[1024]; strcpy(dataset_address, argv[1]), strcpy(query_address, argv[2]);
	if (strcmp(argv[3], "-e") == 0) algorithm_used = ENUMERATE;
	else if (strcmp(argv[3], "-b") == 0) algorithm_used = BINARY;
	else goto argument_error;

	double runtime;
	printf("----------Now processing %s----------\n", dataset_address);

	timer.start();
	G.read_graph_from_dataset(dataset_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d, %d\n", "|E|, |V|", G.M, G.N);
	printf("- %-20s: %lf\n", "Read graph time", runtime);

	timer.start();
	read_query_from_file(query_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d\n", "Number of queries", Query.size());

	for (now_query = 0; now_query < Query.size(); now_query++) {
		printf("----- Query #%d -----\n", now_query);

		timer.start();
		G.initialize_orientation();
		timer.end(); runtime = timer.time();
		printf("- %-20s: %lf\n", "Initialization time", runtime);

		timer.start();
		if (algorithm_used == ENUMERATE) G.solve_enumerate();
		else if (algorithm_used == BINARY) G.solve_binary();
		timer.end(); runtime = timer.time();
		printf("- %-20s: %lf\n", "Get D time", runtime);

		G.check_correctness();
	}
	return 0;
}
