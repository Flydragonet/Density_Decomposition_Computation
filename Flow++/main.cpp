#include <bits/stdc++.h>
using namespace std;
typedef int type;
const type INF = 2000000000;

double timer;
inline void timer_start() { timer = clock(); }
inline double timer_end() { double t = (clock() - timer) / CLOCKS_PER_SEC; return t; }
struct Queue
{
	type* node; type head, tail;
	Queue() {}
	Queue(type size) { node = (type*)malloc(size * sizeof(type)); head = tail = 0; }
	~Queue() { free(node); }
	bool empty() { return head == tail; }
	type pop() { return node[head++]; }
	void push(type x) { node[tail++] = x; }
};
struct Set
{
	type* nodes;
	bool* in;
	type size, E_number;
	bool have_found;
	Set() {}
	Set(type sz) { size = 0; nodes = (type*)malloc(sz * sizeof(type)); in = (bool*)malloc(sz * sizeof(type)); memset(in, 0, sz * sizeof(type)); E_number = 0; }
	void reset(type sz) { size = 0; nodes = (type*)malloc(sz * sizeof(type)); in = (bool*)malloc(sz * sizeof(type)); memset(in, 0, sz * sizeof(type)); E_number = 0; }
	void free_memory() { free(nodes), free(in); }
	void push(type x) { nodes[size++] = x; in[x] = true; }
};

struct Edge { type t1, t2, to; };
struct Graph
{
	Graph() {}
	void read_edge(char* graph_name);
	void init_orientation();
	void decomposition();
	void output_distribution(char* graph_name, FILE* distribution);

	Edge* e;
	type n, m;
	type* d, * r, * c;
	type** adj;
	type* adj_length;
	type test_value;

	type get_max_d();
	void dfs_layer(type X_rank, type Y_rank); //Give R[X_rank] and R[Y_rank], find all R between range (Y_rank, X_rank)
	void ReTest(type Z_rank); //find R[Z_rank]
	Set* R;
	type* dist;
	type* cur;
	type* p;
	type* high_indegree, high_num;
	bool DinicBFS(type X_rank, type Y_rank);
	bool DinicDFS(type);
	void update_r(type X_rank, type Y_rank, type Z_rank);

	bool check();
	void display_every_r();

	~Graph()
	{
		free(e); free(d); free(r); free(adj); free(adj_length); free(cur); free(dist); free(p); free(high_indegree); free(c);
	}
};
Graph G;
void Graph::read_edge(char* graph_name)
{
	const type ARRAY_SIZE_INCREMENT = 500000;
	char prefix[100] = "../Graphs/";
	strcat(prefix, graph_name);
	FILE* dataset, * file;
	file = fopen(prefix, "r");

	type Emax = ARRAY_SIZE_INCREMENT;
	n = m = 0;
	e = (Edge*)malloc(Emax * sizeof(Edge));

	char line[200];
	while (fgets(line, 200, file))
	{
		type i = 0;
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
	e = (Edge*)realloc(e, m * sizeof(Edge));

	d = (type*)malloc(n * sizeof(type)); memset(d, 0, n * sizeof(type));
	r = (type*)malloc(n * sizeof(type)); memset(r, 0, n * sizeof(type));
	c = (type*)malloc(n * sizeof(type)); memset(c, -1, n * sizeof(type));
	dist = (type*)malloc(n * sizeof(type));
	cur = (type*)malloc(n * sizeof(type));
	p = (type*)malloc(n * sizeof(type));
	high_indegree = (type*)malloc(n * sizeof(type)); high_num = 0;
	adj = (type**)malloc(n * sizeof(type*));
	adj_length = (type*)malloc(n * sizeof(type));
	type* now_edge_num = (type*)malloc(n * sizeof(type));
	for (type i = 0; i < n; i++) adj_length[i] = now_edge_num[i] = 0;
	for (type i = 0; i < m; i++) { adj_length[e[i].t1]++, adj_length[e[i].t2]++; }
	for (type i = 0; i < n; i++) adj[i] = (type*)malloc(adj_length[i] * sizeof(type));
	for (type i = 0; i < m; i++)
	{
		type t1 = e[i].t1, t2 = e[i].t2;
		adj[t1][now_edge_num[t1]++] = i;
		adj[t2][now_edge_num[t2]++] = i;
	}
	free(now_edge_num);

	printf("|V|\t\t\t%ld\n|E|\t\t\t%ld\n", n, m);
	return;
}
bool cmp(type a, type b) { return G.d[a] < G.d[b]; }
void Graph::init_orientation()
{
	type* p2node = (type*)malloc(n * sizeof(type)), * node2p = (type*)malloc(n * sizeof(type)), * d_start = (type*)malloc(n * sizeof(type));
	for (type i = 0; i < n; i++)
		d[i] = adj_length[i], p2node[i] = i;
	sort(p2node, p2node + n, &cmp);
	type nowr = -1, pointer = 0;
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
		type now = p2node[pointer];
		if (c[now] != -1) continue;
		nowr = max(nowr, d[now]);
		c[now] = nowr;
		for (type j = 0; j < adj_length[now]; j++)
		{
			Edge& ne = e[adj[now][j]];
			type tar = ne.t1 == now ? ne.t2 : ne.t1;
			if (c[tar] != -1) continue;
			ne.to = now;
			type lp = d_start[d[now]], rp = node2p[now], ln = p2node[lp], rn = now;
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
	memset(d, 0, n * sizeof(type));
	for (type i = 0; i < m; i++)
		d[e[i].to]++;
}
type Graph::get_max_d()
{
	type maxd = 0;
	for (type i = 0; i < n; i++)
		maxd = maxd > d[i] ? maxd : d[i];
	return maxd;
}
void Graph::decomposition()
{
	type max_d = get_max_d();
	//printf("Degeneracy: %ld\n", max_d);
	R = (Set*)malloc((max_d + 2) * sizeof(Set));
	for (type i = 0; i <= max_d + 1; i++)
	{
		R[i].reset(n);
		R[i].have_found = false;
	}
	//R[0] is V
	for (type i = 0; i < n; i++)
		R[0].push(i);
	R[0].E_number = m;
	R[0].have_found = true;
	//R[max_d+1] is empty set
	R[max_d + 1].have_found = true;
	//begin recursion
	dfs_layer(max_d + 1, 0);
}
void Graph::dfs_layer(type X_rank, type Y_rank)
{
	Set& X = R[X_rank], & Y = R[Y_rank];
	//get the number of edges
	type XY_E = R[Y_rank].E_number - R[X_rank].E_number;
	//define variables of the binary search
	type l = Y_rank, r = X_rank, l_layer = Y_rank, r_layer = X_rank;
	//first find the R[i], such that i is the maximum number satisfying |XY_E - E(R[i])| < |XY_E|/2
	while (r > l)
	{
		type mid = (r + l + 1) / 2;
		ReTest(mid);
		if (R[Y_rank].E_number - R[mid].E_number < XY_E / 2)
			l = mid, l_layer = mid;
		else
			r = mid - 1, r_layer = mid;
	}
	type k = l;
	//begin deeper recursion
	//cout << "X_rank, k, Y_rank: " << X_rank << ", " << k << ", " << Y_rank << "\n";
	if (k - Y_rank >= 2 && R[Y_rank].size != R[k].size)
		dfs_layer(k, Y_rank);
	k++;
	if (X_rank - k >= 2 && R[X_rank].size != R[k].size)
		dfs_layer(X_rank, k);
	return;
}
void Graph::ReTest(type Z_rank)
{
	//pruning: if R[tem] have been found, then return directly
	if (R[Z_rank].have_found)
		return;
	//find R[Z_rank] within the range R[Y_rank] - R[X_rank], now find the integer X_rank and Y_rank
	type X_rank, Y_rank;
	X_rank = Z_rank + 1;
	while (!R[X_rank].have_found)
		X_rank++;
	Y_rank = Z_rank - 1;
	while (!R[Y_rank].have_found)
		Y_rank--;
	Set& X = R[X_rank], & Y = R[Y_rank];
	test_value = Z_rank - 1;
	while (DinicBFS(X_rank, Y_rank))
	{
		for (type p = 0; p < Y.size; p++)
			cur[Y.nodes[p]] = 0;
		for (type i = 0; i < high_num; i++)
		{
			p[high_indegree[i]] = INF + 1;
			DinicDFS(high_indegree[i]);
		}
	}
	update_r(X_rank, Y_rank, Z_rank);
	R[Z_rank].have_found = true;
	return;
}
bool Graph::DinicBFS(type X_rank, type Y_rank)
{
	type t_dist = INF;
	high_num = 0;
	Set& X = R[X_rank], & Y = R[Y_rank];
	Queue Q(Y.size);
	for (type p = 0; p < Y.size; p++)
	{
		type i = Y.nodes[p];
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
		type x = Q.pop();
		for (type i = 0; i < adj_length[x]; i++)
		{
			Edge& ne = e[adj[x][i]];
			if (ne.to != x) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (dist[from] != INF || X.in[from]) continue;
			dist[from] = dist[x] + 1;
			if (d[from] < test_value) t_dist = dist[from] + 1;
			Q.push(from);
		}
	}
	return t_dist != INF;
}
bool Graph::DinicDFS(type x)
{
	if (d[x] < test_value)
	{
		d[e[p[x]].to]--;
		d[x]++;
		e[p[x]].to = x;
		return true;
	}
	for (type& i = cur[x]; i < adj_length[x]; i++)
	{
		Edge& ne = e[adj[x][i]];
		if (ne.to != x) continue;
		type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
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
void Graph::update_r(type X_rank, type Y_rank, type Z_rank)
{
	Set& X = R[X_rank], & Y = R[Y_rank], & Z = R[Z_rank];
	bool* vis = (bool*)malloc(n * sizeof(bool));
	Queue Q(Y.size);
	for (type p = 0; p < Y.size; p++)
	{
		type i = Y.nodes[p];
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
		type x = Q.pop();
		Z.push(x);
		for (type i = 0; i < adj_length[x]; i++)
		{
			Edge& ne = e[adj[x][i]];
			if (ne.to != x) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (vis[from] || X.in[from]) continue;
			vis[from] = true;
			r[from] = max(r[from], test_value + 1);
			Q.push(from);
		}
	}
	for (type p = 0; p < X.size; p++)
		Z.push(X.nodes[p]);
	Z.E_number = 0;
	for (type p = 0; p < Z.size; p++)
		Z.E_number += d[Z.nodes[p]];
	return;
}
void Graph::display_every_r()
{
	for (type i = 0; i < m; i++)
	{
		type from = e[i].to == e[i].t1 ? e[i].t2 : e[i].t1;
		printf("%ld %ld\n", from, e[i].to);
	}
	for (type i = 0; i < n; i++)
	{
		printf("r[%ld] = %ld, d[%ld] = %ld\n", i, r[i], i, d[i]);
	}
}
void Graph::output_distribution(char* graph_name, FILE* distribution)
{
	fprintf(distribution, "%s,", graph_name);
	type max_d = get_max_d();
	max_d++;
	type* R = (type*)malloc(max_d * sizeof(type)); memset(R, 0, max_d * sizeof(type));
	for (type i = 0; i < n; i++)
		R[r[i]]++;
	for (type i = 0; i < max_d; i++)
		fprintf(distribution, "%ld,", R[i]);
	fprintf(distribution, "\n");
}
bool Graph::check()
{
	type now_d = get_max_d();
	while (now_d >= 0)
	{
		memset(dist, 0, n * sizeof(type));
		Queue Q(n);
		for (type i = 0; i < n; i++)
		{
			if (d[i] >= now_d)
			{
				Q.push(i);
				dist[i] = 1;
			}
		}
		while (!Q.empty())
		{
			type x = Q.pop();
			if (r[x] < now_d || d[x] < now_d - 1)
				return false;
			for (type i = 0; i < adj_length[x]; i++)
			{
				Edge& ne = e[adj[x][i]];
				if (ne.to != x) continue;
				type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (dist[from]) continue;
				dist[from] = 1;
				Q.push(from);
			}
		}
		now_d--;
	}
	return true;
}

int main()
{
	FILE* dataset = fopen("./dataset.txt", "r"), * distribution = fopen("./distribution.csv", "w");
	char graph_name[100];
	double runtime;

	fscanf(dataset, "%s", &graph_name);

	printf("----------Now processing graph: %s----------\n", graph_name);


	timer_start(); G.read_edge(graph_name); runtime = timer_end();
	printf("Reading edges runtime:\t%.2lf seconds\n", runtime);

	//getchar();


	timer_start(); G.init_orientation(); double appro_runtime = timer_end();

	timer_start(); G.decomposition(); runtime = timer_end();
	printf("Initialize orientation:\t%.2lf seconds\n", appro_runtime);
	printf("Get the decomposition:\t%.2f seconds\n", runtime);

	G.output_distribution(graph_name, distribution);

	printf("Pseudoarboricity:\t%ld\n", G.get_max_d());
	//G.display_every_r();
	assert(G.check());

	fclose(dataset);
	fclose(distribution);
	return 0;
}
