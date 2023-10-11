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
	void restart() { head = tail = 0; }
	void restart(type size) { node = (type*)malloc(size * sizeof(type)); head = tail = 0; }
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
	bool ReTest(type test_value);
	type* dist;
	type* cur;
	type* p;
	type* high_indegree, high_num;
	type* remain; type remain_size;
	bool DinicBFS();
	bool DinicDFS(type);
	void update_r();

	bool check();
	void display_every_r();

	~Graph()
	{
		free(e); free(d); free(r); free(adj); free(adj_length); free(cur); free(dist); free(p); free(high_indegree); free(remain);
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
	remain = (type*)malloc(n * sizeof(type)); remain_size = 0;
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
	for (int i = 0; i < n; i++) remain[i] = i;
	remain_size = n;
	type now_rank = 1;
	while (!ReTest(now_rank))
	{
		now_rank++;
		//printf("%ld\n", now_rank);
	}
}
bool Graph::ReTest(type tem)
{
	test_value = tem;
	while (DinicBFS())
	{
		memset(cur, 0, n * sizeof(type));
		for (type i = 0; i < high_num; i++)
		{
			p[high_indegree[i]] = INF + 1;
			DinicDFS(high_indegree[i]);
		}
	}
	update_r();
	return get_max_d() <= test_value;
}
bool Graph::DinicBFS()
{
	type t_dist = INF;
	high_num = 0;
	Queue Q(remain_size);
	for (type p = 0; p < remain_size; p++)
	{
		type i = remain[p];
		if (d[i] > test_value)
		{
			dist[i] = 0;
			Q.push(i);
			high_indegree[high_num++] = i;
		}
		else
			dist[i] = INF;
	}
	while (!Q.empty())
	{
		type x = Q.pop();
		if (dist[x] >= t_dist) break;
		for (type i = 0; i < adj_length[x]; i++)
		{
			Edge& ne = e[adj[x][i]];
			if (ne.to != x) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (dist[from] != INF) continue;
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
void Graph::update_r()
{
	Queue Q(remain_size);
	for (type p = 0; p < remain_size; p++)
	{
		type i = remain[p];
		if (d[i] >= test_value)
		{
			r[i] = test_value;
			Q.push(i);
		}
	}
	remain_size = 0;
	while (!Q.empty())
	{
		type x = Q.pop();
		remain[remain_size++] = x;
		for (type i = 0; i < adj_length[x]; i++)
		{
			Edge& ne = e[adj[x][i]];
			if (ne.to != x) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (r[from] == test_value) continue;
			r[from] = test_value;
			Q.push(from);
		}
	}
	return;
}
void Graph::display_every_r()
{
	for (type i = 0; i < n; i++)
	{
		printf("r[%ld] = %ld\n", i, r[i]);
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
			if (r[x] < now_d || d[x] < now_d - 1) return false;
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

	fscanf(dataset, "%s", &graph_name) == 1;

	printf("----------Now processing graph: %s----------\n", graph_name);

	timer_start(); G.read_edge(graph_name); runtime = timer_end();
	printf("Reading edges runtime:\t%.2lf seconds\n", runtime);

	timer_start(); G.init_orientation(); double appro_runtime = timer_end();

	timer_start(); G.decomposition(); runtime = timer_end();
	printf("Initialize orientation:\t%.2lf seconds\n", appro_runtime);
	printf("Get the decomposition:\t%.2f seconds\n", runtime);

	G.output_distribution(graph_name, distribution);

	printf("Pseudoarboricity:\t%ld\n", G.get_max_d());
	assert(G.check());
	//G.display_every_r();
	fclose(dataset);
	fclose(distribution);
	return 0;
}
