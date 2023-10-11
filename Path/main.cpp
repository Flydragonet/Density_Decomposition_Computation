#include <bits/stdc++.h>
#include <unordered_set>
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
	type* d, *r;
	type** adj;
	type* adj_length;
	type test_value;

	type get_max_d();
	bool ReTest(type test_value);
	type* dist;
	type* p;
	void update_r();

	bool check();
	void display_every_r();

	~Graph()
	{
		free(e); free(d); free(r); free(adj); free(adj_length); free(dist);
	}
};
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
	dist = (type*)malloc(n * sizeof(type));
	p = (type*)malloc(n * sizeof(type));
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
void Graph::init_orientation()
{
	for (type i = 0; i < m; i++)
		e[i].to = e[i].t2, d[e[i].to]++;
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
	type now_rank = 1;
	while (!ReTest(now_rank))
	{
		//printf("%ld\n", now_rank);
		now_rank++;
	}
}
bool Graph::ReTest(type tem)
{
	test_value = tem;
	unordered_set<type> pend;
	for (type i = 0; i < n; i++)
		if (d[i] == tem - 1)
			pend.insert(i);
	while (true)
	{
		Queue Q(n);
		for (type i = 0; i < n; i++)
			p[i] = INF;
		for (auto i : pend)
		{
			Q.push(i);
			p[i] = INF + 1;
		}
		type tar = INF;
		while (!Q.empty())
		{
			type x = Q.pop();
			if (d[x] >= test_value + 1)
			{
				tar = x;
				break;
			}
			for (type i = 0; i < adj_length[x]; i++)
			{
				Edge& ne = e[adj[x][i]];
				if (ne.to == x) continue;
				if (p[ne.to] != INF) continue;
				p[ne.to] = adj[x][i];
				Q.push(ne.to);
			}
		}
		if (tar != INF)
		{
			type x = tar;
			while (p[x] != INF + 1)
			{
				type from = e[p[x]].t1 == x ? e[p[x]].t2 : e[p[x]].t1;
				d[from]++;
				d[x]--;
				e[p[x]].to = from;
				x = from;
			}
			pend.erase(x);
		}
		else
			break;
	}
	update_r();
	return get_max_d() <= test_value;
}
void Graph::update_r()
{
	Queue Q(n);
	for (type i = 0; i < n; i++)
	{
		if (d[i] >= test_value)
		{
			r[i] = test_value;
			Q.push(i);
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
	FILE* dataset = fopen("./dataset.txt", "r"), *distribution = fopen("./distribution.csv", "w");
	char graph_name[100];
	double runtime;
	while (fscanf(dataset, "%s", &graph_name) == 1)
	{
		printf("----------Now processing graph: %s----------\n", graph_name);
		Graph G;
		
		timer_start(); G.read_edge(graph_name); runtime = timer_end();
		printf("Reading edges runtime:\t%.2lf seconds\n", runtime);

		timer_start(); G.init_orientation(); runtime = timer_end();
		printf("Initialize orientation:\t%.2lf seconds\n", runtime);

		timer_start(); G.decomposition(); runtime = timer_end();
		printf("Get the decomposition:\t%.2f seconds\n", runtime);

		G.output_distribution(graph_name, distribution);

		printf("Pseudoarboricity:\t%ld\n", G.get_max_d());
		assert(G.check());
		//G.display_every_r();
	}
	fclose(dataset);
	fclose(distribution);
	return 0;
}
