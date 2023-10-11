#include <bits/stdc++.h>
#include <unordered_set>
#include <unordered_map>
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
struct Stack
{
	type* node; type head;
	Stack() {}
	void init(type size) { node = (type*)malloc(size * sizeof(type)); head = 0; }
	~Stack() { free(node); }
	bool empty() { return head == 0; }
	type pop() { return node[--head]; }
	void push(type x) { node[head++] = x; }
	void clear() { head = 0; }
};
struct Set
{
	type* nodes;
	bool* in;
	type size, rank;
	Set() {}
	Set(type sz) { size = 0; nodes = (type*)malloc(sz * sizeof(type)); in = (bool*)malloc(sz * sizeof(type)); memset(in, 0, sz * sizeof(type)); }
	void free_memory() { free(nodes), free(in); }
};

struct Edge { type t1, t2, to, color; bool deleted; };
struct Graph
{
	Graph() {}
	void read_edge(char* graph_name);
	void init_orientation();
	void decomposition();

	Edge* e;
	type n, m;
	type* d, * r;
	type** adj;
	type* adj_length;
	type test_value;

	type get_max_d();
	void dfs_layer(Set X, Set Y);
	bool ReTest(type test_value, Set X, Set Y, Set& Z);
	type* dist;
	type* cur;
	type* p;
	type* high_indegree, high_num;
	bool DinicBFS(Set X, Set Y);
	bool DinicDFS(type);
	void update_r(Set X, Set Y, Set& Z);

	bool check();
	void display();

	void init_color();
	type delete_edge_number = 10000;
	type* to_delete;
	void generate_delete_edge();
	type get_max_color_edge(type x);
	type r0;
	void deletion();
	void insertion();
	void delete_edge(type id);
	void insert_edge(type id);

	type dfs_clock, scc_cnt;
	type scc_dfs(type now);
	unordered_set<type> can_get_r0;
	unordered_map<type, type> dfsn;
	unordered_map<type, type> sccno;
	vector<vector<type> > scc;
	Stack S;

	~Graph()
	{
		free(e); free(d); free(r); free(adj); free(adj_length); free(cur); free(dist); free(p); free(high_indegree); free(to_delete);
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
		e[m].t1 = e[m].t2 = 0; e[m].deleted = false;
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
	S.init(n);
	e = (Edge*)realloc(e, m * sizeof(Edge));

	d = (type*)malloc(n * sizeof(type)); memset(d, 0, n * sizeof(type));
	r = (type*)malloc(n * sizeof(type)); memset(r, 0, n * sizeof(type));
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
void Graph::init_orientation()
{
	enum how_to_init_orientation { ARBITRARY, APPROXIMATION };
	type your_choice = APPROXIMATION;
	if (your_choice == ARBITRARY)
	{
		for (type i = 0; i < m; i++)
			e[i].to = e[i].t2, d[e[i].to]++;
	}
	else if (your_choice == APPROXIMATION)
	{
		for (type i = 0; i < m; i++)
		{
			Edge& ne = e[i];
			if (d[ne.t1] < d[ne.t2]) { ne.to = ne.t1; d[ne.t1]++; }
			else { ne.to = ne.t2; d[ne.t2]++; }
		}
		type iteration_number = 30;
		for (type t = 0; t < iteration_number; t++)
		{
			for (type i = 0; i < m; i++)
			{
				Edge& ne = e[i];
				type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (d[ne.to] >= d[from] + 2) { d[ne.to]--; d[from]++; ne.to = from; }
			}
		}
	}
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
	Set X(n), Y(n);
	for (type i = 0; i < n; i++)
		Y.nodes[Y.size++] = i, Y.in[i] = true;
	X.rank = INF, Y.rank = 0;
	dfs_layer(X, Y);
}
void Graph::dfs_layer(Set X, Set Y)
{
	if (X.rank <= Y.rank + 1)
		return;
	type XY_E = 0, XY_V = Y.size - X.size;
	for (type p = 0; p < Y.size; p++)
	{
		type i = Y.nodes[p];
		if (X.in[i]) continue;
		for (type j = 0; j < adj_length[i]; j++)
		{
			Edge& ne = e[adj[i][j]];
			type from = ne.t1 == i ? ne.t2 : ne.t1;
			if (Y.in[from] && !X.in[from] && from > i)
				XY_E++; //edge in Y-X
			else if (X.in[from])
				XY_E++; //edge between Y-X and X
		}
	}
	type k = ceil(1.0 * XY_E / XY_V); if (XY_E == 0) k = 0;
	if (k == X.rank) k--;
	if (k == Y.rank) k++;
	Set Z(n); Z.rank = k;
	ReTest(k - 1, X, Y, Z);
	//cout << "k = " << k << '\n';
	//display_every_r();
	if (Z.size == X.size) return;
	if (X.rank - Z.rank >= 2)
		dfs_layer(X, Z);
	if (Z.rank - Y.rank >= 2)
		dfs_layer(Z, Y);
	Z.free_memory();
}
bool Graph::ReTest(type tem, Set X, Set Y, Set& Z)
{
	test_value = tem;
	while (DinicBFS(X, Y))
	{
		for (type p = 0; p < Y.size; p++)
			cur[Y.nodes[p]] = 0;
		for (type i = 0; i < high_num; i++)
		{
			p[high_indegree[i]] = INF + 1;
			DinicDFS(high_indegree[i]);
		}
	}
	update_r(X, Y, Z);
	return get_max_d() <= test_value;
}
bool Graph::DinicBFS(Set X, Set Y)
{
	type t_dist = INF;
	high_num = 0;
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
void Graph::update_r(Set X, Set Y, Set& Z)
{
	Queue Q(Y.size);
	for (type p = 0; p < Y.size; p++)
	{
		type i = Y.nodes[p];
		if (X.in[i]) continue;
		if (d[i] > test_value)
		{
			r[i] = test_value + 1;
			Q.push(i);
		}
	}
	while (!Q.empty())
	{
		type x = Q.pop();
		Z.nodes[Z.size++] = x, Z.in[x] = true;
		for (type i = 0; i < adj_length[x]; i++)
		{
			Edge& ne = e[adj[x][i]];
			if (ne.to != x) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (r[from] == test_value + 1 || X.in[from]) continue;
			r[from] = test_value + 1;
			Q.push(from);
		}
	}
	for (type p = 0; p < X.size; p++)
		Z.nodes[Z.size++] = X.nodes[p], Z.in[X.nodes[p]] = true;
	return;
}
void Graph::init_color()
{
	type* color_cnt = (type*)malloc(n * sizeof(type));
	memset(color_cnt, 0, n * sizeof(type));
	for (type i = 0; i < n; i++)
		for (type j = 0; j < adj_length[i]; j++)
			if (e[adj[i][j]].to == i)
				e[adj[i][j]].color = color_cnt[i]++;
	free(color_cnt);
}
void Graph::generate_delete_edge()
{
	to_delete = (type*)malloc(delete_edge_number * sizeof(type));
	set<type> unique;
	//srand(time(NULL));
	srand(0);
	while (unique.size() < delete_edge_number)
		unique.insert(rand() % m);
	type cnt = 0;
	for (auto i : unique)
		to_delete[cnt++] = i;
}
type Graph::get_max_color_edge(type x)
{
	for (type i = adj_length[x] - 1; i >= 0; i--)
		if (e[adj[x][i]].to == x && e[adj[x][i]].color == d[x] - 1 && !e[adj[x][i]].deleted)
			return adj[x][i];
	assert(0);
}
void Graph::deletion()
{
	for (type i = 0; i < delete_edge_number; i++)
	{
		//if (i % 10 == 0)
			//printf("%d %d (%d, %d)\n", i, to_delete[i], e[to_delete[i]].t1, e[to_delete[i]].t2);
		delete_edge(to_delete[i]);
		//assert(check());
		//display();
	}
}
void Graph::delete_edge(type id)
{
	enum DEL_or_DELpp { DEL, DELpp }; //choose, DEL or DEL++
	DEL_or_DELpp algorithm = DEL;
	type u, v;
	unordered_set<type> pend;
	if (e[id].to == e[id].t1) u = e[id].t2, v = e[id].t1;
	if (e[id].to == e[id].t2) u = e[id].t1, v = e[id].t2;
	r0 = r[v];
	if (d[v] == r[v] - 1)
	{
		type w = -1;
		Queue Q(n); Q.push(v);
		unordered_map<type, type> p; p[v] = -1;
		while (!Q.empty())
		{
			type x = Q.pop();
			if (d[x] == r[v]) { w = x; break; }
			for (type i = 0; i < adj_length[x]; i++)
			{
				Edge& ne = e[adj[x][i]];
				if (ne.to == x || ne.deleted) continue;
				if (r[ne.to] != r[v] || p.find(ne.to) != p.end()) continue;
				p[ne.to] = adj[x][i];
				Q.push(ne.to);
			}
		}
		assert(w != -1);
		if (algorithm == DELpp)
		{
			Queue Q(n);
			Q.push(w); pend.insert(w);
			while (!Q.empty())
			{
				type x = Q.pop();
				for (type i = 0; i < adj_length[x]; i++)
				{
					Edge& ne = e[adj[x][i]];
					if (ne.to != x || ne.deleted) continue;
					type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
					if (r[from] != r0 || d[from] != r0 - 1 || pend.find(from) != pend.end()) continue;
					pend.insert(from);
					Q.push(from);
				}
			}
		}
		e[get_max_color_edge(w)].color = e[p[w]].color;
		while (true)
		{
			type edge_id = p[w];
			d[e[edge_id].to]--;
			e[edge_id].to = e[edge_id].to == e[edge_id].t1 ? e[edge_id].t2 : e[edge_id].t1;
			d[e[edge_id].to]++;
			w = e[edge_id].to;
			if (w == v)
			{
				e[edge_id].color = e[id].color;
				break;
			}
			e[edge_id].color = e[p[w]].color;
		}
	}
	else if (d[v] == r[v])
	{
		e[get_max_color_edge(v)].color = e[id].color;
		if (algorithm == DELpp)
		{
			Queue Q(n);
			Q.push(v); pend.insert(v);
			while (!Q.empty())
			{
				type x = Q.pop();
				for (type i = 0; i < adj_length[x]; i++)
				{
					Edge& ne = e[adj[x][i]];
					if (ne.to != x || ne.deleted) continue;
					type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
					if (r[from] != r0 || d[from] != r0 - 1 || pend.find(from) != pend.end()) continue;
					pend.insert(from);
					Q.push(from);
				}
			}
		}
	}
	e[id].deleted = true; d[v]--;
	//orientation has updated, then update the r[] array
	if (algorithm == DEL) //DEL or DEL++
	{
		Queue Q(n);
		memset(dist, 0, n * sizeof(type));
		for (type i = 0; i < n; i++)
			if (r[i] == r0 && d[i] == r0)
			{
				Q.push(i);
				dist[i] = 1;
			}
		while (!Q.empty())
		{
			type x = Q.pop();
			for (type i = 0; i < adj_length[x]; i++)
			{
				Edge& ne = e[adj[x][i]];
				if (ne.to != x || ne.deleted) continue;
				type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (dist[from] || r[from] != r0) continue;
				dist[from] = 1;
				Q.push(from);
			}
		}
		for (type i = 0; i < n; i++)
			if (r[i] == r0 && !dist[i])
				r[i] = r0 - 1;
	}
	else if (algorithm == DELpp)
	{
		dfs_clock = scc_cnt = 0;
		can_get_r0.clear();
		dfsn.clear();
		scc.clear();
		sccno.clear();
		S.clear();
		for (auto now_pend : pend)
		{
			if (sccno.find(now_pend) == sccno.end())
				scc_dfs(now_pend);
			if (can_get_r0.find(now_pend) == can_get_r0.end())
				r[now_pend]--;
		}
	}
}
type Graph::scc_dfs(type now)
{
	if (d[now] == r0)
		can_get_r0.insert(now);
	type lownow = dfsn[now] = ++dfs_clock;
	S.push(now);
	for (type i = 0; i < adj_length[now]; i++)
	{
		Edge& ne = e[adj[now][i]];
		if (ne.to == now || ne.deleted) continue;
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
			type x = S.pop();
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
void Graph::insertion()
{
	for (type i = delete_edge_number - 1; i >= 0; i--)
	{
		insert_edge(to_delete[i]);
		//assert(check());
		//display();
	}
}
void Graph::insert_edge(type id)
{
	type u, v;
	if (r[e[id].t1] <= r[e[id].t2]) v = e[id].t1, u = e[id].t2;
	else v = e[id].t2, u = e[id].t1;
	if (d[v] == r[v] - 1)
	{
		e[id].color = r[v] - 1;
		e[id].to = v; e[id].deleted = false; d[v]++;
	}
	else if (d[v] == r[v])
	{
		e[id].to = v; e[id].deleted = false; d[v]++;
		type w = -1;
		Queue Q(n); Q.push(v);
		unordered_map<type, type> p; p[v] = -1;
		while (!Q.empty())
		{
			type x = Q.pop();
			if (d[x] == r[v] - 1) { w = x; break; }
			for (type i = 0; i < adj_length[x]; i++)
			{
				Edge& ne = e[adj[x][i]];
				if (ne.to != x || ne.deleted) continue;
				type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (r[from] != r[v] || p.find(from) != p.end()) continue;
				p[from] = adj[x][i];
				Q.push(from);
			}
		}
		if (w != -1)
		{
			type pre_color = e[p[w]].color;
			e[p[w]].color = r[v] - 1;
			while (true)
			{
				type edge_id = p[w];
				w = e[edge_id].to;
				d[e[edge_id].to]--;
				e[edge_id].to = e[edge_id].to == e[edge_id].t1 ? e[edge_id].t2 : e[edge_id].t1;
				d[e[edge_id].to]++;
				if (w == v)
				{
					if (e[id].to == v)
						e[id].color = pre_color;
					break;
				}
				swap(e[p[w]].color, pre_color);
			}
		}
		else
		{
			r0 = r[v];
			r[v]++;
			e[id].color = r[v] - 1;
			Queue Q(n); Q.push(v);
			while (!Q.empty())
			{
				type x = Q.pop();
				for (type i = 0; i < adj_length[x]; i++)
				{
					Edge& ne = e[adj[x][i]];
					if (ne.to != x || ne.deleted) continue;
					type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
					if (r[from] != r0) continue;
					r[from]++;
					Q.push(from);
				}
			}
		}
	}
}
void Graph::display()
{
	for (type i = 0; i < n; i++)
	{
		//printf("r[%ld] = %ld\n", i, r[i]);
		for (type j = 0; j < adj_length[i]; j++)
		{
			if (e[adj[i][j]].to == i && !e[adj[i][j]].deleted)
			{
				type from = e[adj[i][j]].to == e[adj[i][j]].t1 ? e[adj[i][j]].t2 : e[adj[i][j]].t1;
				printf("%ld %ld\n", from, e[adj[i][j]].to);
			}
		}
	}
	printf("\n");
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
				if (ne.to != x || ne.deleted) continue;
				type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (dist[from]) continue;
				dist[from] = 1;
				Q.push(from);
			}
		}
		for (type i = 0; i < n; i++)
			if (r[i] == now_d && dist[i] == 0)
				return false;
		now_d--;
	}
	for (type i = 0; i < n; i++)
	{
		unordered_set<type> s;
		for (type j = 0; j < d[i]; j++)
			s.insert(j);
		for (type j = 0; j < adj_length[i]; j++)
		{
			Edge& ne = e[adj[i][j]];
			if (ne.to != i || ne.deleted) continue;
			if (s.find(ne.color) == s.end())
				return false;
			s.erase(ne.color);
		}
		if (!s.empty())
			return false;
	}
	return true;
}

int main()
{
	FILE* dataset = fopen("./dataset.txt", "r"), * distribution = fopen("./distribution.csv", "w");
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

		printf("Pseudoarboricity:\t%ld\n", G.get_max_d());

		G.generate_delete_edge();
		G.init_color();

		//assert(G.check());

		timer_start(); G.deletion(); runtime = timer_end();
		printf("Deletion time:\t%.2f seconds\n", runtime);

		//assert(G.check());

		timer_start(); G.insertion(); runtime = timer_end();
		printf("Insertion time:\t%.2f seconds\n", runtime);

		assert(G.check());
		//G.display_every_r();
	}
	fclose(dataset);
	fclose(distribution);
	return 0;
}
