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
	type size, E_number;
	bool have_found;
	Set() {}
	Set(type sz) { size = 0; nodes = (type*)malloc(sz * sizeof(type)); in = (bool*)malloc(sz * sizeof(type)); memset(in, 0, sz * sizeof(type)); E_number = 0; }
	void reset(type sz) { size = 0; nodes = (type*)malloc(sz * sizeof(type)); in = (bool*)malloc(sz * sizeof(type)); memset(in, 0, sz * sizeof(type)); E_number = 0; }
	void free_memory() { free(nodes), free(in); }
	void push(type x) { nodes[size++] = x; in[x] = true; }
};

struct Edge { type t1, t2, to; bool deleted; };
struct Graph
{
	Graph() {}
	void read_edge(char* graph_name);
	void init_orientation();
	void decomposition();

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
	void display();

	type delete_edge_number = 10000;
	type* to_delete;
	void generate_delete_edge();
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
	DEL_or_DELpp algorithm = DELpp;
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
		while (true)
		{
			type edge_id = p[w];
			d[e[edge_id].to]--;
			e[edge_id].to = e[edge_id].to == e[edge_id].t1 ? e[edge_id].t2 : e[edge_id].t1;
			d[e[edge_id].to]++;
			w = e[edge_id].to;
			if (w == v)
				break;
		}
	}
	else if (d[v] == r[v])
	{
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
			while (true)
			{
				type edge_id = p[w];
				w = e[edge_id].to;
				d[e[edge_id].to]--;
				e[edge_id].to = e[edge_id].to == e[edge_id].t1 ? e[edge_id].t2 : e[edge_id].t1;
				d[e[edge_id].to]++;
				if (w == v)
					break;
			}
		}
		else
		{
			r0 = r[v];
			r[v]++;
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

	timer_start(); G.init_orientation(); runtime = timer_end();
	printf("Initialize orientation:\t%.2lf seconds\n", runtime);

	timer_start(); G.decomposition(); runtime = timer_end();
	printf("Get the decomposition:\t%.2f seconds\n", runtime);

	printf("Pseudoarboricity:\t%ld\n", G.get_max_d());

	G.generate_delete_edge();

	//assert(G.check());

	timer_start(); G.deletion(); runtime = timer_end();
	printf("Deletion time:\t%.2f seconds\n", runtime);

	//assert(G.check());

	timer_start(); G.insertion(); runtime = timer_end();
	printf("Insertion time:\t%.2f seconds\n", runtime);

	assert(G.check());
	//G.display_every_r();

	fclose(dataset);
	fclose(distribution);
	return 0;
}
