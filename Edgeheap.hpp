#pragma once
#include "utils.hpp"
#include "Linear4.hpp"

class Edge
{
public:
	int id; //边序号
	static int count; //边总数
	std::pair<int, int> edge; //边表示
	Vec3 newvt; //生成的新顶点
	double error; //边对应误差
	bool realedge; //是否是实际存在的边(或者是临近点的pair)

	Edge() {}
	Edge(std::pair<int, int> _edge, Vec3 _newvt, double _err, bool _re) : 
		edge(_edge), error(_err), newvt(_newvt), realedge(_re) {}
	Edge(const Edge& e) : id(e.id), edge(e.edge), newvt(e.newvt), error(e.error), realedge(e.realedge) {}
	void setID(int _id) { id = _id; }
	bool operator<(const Edge& e) const { return error > e.error; }
};
int Edge::count = 0;

class EdgeHeap
{
	std::priority_queue<Edge> heap; //边堆集
	std::set<int> lazydel;
	std::map<std::pair<int, int>, int> edge2id;

public:
	EdgeHeap() {}
	int size() { return heap.size(); }
	void pushEdge(Edge e) //assert: edge(u, v) u < v
	{
		edge2id[e.edge] = ++Edge::count;
		e.setID(Edge::count);
		heap.push(e);
	}
	void delEdge(std::pair<int, int> ep) 
	{
		if (ep.first > ep.second) std::swap(ep.first, ep.second);
		lazydel.insert(edge2id[ep]);
	}
	Edge getTop()
	{
		int pop = 0;
		while (lazydel.count(heap.top().id) > 0) heap.pop(), pop++;
		Edge e = heap.top(); heap.pop();
		return e;
	}
};