#pragma once
#include "Edgeheap.hpp"
#include "Vertexlist.hpp"

bool vcmp(const Vertex& v1, const Vertex& v2)
{
	return v1.crd.x < v2.crd.x;
}

class Simplifier
{
	double simprate, disthres;
	int facecount;
	EdgeHeap* edgeheap;
	VertexList* vtlist;
	//std::set<std::tuple<int, int, int>> real;

public:
	Simplifier(double _rate = 0.1, double _dist = 1e-4) :
		facecount(0), simprate(_rate), disthres(_dist), edgeheap(new EdgeHeap()), vtlist(new VertexList()) {}
	~Simplifier() { delete edgeheap, vtlist; }
	void read(const char* filename) //网格读入
	{
		FILE* fp = fopen(filename, "r+");
		char type = fgetc(fp);
		double x, y, z;
		int id1, id2, id3;
		int lines = 0;
		while (type != EOF)
		{
			fgetc(fp);
			if (type == 'v')
			{
				fscanf(fp, "%lf %lf %lf", &x, &y, &z);
				vtlist->addVertex(Vec3(x, y, z));
			}
			else if (type == 'f')
			{
				facecount++;
				fscanf(fp, "%d %d %d", &id1, &id2, &id3);
				vtlist->addFace(std::make_tuple(id1, id2, id3));
			}
			else
			{
				char tmp[256];
				fgets(tmp, 256, fp);
			}
			type = fgetc(fp);
			while(type != '#' && type != EOF && type != 'v' && type != 'f')
				type = fgetc(fp);
		}
		fclose(fp);
	}
	void write(const char* filename) //网格输出
	{
		FILE* fp = fopen(filename, "w+");
		int count = 0;
		for (auto it = vtlist->vlist.begin() + 1; it != vtlist->vlist.end(); ++it)
		{
			if (it->lazydel) continue;
			it->id = ++count; //此处重写没有被删除的各个点id
			double x, y, z;
			std::tie(x, y, z) = it->crd.unpack();
			fprintf(fp, "v %f %f %f\n", x, y, z);
		}
		std::set<std::tuple<int, int, int>> output;
		for (int i = 1; i <= vtlist->size(); i++)
		{
			if (vtlist->vlist[i].lazydel) continue;
			auto faces = vtlist->getTrifaces(i);
			output.insert(faces.begin(), faces.end());
		}
		for (auto it = output.begin(); it != output.end(); ++it)
		{
			int id1, id2, id3;
			std::tie(id1, id2, id3) = *it;
			//printf("orig(%d, %d, %d)\n", id1, id2, id3);
			fprintf(fp, "f %d %d %d\n", vtlist->vlist[id1].id, vtlist->vlist[id2].id, vtlist->vlist[id3].id);
		}
		//printf("facecount = %d\noutput contains %d\n", facecount, output.size());
		fclose(fp);
	}
	void addEdge(int u, int v, bool isedge) //包装一个添加边的函数
	{
		if (u > v) std::swap(u, v);
		Vec3 newvt; double err;
		auto ep = std::make_pair(u, v);
		std::tie(newvt, err) = vtlist->getMinErrPoint(ep);
		edgeheap->pushEdge(Edge(ep, newvt, err, isedge));
	}
	void replaceVertex(int oldid1, int oldid2, int newid) //替换原来的两个点为新的一个点
	{
		auto adjv1 = vtlist->vlist[oldid1].getAdjVertex();
		auto adjv2 = vtlist->vlist[oldid2].getAdjVertex();
		adjv1.insert(adjv2.begin(), adjv2.end());
		std::set<std::tuple<int, int, int>> newfaces;
		for (auto it = adjv1.begin(); it != adjv1.end(); ++it)
		{
			if (*it == oldid1 || *it == oldid2) continue;
			edgeheap->delEdge(std::make_pair(oldid1, *it));
			edgeheap->delEdge(std::make_pair(oldid2, *it));
			auto faces = vtlist->vlist[*it].replaceFace(oldid1, oldid2, newid);
			newfaces.insert(faces.begin(), faces.end());
		}
		Vertex& newvt = vtlist->vlist[newid];
		for (auto it = newfaces.begin(); it != newfaces.end(); ++it)
			newvt.addFace(*it);
		for (auto it = adjv1.begin(); it != adjv1.end(); ++it)
			if (*it != oldid1 && *it != oldid2)
				addEdge(newid, *it, true);
	}
	void runSimp()
	{
		int tarfacecount = int(facecount * simprate);
		int step = 0;
		while (facecount > tarfacecount)
		{
			Edge se = edgeheap->getTop();
			if (se.realedge) facecount -= 2;
			vtlist->addVertex(se.newvt);
			int newid = vtlist->size();
			int v1n, v2n;
			std::tie(v1n, v2n) = se.edge;
			replaceVertex(v1n, v2n, newid);
			//printf("merge %d %s and %d %s to %d %s\n", 
				//v1n, vtlist->vlist[v1n].crd.print().c_str(), v2n, vtlist->vlist[v2n].crd.print().c_str(), newid, vtlist->vlist[newid].crd.print().c_str());
			vtlist->delVertex(v1n);
			vtlist->delVertex(v2n);
			step++;
			if (step % 100 == 0)
			{
				printf("\rSimplifying: Faces %d, Vertexes %d, HeapSize %d       ",
					facecount, vtlist->vertexcount, edgeheap->size());
				/*auto realvertex = 0;
				std::set<Triface> rtf;
				for (auto it = vtlist->vlist.begin() + 1; it != vtlist->vlist.end(); ++it)
				{
					if (it->lazydel) continue;
					rtf.insert(it->relfaces.begin(), it->relfaces.end());
					realvertex += 1;
				}
				//printf(", Realfaces %d, RealVertexes %d\n", rtf.size(), realvertex);
				assert(rtf.size() == facecount);*/
			}
		}
		printf("\n");
	}
	void buildHeap()
	{
		std::vector<Vertex> cvlist = vtlist->vlist;
		cvlist.erase(cvlist.begin());
		//std::sort(cvlist.begin(), cvlist.end(), vcmp);
		for (int i = 0; i < cvlist.size(); i++)
		{
			if (i % 1000 == 0)
				printf("\rBuilding: Vertex %d / %d", i, cvlist.size());
			Vertex& vt = cvlist[i];
			auto adjv = vtlist->vlist[vt.id].getAdjVertex();
			for (auto it = adjv.begin(); it != adjv.end(); ++it) //推邻接点入堆
			{
				//printf("this = %d, adj = %d\n", vt.id, *it);
				if (vt.id > *it) continue;
				addEdge(vt.id, *it, true);
			}

			//推临近点入堆
			//TODO: 巨大bug
			/*int j = i + 1;
			while (j < cvlist.size() && cvlist[j].crd.x - vt.crd.x < disthres)
			{
				if ((cvlist[j].crd - vt.crd).sqrlen() < sqr(disthres) && adjv.count(cvlist[j].id) == 0)
					addEdge(vt.id, cvlist[j].id, false);
				j++;
			}*/
		}
		printf("\n");
	}
};