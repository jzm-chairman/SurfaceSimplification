#pragma once
#include "utils.hpp"
#include "Vec3.hpp"
#include "Linear4.hpp"

#define DISTPUN 0.01
enum borderType {UNDETERMINE, BORDER, INNER};

class Triface
{
public:
	int pid[3];
	Triface(int id1, int id2, int id3) { pid[0] = id1, pid[1] = id2, pid[2] = id3; }
	Triface(std::tuple<int, int, int> pt = std::tuple<int, int, int>()) { int id1, id2, id3; std::tie(id1, id2, id3) = pt; pid[0] = id1, pid[1] = id2, pid[2] = id3; }
	Triface(const Triface& t) { for (int i = 0; i < 3; i++) pid[i] = t.pid[i]; }
	~Triface() {}
	auto toTuple() { return std::make_tuple(pid[0], pid[1], pid[2]); }
	bool searchByID(int id1, int id2) //delete face if containing 2 old vertexes
	{
		int cnt = 0;
		for (int i = 0; i < 3; i++)
			if (pid[i] == id1 || pid[i] == id2) cnt++;
		return cnt >= 2;
	}
	bool replace(int oldid, int newid) //do nothing if oldid not exist and return true if oldid exist
	{
		for (int i = 0; i < 3; i++)
			if (pid[i] == oldid)
			{
				//printf("modify (%d, %d, %d) ", pid[0], pid[1], pid[2]);
				pid[i] = newid;
				//printf("to (%d, %d, %d)\n", pid[0], pid[1], pid[2]);
				return true;
			}
		return false;
	}
	bool operator<(const Triface& t) const { return pid[0] == t.pid[0] ? (pid[1] == t.pid[1] ? pid[2] < t.pid[2] : pid[1] < t.pid[1]) : pid[0] < t.pid[0]; }
	bool operator==(const Triface& t) const { return pid[0] == t.pid[0] && pid[1] == t.pid[1] && pid[2] == t.pid[2]; }
	void print() { printf("face(%d, %d, %d)\n", pid[0], pid[1], pid[2]); }
};

class Vertex
{
public:
	int id; //顶点标准序号(输出前等于数组下标)
	Vec3 crd; //顶点坐标
	//std::set<int> nbr; //邻接表
	std::vector<Triface> relfaces;
	borderType btype; //记录是否是边界
	bool lazydel;

	Vertex(int _id, Vec3 _crd = Vec3()) : id(_id), crd(_crd), lazydel(false), btype(UNDETERMINE) {}
	Vertex(const Vertex& ov) : id(ov.id), crd(ov.crd), relfaces(ov.relfaces), lazydel(ov.lazydel) {}
	void del() { lazydel = true; }
	void addFace(std::tuple<int, int, int> face) { relfaces.push_back(Triface(face)); }
	auto replaceFace(int oldid1, int oldid2, int newid) // replace vertex(oldid1, oldid2) by vertex(newid)
	{ // return a set including all faces with Vertex(newid)
		//printf("id(%d) replace(%d, %d) by %d\n", id, oldid1, oldid2, newid);
		std::set<std::tuple<int, int, int>> newfaces;
		for (auto it = relfaces.begin(); it != relfaces.end(); )
		{
			if (it->searchByID(oldid1, oldid2))
			{
				//printf("id(%d) remove ", id); it->print();
				it = relfaces.erase(it);
			}
			else
			{
				if (it->replace(oldid1, newid) || it->replace(oldid2, newid))
					newfaces.insert(it->toTuple());
				++it;
			}
		}
		return newfaces;
	}
	std::set<int> getAdjVertex() //取所有邻接点
	{
		std::set<int> adjv;
		//printf("id(%d), relfacecount = %d\n", id, relfaces.size());
		for (auto it = relfaces.begin(); it != relfaces.end(); ++it)
		{
			//it->print();
			for (int i = 0; i < 3; i++)
				if (id != it->pid[i]) adjv.insert(it->pid[i]);
		}
		return adjv;
	}
};

class VertexList
{
public:
	std::vector<Vertex> vlist;
	static int vertexcount;

	VertexList() : vlist(std::vector<Vertex>()) { vlist.push_back(Vertex(0)); }
	int size() { return vlist.size() - 1; }
	void addVertex(Vec3 pos) { vlist.push_back(Vertex(vlist.size(), pos)); vertexcount++; }
	void delVertex(int n) { vlist[n].del(); vertexcount--; }
	void addFace(std::tuple<int, int, int> face)
	{
		int id1, id2, id3;
		std::tie(id1, id2, id3) = face;
		vlist[id1].addFace(face);
		vlist[id2].addFace(face);
		vlist[id3].addFace(face);
	}
	auto getTrifaces(int u) //返回值为数组下标的vector
	{
		std::vector<std::tuple<int, int, int>> faces;
		auto origfaces = vlist[u].relfaces;
		for (auto it = origfaces.begin(); it != origfaces.end(); ++it)
			faces.push_back(it->toTuple());
		return faces;
	}
	Matrix4 quadErrMat(int u) 
	//指定下标求基本误差矩阵并返回对应面片以备其他操作
	{
		Matrix4 err;
		auto faces = getTrifaces(u);
		for (auto it = faces.begin(); it != faces.end(); ++it)
		{
			int u1, u2, u3;
			std::tie(u1, u2, u3) = *it;
			Vec3& v1c = vlist[u1].crd;
			Vec3& v2c = vlist[u2].crd;
			Vec3& v3c = vlist[u3].crd;
			//叉积法求解平面法方向(a, b, c)
			Vec3 normdir = (v2c - v1c).cross(v3c - v1c); 
			if(normdir.isNormal()) normdir.norm();
			double d = -v1c.dot(normdir); //平面的参数d
			Vec4 plane(normdir, d);
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					err[i][j] += plane[i] * plane[j]; //K=p*p'
			/*if (!err.isNormal())
			{
				printf("id(%d): %s\n", u1, v1c.print().c_str());
				printf("id(%d): %s\n", u2, v2c.print().c_str());
				printf("id(%d): %s\n", u3, v3c.print().c_str());
				printf("Plane:"); plane.print();
				err.print();
			}
			assert(err.isNormal());*/
		}
		return err;
	}
	int checkIntersec(int u, int v)
	{
		auto faces1 = vlist[u].relfaces;
		auto faces2 = vlist[v].relfaces;
		std::vector<Triface> itsc(10);
		std::sort(faces1.begin(), faces1.end());
		std::sort(faces2.begin(), faces2.end());
		auto it = std::set_intersection(faces1.begin(), faces1.end(), faces2.begin(), faces2.end(), itsc.begin());
		return it - itsc.begin();
	}
	std::tuple<Vec3, double> getMinErrPoint(std::pair<int, int> vp) //指定一对点求解合适的新顶点位置和对应误差
	{
		int u, v; std::tie(u, v) = vp;
		//assert(vlist[u].btype != UNDETERMINE && vlist[v].btype != UNDETERMINE);
		//if (vlist[u].btype == BORDER || vlist[v].btype == BORDER) return std::make_tuple(Vec3(), INF);
		if (checkIntersec(u, v) < 2) return std::make_tuple(Vec3(), INF);
		Vec3 uc = vlist[u].crd; Vec3 vc = vlist[v].crd;
		double distpun = sqrt((vlist[u].crd - vlist[v].crd).sqrlen()) * DISTPUN;
		Matrix4 errmatu, errmatv;
		errmatu = quadErrMat(u);
		errmatv = quadErrMat(v);
		Matrix4 errmat = errmatu + errmatv;
		Matrix4 A = errmat;
		A[3][0] = A[3][1] = A[3][2] = 0, A[3][3] = 1; //误差矩阵转换为系数矩阵
		Vec4 solu = solve4(A, Vec4(Vec3())); //注意这几步有隐式转换
		double suvdist = sqrt((Vec3(solu) - uc).sqrlen()) + sqrt((Vec3(solu) - vc).sqrlen());
		double uvdist = sqrt((uc - vc).sqrlen()); //s到uv距离和不能超过uv距离太多防止出现不合理解
		double errval = quadFormValue(errmat, solu);
		//if (solu.norminf() > eps)
			//printf("(%d, %d) valid, error = (%f, %f), minerror = %f\n", vp.first, vp.second, quadFormValue(errmat, vlist[u].crd), quadFormValue(errmat, vlist[v].crd), errval);
		if (solu.norminf() > eps && suvdist < 2 * uvdist) //解有效，选取即可
			return std::make_tuple(solu, errval + distpun);
		//系数矩阵不可逆的时候用简单策略
		//errmat.print();
		//printf("%d ", u); Vec4(vlist[u].crd).print();
		//printf("%d ", v); Vec4(vlist[v].crd).print();
		double erru = quadFormValue(errmat, vlist[u].crd);
		double errv = quadFormValue(errmat, vlist[v].crd);
		Vec3 mid = (vlist[u].crd + vlist[v].crd) / 2;
		double errmid = quadFormValue(errmat, mid);
		//printf("(%d, %d) invalid, error = (%f, %f, %f)\n", vp.first, vp.second, erru, errv, errmid);
		if (errmid < erru && errmid < errv) return std::make_tuple(mid, errmid + distpun);
		else if (erru < errv) return std::make_tuple(vlist[u].crd, erru + distpun);
		else return std::make_tuple(vlist[v].crd, errv + distpun);
	}
};
int VertexList::vertexcount = 0;