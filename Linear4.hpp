#pragma once
#include "utils.hpp"
#include "Vec3.hpp"

class Vec4
{
public:
	double vec[4];

	Vec4() { memset(vec, 0, sizeof(double) * 4); }
	Vec4(const Vec4& ov) { for (int i = 0; i < 4; i++) vec[i] = ov.vec[i]; }
	Vec4(double a, double b, double c, double d) { vec[0] = a, vec[1] = b, vec[2] = c, vec[3] = d; }
	Vec4(Vec3 base, double ext = 1.0) : Vec4(base.x, base.y, base.z, ext) {}
	Vec4(double _vec[4]) { for (int i = 0; i < 4; i++) vec[i] = _vec[i]; }
	double norminf() { double nf = 0; for (int i = 0; i < 4; i++) nf += abs(vec[i]); return nf; }
	Vec4 operator+(const Vec4& o)
	{
		double tmp[4];
		for (int i = 0; i < 4; i++)
			tmp[i] = vec[i] + o.vec[i];
		return Vec4(tmp);
	}
	double& operator[](int i) { return vec[i]; }
	void print()
	{
		printf("Vector:");
		for (int i = 0; i < 4; i++)
			printf("%f;", vec[i]);
		printf("\n");
	}
	operator Vec3() { return Vec3(vec[0], vec[1], vec[2]); } //齐次坐标截断
};

class Matrix4
{
public:
	double mat[4][4];

	Matrix4() { memset(mat, 0, sizeof(double) * 16); }
	Matrix4(double _mat[4][4])
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				mat[i][j] = _mat[i][j];
	}
	Matrix4(const Matrix4& om)
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				mat[i][j] = om.mat[i][j];
	}
	Matrix4 operator+(const Matrix4& o)
	{
		double tmp[4][4];
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				tmp[i][j] = mat[i][j] + o.mat[i][j];
		return Matrix4(tmp);
	}
	void print()
	{
		printf("Matrix:\n");
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				printf("%f ", mat[i][j]);
			printf(";");
		}
		printf("\n");
	}
	double* operator[](int i) { return mat[i]; }
};

double quadFormValue(Matrix4 K, Vec4 v) //给定二次型求值
{
	//K.print();
	//v.print();
	double val = 0.0;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			val += K[i][j] * v[i] * v[j];
	}
	return val;
}

Vec4 solve4(Matrix4 A, Vec4 b) //部分主元LU分解解线性方程组
{
	for (int k = 0; k < 3; k++)
	{
		//选主元过程
		int s = k;
		for (int i = k + 1; i < 4; i++)
			if (abs(A[s][k]) < abs(A[i][k])) s = i;
		if (abs(A[s][k]) < eps) return Vec4(); //LU分解中断:主元为0
		if (s != k)
		{
			std::swap(b[s], b[k]);
			for (int j = 0; j < 4; j++)
				std::swap(A[s][j], A[k][j]);
		}

		//LU分解
		for (int i = k + 1; i < 4; i++)
		{
			A[i][k] = A[i][k] / A[k][k];
			for (int j = k + 1; j < 4; j++)
				A[i][j] = A[i][j] - A[i][k] * A[k][j];
		}
	}

	//回代解方程
	Vec4 x, y;
	for (int i = 0; i < 4; i++)
	{
		y[i] = b[i];
		for (int j = 0; j < i; j++)
			y[i] = y[i] - A[i][j] * y[j];
	}
	for (int i = 3; i >= 0; i--)
	{
		x[i] = y[i];
		for (int j = 3; j > i; j--)
			x[i] = x[i] - A[i][j] * x[j];
		x[i] = x[i] / A[i][i];
	}
	return x;
}