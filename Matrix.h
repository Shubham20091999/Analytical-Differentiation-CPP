#pragma once
#include<map>
#include <iostream>

using namespace std;

template<typename T>
class Matrix {
protected:
	std::map<unsigned int, T> mat;
	unsigned int m;
	unsigned int n;

public:
	Matrix(unsigned int m, unsigned int n) :
		m(m), n(n) {}

	T operator()(unsigned int i, unsigned int j) const {
		return Peek(i, j);
	}

	T Peek(unsigned int i, unsigned int j) const {
		if (i >= m) {
			throw IndexOutOfRange(i, m);
		}
		if (j >= n) {
			throw IndexOutOfRange(j, n);
		}
		auto it = mat.find(getIndex(i, j));
		if (it == mat.end()) {
			return T(0);
		}
		return it->second;
	}

	bool operator()(unsigned int i, unsigned int j, const T& val) {
		return Populate(i, j,val);
	}

	bool Populate(unsigned int i, unsigned int j, const T& val) {
		if (i >= m) {
			throw IndexOutOfRange(i, m);
		}
		if (j >= n) {
			throw IndexOutOfRange(j, n);
		}
		auto it = mat.find(getIndex(i, j));
		if (val == 0) {
			if (it != mat.end()) mat.erase(it->first);
			return false;
		}
		else {
			if (it == mat.end()) mat.insert({ getIndex(i,j),val });
			else it->second = val;
			return true;
		}
	}

	friend Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B) {
		if (A.n != B.n or A.m != B.m) {
			throw GenException("Matrices cannot be added");
		}

		Matrix<T> ret(A.m, A.n);

		for (unsigned int i = 0; i < A.m; i++) {
			for (unsigned int j = 0; j < A.n; j++) {
				re(i, j, A(i, j) + B(i, j));
			}
		}

		return ret;
	}

	Matrix<T> operator/(double b) {
		//TODO: or change the value of A matrix and return it
		Matrix<T> ret(m, n);
		for (unsigned int i = 0; i < m; i++) {
			for (unsigned int j = 0; j < n; j++) {
				ret(i, j, this->operator()(i, j) / b);
			}
		}
		return ret;
	}

	Matrix<T> operator+(const Matrix<T>& B) {
		if (n != B.n or m != B.m) {
			throw GenException("Matrices cannot be added");
		}

		Matrix<T> ret(m, n);

		for (unsigned int i = 0; i < m; i++) {
			for (unsigned int j = 0; j < n; j++) {
				ret(i, j, this->operator()(i, j) + B(i, j));
			}
		}
		return ret;
	}

	void removeAll() {
		mat.clear();
	}

	unsigned int getIndex(unsigned int i, unsigned int j) const {
		return i * n + j;
	}
};

template<typename T>
class Vec :public Matrix<T> {
public:
	Vec(unsigned int m) :
		Matrix<T>(m, 1) {

	}

	void operator()(unsigned int i, const T& val) {
		Matrix<T>::Populate(i, 0, val);
	}

	T operator()(unsigned int i) const {
		return  Matrix<T>::Peek(i, 0);
	}

	T magnitude() const {
		T ret(0);
		for (const auto& a : this->mat) {
			ret += pow(a.second, 2);
		}
		return pow(ret, 1.0 / 2.0);
	}

	unsigned int dim()const {
		return this->m;
	}

	static Matrix<T> multT(const Vec<T>& a, const Vec<T>& b) {
		if (a.dim() != b.dim())
			throw GenException("these vectors cannot be multiplied");

		Matrix<T> ret(a.dim(), a.dim());

		for (unsigned int i = 0; i < a.dim(); i++) {
			for (unsigned int j = 0; j < a.dim(); j++) {
				ret(i, j, a(i) * b(j));
			}
		}
		return ret;
	}
};


namespace linearSolvers {
	static Vec<double> GaussElimination(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.dim();
		Vec<double> x(N);
		double temp;
		int i, j, k;
		for (i = 0; i < N; i++)
		{
			temp = fabs(M(i, i));
			k = i;
			for (j = i + 1; j < N; j++)
				if (temp < fabs(M(j, i)))
				{
					temp = fabs(M(j, i));
					k = j;
				}
			if (fabs(M(k, i)) < 0.00001)
			{
				cout << "The matrix is singular: The system has either no solution or infinitely many solution";
				exit(0);
			}
			if (k != i)
			{
				for (j = 0; j < N; j++)
				{
					temp = M(k, j);
					M(k, j, M(i, j));
					M(i, j, temp);
				}
				temp = u(k);
				u(k, u(i));
				u(i, temp);
			}
			for (j = i + 1; j < N; j++)
			{
				temp = M(j, i) / M(i, i);
				for (k = 0; k < N; k++)
					M(j, k, M(j, k) - temp * M(i, k));
				u(j, u(j) - temp * u(i));
			}
		}
		x(N - 1, u(N - 1) / M(N - 1, N - 1));
		for (i = N - 2; i >= 0; i--)
		{
			x(i, u(i));
			for (j = i + 1; j < N; j++)
				x(i, x(i) - M(i, j) * x(j));
			x(i, x(i) / M(i, i));
		}
		return x;
	}

	static Vec<double> LU_Decomposition(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.dim();
		Vec<double> x(N), z(N);
		Matrix<double> l(N, N), u1(N, N);
		int i, j, k;

		for (i = 0; i < N; i++)
			l(i, 0, M(i, 0));
		for (j = 1; j < N; j++)
			u1(0, j, M(0, j) / l(0, 0));
		for (i = 0; i < N; i++)
			u1(i, i, 1);
		for (i = 1; i < N; i++)
			for (j = 1; j < N; j++)
				if (i >= j)
				{
					l(i, j, M(i, j));
					for (k = 0; k <= j - 1; k++) {

						l(i, j, l(i, j) - l(i, k) * u1(k, j));
					}
				}
				else
				{
					u1(i, j, M(i, j));
					for (k = 0; k <= i - 1; k++)
						u1(i, j, u1(i, j) - l(i, k) * u1(k, j));
					u1(i, j, u1(i, j) / l(i, i));
				}

		z(0, u(0) / l(0, 0));
		for (i = 1; i < N; i++)
		{
			z(i, u(i));
			for (j = 0; j <= i - 1; j++)
				z(i, z(i) - l(i, j) * z(j));
			z(i, z(i) / l(i, i));
		}
		x(N - 1, z(N - 1));
		for (i = N - 2; i >= 0; i--)
		{
			x(i, z(i));
			for (j = i + 1; j < N; j++)
				x(i, x(i) - u1(i, j) * x(j));
		}
		return x;
	}

	static Vec<double> TriDiagonal(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.dim();
		Vec<double> x(N);
		vector<double> A(N - 1), B(N), C(N - 1), D(N);
		vector<double> C_star(N - 1), D_star(N);
		int i, j;

		if (N < 2) {
			cout << "Method can't be applied";
			exit(0);
		}
		for (i = 0; i < N - 2; i++)
			for (j = i; j < N - 2; j++)
			{
				if (M(i, j + 2) != 0)
				{
					cout << "Method can't be applied";
					exit(0);
				}

			}
		for (i = N - 1; i > 1; i--)
			for (j = i; j > 1; j--)
			{
				if (M(i, j - 2) != 0)
				{
					cout << "Method can't be applied";
					exit(0);
				}
			}
		for (i = 0; i < N; i++)
		{
			B[i] = M(i, i);
			D[i] = u(i);
		}
		for (i = 1; i < N; i++)
		{
			A[i] = M(i, i - 1);
		}
		for (i = 0; i < N - 1; i++)
		{
			C[i] = M(i, i + 1);
		}

		C_star[0] = C[0] / B[0];
		D_star[0] = D[0] / B[0];
		for (i = 1; i < N; i++)
		{
			C_star[i] = C[i] / (B[i] - A[i] * C_star[i - 1]);
			D_star[i] = (D[i] - A[i] * D_star[i - 1]) / (B[i] - A[i] * C_star[i - 1]);
		}
		x(N - 1, D_star[N - 1]);
		for (i = N - 2; i >= 0; i--)
		{
			x(i, D_star[i] - C_star[i] * x(i + 1));
		}
		return x;
	}

	static Vec<double> Gauss_Jacobi(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.dim();
		Vec<double> x(N), xn(N);
		int i, j, flag;
		double sum, eps = 0.0001;

		flag = 0;
		for (i = 0; i < N; i++)
		{
			sum = 0;
			for (j = 0; j < N; j++)
				if (i != j)
					sum += fabs(M(i, j));
			if (sum > fabs(M(i, i)))
				flag = 1;
		}
		if (flag == 1)
		{
			flag = 0;
			for (j = 0; j < N; j++)
			{
				sum = 0;
				for (i = 0; i < N; i++)
					if (i != j)
						sum += fabs(M(i, j));
				if (sum > fabs(M(j, j)))
					flag = 1;
			}
		}
		if (flag == 1)
		{
			cout << "The co-efficient matrix is not diagonally dominant\n";
			cout << "The Gauss-Jacobi method doesn't converge surely\n";
			exit(0);
		}
		do
		{
			for (i = 0; i < N; i++)
			{
				sum = u(i);
				for (j = 0; j < N; j++)
					if (j != i)
						sum -= M(i, j) * x(j);
				xn(i, sum / M(i, i));
			}
			flag = 0;
			for (i = 0; i < N; i++)
				if (fabs(x(i) - xn(i)) > eps)
					flag = 1;
			if (flag == 1)
				for (i = 0; i < N; i++)
					x(i, xn(i));
		} while (flag == 1);
		return xn;
	}

	static Vec<double> Gauss_Seidal(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.dim();
		int test, i, j;
		double sum, eps = 0.0001;

		Vec<double> x(N), xn(N);

		test = 0;
		for (i = 0; i < N; i++)
		{
			sum = 0;
			for (j = 0; j < N; j++)
				if (i != j)
					sum += fabs(M(i, j));
			if (sum > fabs(M(i, i)))
				test = 1;
		}
		if (test == 1)
		{
			test = 0;
			for (j = 0; j < N; j++)
			{
				sum = 0;
				for (i = 0; i < N; i++)
					if (i != j)
						sum += fabs(M(i, j));
				if (sum > fabs(M(j, j)))
					test = 1;
			}
		}

		if (test == 1)
		{
			cout << "The co-efficient matrix is not diagonally dominant\n";
			cout << "The Gauss-Seidel method doesn't converge surely\n";
			exit(0);
		}
		do
		{
			for (i = 0; i < N; i++)
			{
				sum = u(i);
				for (j = 0; j < N; j++)
				{
					if (j < i)
						sum -= M(i, j) * xn(j);
					else if (j > i)
						sum -= M(i, j) * x(j);
				}
				xn(i, sum / M(i, i));
			}
			test = 0;
			for (i = 0; i < N; i++)
				if (fabs(x(i) - xn(i)) > eps)
					test = 1;
			if (test == 1)
				for (i = 0; i < N; i++)
					x(i, xn(i));
		} while (test == 1);
		return xn;
	}

	static Vec<double> SOR(Matrix<double>& M, Vec<double>& u) {
		Vec<double> x(u.dim()), xn(u.dim());
		unsigned int i, j, flag;
		double sum, eps = 0.001, w = 1;

		do
		{
			for (i = 0; i < u.dim(); i++)
			{
				sum = u(i) * w + M(i, i) * x(i);
				for (j = 0; j < u.dim(); j++)
				{
					if (j < i)
						sum -= M(i, j) * xn(j) * w;
					else if (j >= i)
						sum -= M(i, j) * x(j) * w;
					xn(i, sum / M(i, i));
				}
			}
			flag = 0;
			for (i = 0; i < u.dim(); i++)
				if (fabs(x(i) - xn(i)) > eps)
					flag = 1;
			if (flag == 1)
				for (i = 0; i < u.dim(); i++)
					x(i, xn(i));
		} while (flag == 1);
		return xn;
	}
};