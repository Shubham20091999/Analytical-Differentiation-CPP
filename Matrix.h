#pragma once
#include<vector>
#include <iostream>

using namespace std;

template<typename T>
class Matrix {
public:
	std::vector<T> mat;
	unsigned int m;//number of rows
	unsigned int n;//number of colunms

public:
	Matrix(unsigned int m, unsigned int n) :
		m(m), n(n) {
		mat = std::vector<T>(m * n);
		for (unsigned int i = 0; i < m * n; i++) {
			mat[i] = T(0);
		}
	}

	virtual void Populate(unsigned int i, unsigned int j, T val) {
		mat[i * n + j] = val;
	}


	virtual T Peek(unsigned int i, unsigned int j) const {
		return mat[i * n + j];
	}

	virtual T& operator()(unsigned int i, unsigned int j) {
		return mat[i * n + j];
	}

	unsigned int nCols() const {
		return n;
	}

	unsigned int nRows() const {
		return m;
	}

	friend Matrix<T> operator *(const Matrix<T>& a, const Matrix<T>& b) {
		if (a.n != b.m) {
			throw exception("Matrix cannot be multiplied\n");
		}
		Matrix<T> ret(a.m, b.n);
		for (unsigned int i = 0; i < a.m; i++) {
			for (unsigned int j = 0; j < b.n; j++) {
				ret(i, j) = T(0);
				for (unsigned int k = 0; k < a.n; k++) {
					ret(i, j) += a.Peek(i, k) * b.Peek(k, j);
				}
			}
		}
		return ret;
	}

	friend Matrix<T> operator +(const Matrix<T>& a, const Matrix<T>& b) {
		if (a.m != b.m or a.n != b.n) {
			throw exception("Matrix cannot be added\n");
		}
		Matrix<T> ret(a.m, b.n);
		for (unsigned int i = 0; i < a.m; i++) {
			for (unsigned int j = 0; j < b.n; j++) {
				ret(i, j) = a.Peek(i, j) + b.Peek(i, j);
			}
		}
		return ret;
	}

	Matrix<T> Transpose() {
		Matrix<T> ret(n, m);
		for (unsigned int i = 0; i < m; i++) {
			for (unsigned int j = 0; j < n; j++) {
				ret(j, i) = Peek(i, j);
			}
		}
		return ret;
	}

	Matrix<T> operator /(T a) {
		Matrix<T> ret(n, m);
		for (unsigned int i = 0; i < m; i++) {
			for (unsigned int j = 0; j < n; j++) {
				ret(i, j) = Peek(i, j) / a;
			}
		}
		return ret;
	}
};

template<typename T>
class Vec :public Matrix <T> {

public:
	Vec(int m) :
		Matrix<T>(m, 1) {

	}


	T magnitude() {
		T ret = T(0);
		for (unsigned int i = 0; i < this->m; i++) {
			ret += pow(this->mat[i], 2);
		}
		return ret;
	}

	void Populate(unsigned int i, T val) {
		this->mat[i] = val;
	}

	T Peek(unsigned int i) const {
		return this->mat[i];
	}

	T& operator()(unsigned int i) {
		return this->mat[i];
	}

	T& operator[](unsigned int i) {
		return this->mat[i];
	}

	unsigned int size() {
		return this->m;
	}
};