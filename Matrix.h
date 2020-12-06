#pragma once
#include<vector>
#include <iostream>

using namespace std;

template<typename T>
class Matrix {
	std::vector<T> mat;
	unsigned int m;//number of rows
	unsigned int n;//number of colunms

public:
	Matrix(unsigned int m, unsigned int n) :
		m(m), n(n) {
		mat = std::vector<T>(m * n);
		for (unsigned int i = 0; i < m * n; i++) {
			mat[i] = T();
		}
	}

	void Populate(unsigned int i, unsigned int j, T val) {
		mat[i * n + j] = val;
	}


	T Peek(unsigned int i, unsigned int j) {
		return mat[i * n + j];
	}

	T& operator()(unsigned int i, unsigned int j = 0) {
		return mat[i * n + j];
	}

	unsigned int nCols() const {
		return n;
	}

	unsigned int nRows() const {
		return m;
	}
};
