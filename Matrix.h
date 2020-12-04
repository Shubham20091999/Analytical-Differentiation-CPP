#pragma once
#include<vector>

template<typename T>
class Matrix {
	std::vector<T> mat;
	unsigned int m;//number of rows
	unsigned int n;//number of colunms

public:
	Matrix(unsigned int m, unsigned int n) :
		m(m), n(n) {
		mat = std::vector<T>(m * n);
	}

	void Populate(unsigned int i, unsigned int j, T val) {
		mat[i * m + j] = val;
	}

	T& Glimpse(unsigned int i, unsigned int j){
		return mat[i * m + j];
	}
	unsigned int nCols() const {
		return m;
	}

	unsigned int nRows() const {
		return n;
	}
};