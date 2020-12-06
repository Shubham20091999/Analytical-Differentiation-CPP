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


	T& Glimpse(unsigned int i, unsigned int j) {
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

	static Matrix<T> rand(unsigned int m, unsigned int n) {
		Matrix<double> ret(m, n);
		for (unsigned int i = 0; i < m; i++) {
			for (unsigned int j = 0; j < n; j++) {
				ret(i, j) = T(rand()) / RAND_MAX * 5;
			}
		}
		return ret;
	}
};



template <class T>
class SparseMatrix {
private:
	int m, n;
	vector<int> row, col;
	vector<T> value;

	void construct(int, int);

public:
	//constructor overloading
	SparseMatrix() {};
	SparseMatrix(int n); // square matrix of dimension nxn
	SparseMatrix(int rows, int cols); // general matrix

	//operator overloading
	vector<T> operator *(vector<T>);
	void operator = (Matrix<T>);

	SparseMatrix<T> operator +(Matrix<T>);
	SparseMatrix<T> operator -(Matrix<T>);

	T operator () (int, int);
	SparseMatrix<T>& operator () (int, int, T);

	//member functions
	//SparseMatrix<T>& set(T, int, int);
	//SparseMatrix<T>& input();
	T get(int, int);
	int getRowCount(void);
	int getColCount(void);
	void disp();
	//void dispCol();
	//void dispPos();
	//void dispValue();

	vector<int> getRow(void);
	vector<int> getCol(void);
	vector<T> getValue(void);

	template<class U>
	friend void disp(vector<U>);

	//template<class U>
	//friend vector<U> operator /(vector<U>, U);

	template<class U>
	friend SparseMatrix<U> operator *(vector<U>, vector<U>);


	template <class U>
	friend U norm(vector<U>);

	unsigned int nCols() {
		return n;
	}

	unsigned int nRows() {
		return m;
	}
};

//initializing using constructor
template<class T>
SparseMatrix<T>::SparseMatrix(int n) { this->construct(n, n); }

template<class T>
SparseMatrix<T>::SparseMatrix(int rows, int cols) { this->construct(rows, cols); }

template<class T>
void SparseMatrix<T>::construct(int rows, int cols) {
	this->m = rows;
	this->n = cols;
}

template<class T>
T SparseMatrix<T>::operator () (int i, int j) { return this->get(i, j); }

template<class T>
SparseMatrix<T>& SparseMatrix<T>::operator () (int i, int j, T k) {
	if (!(k == 0)) {
		if (this->get(i, j) == 0) {
			this->row.push_back(i);
			this->col.push_back(j);
			this->value.push_back(k);
		}
		else {
			for (unsigned int p = 0; p < this->row.size(); p++) {
				if ((this->row[p] == i) && (this->col[p] == j)) {
					this->value[p] = k;
					break;
				}
			}
		}
	}
	return *this;
}

//retrieve value from SM using index i,j
template<class T>
T SparseMatrix<T>::get(int i, int j) {
	//create an exception handler for out of bound
	if (i >= this->m) {
		cout << i << "Must be less than " << this->m << endl;
		exit(0);
	}
	if (j >= this->n) {
		cout << j << "Must be less than " << this->n << endl;
		exit(0);
	}
	T k = 0;
	for (unsigned int p = 0; p < this->row.size(); p++) {
		if ((this->row[p] == i) && (this->col[p] == j)) {
			k = this->value[p];
			break;
		}
	}
	return k;
}

//get row & column count
template<class T>
int SparseMatrix<T>::getRowCount(void) { return this->m; }

template<class T>
int SparseMatrix<T>::getColCount(void) { return this->n; }

template<class T>
void SparseMatrix<T>::disp() {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			cout << this->get(i, j) << '\t';
		cout << endl;
	}
}

template<class T>
vector<T> SparseMatrix<T>::operator *(vector<T> v) {
	if (v.size() != this->n) {
		cout << "MATH ERROR: Invalid dimension, cannot be multiplied!" << endl;
	}
	vector<T> v1(m, 0);
	for (int i = 0; i < this->m; i++)
		for (int j = 0; j < v.size(); j++)
			v1[i] = v1[i] + this->get(i, j) * v[j];
	return v1;
}

template<class T>
void SparseMatrix<T>::operator = (Matrix<T> u) {
	this->m = u.m;
	this->n = u.n;
	this->row = u.row;
	this->col = u.col;
	this->value = u.value;
}

template<class T>
SparseMatrix<T> SparseMatrix<T>::operator +(Matrix<T> A) {
	if (this->m != A.m || this->m != A.m) {
		cout << "MATH ERROR: Matrix Dimensions don't match, Addition Halted!\n";
		exit(1);
	}

	Matrix<T> sum(A.m, A.n);

	for (int i = 0; i < A.m; i++)
		for (int j = 0; j < A.n; j++)
			sum(i, j, (*this)(i, j) + A(i, j));
	return sum;
}

template<class T>
SparseMatrix<T> SparseMatrix<T>::operator -(Matrix<T> A) {
	if (this->m != A.m || this->m != A.m) {
		cout << "MATH ERROR: Matrix Dimensions don't match, Addition Halted!\n";
		exit(1);
	}

	Matrix<T> diff(A.m, A.n);

	for (int i = 0; i < A.m; i++)
		for (int j = 0; j < A.n; j++)
			diff(i, j, (*this)(i, j) - A(i, j));
	return diff;
}

template<class T>
vector<int> SparseMatrix<T>::getRow(void) { return this->row; }

template<class T>
vector<int> SparseMatrix<T>::getCol(void) { return this->col; }

template<class T>
vector<T> SparseMatrix<T>::getValue(void) { return this->value; }

template<class T>
void disp(vector<T> v) {
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << "   ";
	cout << endl;
}

template<class T>
vector<T> operator /(vector<T> v, T k) {
	if (k == 0) {
		cout << "MATH ERROR: Division By Zero!\n";
		exit(1);
	}

	vector<T> u(v.size());
	for (int i = 0; i < v.size(); i++)
		u[i] = v[i] / k;

	return u;
}

template<class T>
Matrix<T> operator *(vector<T> u, vector<T> v) {
	Matrix<T> P(u.size(), v.size());
	for (int i = 0; i < u.size(); i++)
		for (int j = 0; j < v.size(); j++)
			P(i, j, u[i] * v[j]);
	return P;
}

template<class T>
vector<T> operator -(vector<T> u, vector<T> v) {
	if (u.size() != v.size()) {
		cout << "MATH ERROR: Vector Sizes don't Match!\n";
		exit(1);
	}

	vector<T> w(v.size());
	for (int i = 0; i < v.size(); i++)
		w[i] = u[i] - v[i];

	return w;
}

template <class T>
T norm(vector<T> u) {
	T s = 0.0;
	for (int i = 0; i < u.size(); i++)
		s = s + u[i] * u[i];
	return sqrt(s);
}