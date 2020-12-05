#pragma once
#include"AD.h"

class Solver {
public:
	AD::ptr A, B, C, D, E;
	double a = 0, c = 0, b = 1, d = 1;
	unsigned int resx = 8;
	unsigned int resy = 8;

	double dx;
	double dy;

	unsigned int getPos(unsigned int i, unsigned int j) {
		return j * (resx - 2) + i;
	}

	Solver(int lol) {
		A = make_shared<AD>(AD::parse("1"));
		B = make_shared<AD>(AD::parse("1"));
		C = make_shared<AD>(AD::parse("0"));
		D = make_shared<AD>(AD::parse("0"));
		E = make_shared<AD>(AD::parse("-2"));
		dx = 1 / (double)(resx - 1);
		dy = 1 / (double)(resy - 1);
	}

	AD::ptr getVal(unsigned int i, unsigned int j) {
		if (i == 0)
		{
			return make_shared<AD>(AD(0));
		}
		if (i == resx - 1) {
			return make_shared<AD>(AD(0));
		}
		if (j == 0) {
			return make_shared<AD>(AD(0));
		}
		if (j == resy - 1) {
			return make_shared<AD>(AD(0));
		}
		return make_shared<AD>(AD("u" + to_string(getPos(i - 1, j - 1))));
	}

	AD::ptr getUxx(int i, int j) {
		return (getVal(i + 1, j) - AD::getNum(2) * getVal(i, j) + getVal(i - 1, j)) * AD::getNum(pow(resx - 1, 2));
	}

	AD::ptr getUyy(int i, int j) {
		return (getVal(i, j + 1) - AD::getNum(2) * getVal(i, j) + getVal(i, j - 1)) * AD::getNum(pow(resy - 1, 2));
	}

	AD::ptr getUx(int i, int j) {
		return (getVal(i + 1, j) - getVal(i - 1, j)) * (AD::getNum((resx - 1) / 2));
	}

	AD::ptr getUy(int i, int j) {
		return (getVal(i, j + 1) - getVal(i, j - 1)) * (AD::getNum((resy - 1) / 2));
	}

	AD::ptr getU(int i, int j) {
		return getVal(i, j);
	}

	SparseMatrix<AD> getF() {
		SparseMatrix<AD> ret((resx - 2) * (resy - 2), 1);

		for (unsigned int i = 1; i < resx - 1; i++) {
			for (unsigned int j = 1; j < resy - 1; j++) {
				unsigned int pos = getPos(i - 1, j - 1);
				ret(pos, 0, *(A * getUxx(i, j) + B * getUyy(i, j) + C * getUx(i, j) + D * getUy(i, j) + E));
			}
		}

		return ret;
	}

	SparseMatrix<AD> getJ() {
		SparseMatrix<AD> F = getF();
		SparseMatrix<AD> ret((resx - 2) * (resy - 2), (resx - 2) * (resy - 2));

		for (int i = 0; i < (resx - 2) * (resy - 2); i++) {
			for (int j = 0; j < (resx - 2) * (resy - 2); j++) {
				ret(i, j, *F(i, 0).derivative("u" + to_string(j)));
			}
		}
		return ret;
	}

};