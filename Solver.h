#pragma once

#include"AD.h"
#include<functional>


class Solver {
public:
	AD::ptr cUxx, cUyy, cUx, cUy, cc, fXL, fXU, fYL, fYU;
	double bXL = 0, bYL = 0, bXU = 1, bYU = 1;

private:
	unsigned int resx = 4;
	unsigned int resy = 4;
	double dx;
	double dy;

	unsigned int getPos(unsigned int i, unsigned int j) {
		return j * (resx - 2) + i;
	}

public:
	Solver(double a = 0, double bXU = 1, double bYL = 0, double bYU = 1) :
		bXL(a), bXU(bXU), bYL(bYL), bYU(bYU) {
		cUxx = make_shared<AD>(AD(0));
		cUyy = cUxx;
		cUx = cUxx;
		cUy = cUxx;
		cc = cUxx;
		fXL = cUxx;
		fXU = cUxx;
		fYL = cUxx;
		fYU = cUxx;

		dx = (bXU - a) / (double)(resx - 1);
		dy = (bYU - bYL) / (double)(resy - 1);
	}

	void setcUxx(AD&& a) {
		cUxx = make_shared<AD>(a);
	}

	void setcUyy(AD&& a) {
		cUyy = make_shared<AD>(a);
	}
	void setcUx(AD&& a) {
		cUx = make_shared<AD>(a);
	}
	void setcUyx(AD&& a) {
		cUy = make_shared<AD>(a);
	}

	void setcc(AD&& a) {
		cc = make_shared<AD>(a);
	}

private:
	AD::ptr getVal(unsigned int i, unsigned int j) {
		if (i == 0)
			return make_shared<AD>(AD(fXL->evaluate({ {"x",bXL + i * dx},{"y",bYL + j * dy} })));
		if (i == resx - 1)
			return make_shared<AD>(AD(fXU->evaluate({ {"x",bXL + i * dx},{"y",bYL + j * dy} })));
		if (j == 0)
			return make_shared<AD>(AD(fYL->evaluate({ {"x",bXL + i * dx},{"y",bYL + j * dy} })));
		if (j == resy - 1)
			return make_shared<AD>(AD(fYU->evaluate({ {"x",bXL + i * dx},{"y",bYL + j * dy} })));
		return make_shared<AD>(AD("u" + to_string(getPos(i - 1, j - 1))));
	}

	AD::ptr getUxx(int i, int j) {
		return (getVal(i + 1, j) - AD::getNum(2) * getVal(i, j) + getVal(i - 1, j)) * AD::getNum(pow(1 / dx, 2));
	}

	AD::ptr getUyy(int i, int j) {
		return (getVal(i, j + 1) - AD::getNum(2) * getVal(i, j) + getVal(i, j - 1)) * AD::getNum(pow(1 / dy, 2));
	}

	AD::ptr getUx(int i, int j) {
		return (getVal(i + 1, j) - getVal(i - 1, j)) * (AD::getNum(1 / dx / 2));
	}

	AD::ptr getUy(int i, int j) {
		return (getVal(i, j + 1) - getVal(i, j - 1)) * (AD::getNum(1 / dy / 2));
	}

	AD::ptr getU(int i, int j) {
		return getVal(i, j);
	}

	Vec<AD> getF() {
		Vec<AD> ret((resx - 2) * (resy - 2));

		for (unsigned int i = 1; i < resx - 1; i++) {
			for (unsigned int j = 1; j < resy - 1; j++) {
				unsigned int pos = getPos(i - 1, j - 1);
				ret(pos) = *(((cUxx * getUxx(i, j) + cUyy * getUyy(i, j) + cUx * getUx(i, j) + cUy * getUy(i, j) + cc))->copy());
				ret(pos).putVal({ {"x",bXL + i * dx},{"y",bYL + j * dy} });
				ret(pos).replaceUnknown("u", "u" + to_string(pos));
			}
		}

		return ret;
	}

	Matrix<AD> getJ(Matrix<AD> F) {
		Matrix<AD> ret((resx - 2) * (resy - 2), (resx - 2) * (resy - 2));
		unsigned int tmp = (resx - 2) * (resy - 2);
		for (unsigned int i = 0; i < tmp; i++) {
			for (unsigned int j = 0; j < tmp; j++) {
				ret(i, j) = *F(i, 0).derivative("u" + to_string(j));
			}
		}
		return ret;
	}

public:
	Vec<double> SolveNewton(const function<Vec<double>(Matrix<double>&, Vec<double>&)>& linearSolver) {
		auto F = getF();
		auto J = getJ(F);

		unsigned int N = F.nRows();

		map<string, double> U;
		Vec<double> tmpF(N);
		Matrix<double> tmpJ(N, N);
		for (unsigned int i = 0; i < N; i++) {
			U.insert({ "u" + to_string(i),(double)rand() / RAND_MAX * 5 });
		}

		for (int _ = 0; _ < 100; _++) {
			for (unsigned int i = 0; i < N; i++) {
				tmpF(i) = -F(i).evaluate(U);
				for (unsigned int j = 0; j < N; j++) {
					tmpJ(i, j) = J(i, j).evaluate(U);
				}
			}


			Vec<double> tmpu = linearSolver(tmpJ, tmpF);

			double mag = tmpu.magnitude();
			if (mag < 0.00001)
				break;

			for (unsigned int i = 0; i < N; i++) {
				std::map<string, double>::iterator it = U.find("u" + to_string(i));
				if (it != U.end())
					it->second += tmpu(i);
			}
		}

		Vec<double> ret(N);
		for (unsigned int i = 0; i < N; i++) {
			ret(i) = U.find("u" + to_string(i))->second;
		}
		return ret;
	}

	Vec<double> SolveBroyden(const function<Vec<double>(Matrix<double>&, Vec<double>&)>& linearSolver) {
		auto F = getF();
		auto B = getJ(F);

		unsigned int N = F.nRows();

		//
		map<string, double> U;
		Vec<double> tmpF(N);
		Matrix<double> tmpB(N, N);

		//Initializing U
		for (unsigned int i = 0; i < N; i++) {
			U.insert({ "u" + to_string(i),(double)rand() / RAND_MAX * 5 });
		}

		//Initializing B
		for (unsigned int i = 0; i < N; i++) {
			tmpF(i) = -F(i).evaluate(U);
			for (unsigned int j = 0; j < N; j++) {
				tmpB(i, j) = B(i, j).evaluate(U);
			}
		}

		//Loops
		for (int _ = 0; _ < 100; _++) {

			Vec<double> tmpu = linearSolver(tmpB, tmpF);

			double mag = tmpu.magnitude();
			if (mag < 0.00001)
				break;

			for (unsigned int i = 0; i < N; i++) {
				std::map<string, double>::iterator it = U.find("u" + to_string(i));
				if (it != U.end())
					it->second += tmpu(i);
			}

			for (unsigned int i = 0; i < N; i++) {
				tmpF(i) = -F(i).evaluate(U);
			}

			tmpB = tmpB + (tmpF * tmpu.Transpose()) / mag;
		}

		Vec<double> ret(N);
		for (unsigned int i = 0; i < N; i++) {
			ret(i) = U.find("u" + to_string(i))->second;
		}
		return ret;
	}
};