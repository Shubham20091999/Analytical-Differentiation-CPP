#pragma once

#include"AD.h"

class Solver {
public:
	AD::ptr cUxx, cUyy, cUx, cUy, cc, fXL, fXU, fYL, fYU;
	double bXL = 0, bYL = 0, bXU = 1, bYU = 1;
	unsigned int resx = 4;
	unsigned int resy = 4;

private:
	double dx;
	double dy;

	unsigned int getPos(unsigned int i, unsigned int j) {
		return j * (resx - 2) + i;
	}

public:
	Solver(double a = 0, double bXU = 1, double bYL = 0, double bYU = 1) :
		bXL(a), bXU(bXU), bYL(bYL), bYU(bYU) {
		cUxx = make_shared<AD>(AD::parse("0"));
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

	Matrix<AD> getF() {
		Matrix<AD> ret((resx - 2) * (resy - 2), 1);

		for (unsigned int i = 1; i < resx - 1; i++) {
			for (unsigned int j = 1; j < resy - 1; j++) {
				unsigned int pos = getPos(i - 1, j - 1);
				ret(pos, 0) = *(((cUxx * getUxx(i, j) + cUyy * getUyy(i, j) + cUx * getUx(i, j) + cUy * getUy(i, j) + cc))->copy());
				ret(pos, 0).putVal({ {"x",bXL + i * dx},{"y",bYL + j * dy} });
				ret(pos, 0).replaceUnknown("u", "u" + to_string(pos));
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
	static Matrix<double> GaussElimination(Matrix<double>& M, Matrix<double>& u)
	{
		int N = u.nRows();
		Matrix<double> x(N, 1);
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
					M.Populate(k, j, M(i, j));
					M.Populate(i, j, temp);
				}
				temp = u(k, 0);
				u(k, 0) = u(i, 0);
				u(i, 0) = temp;
			}
			for (j = i + 1; j < N; j++)
			{
				temp = M(j, i) / M(i, i);
				for (k = 0; k < N; k++)
					M.Populate(j, k, M(j, k) - temp * M(i, k));
				u(j, 0) = u(j, 0) - temp * u(i, 0);
			}
		}
		x(N - 1, 0) = u(N - 1, 0) / M(N - 1, N - 1);
		for (i = N - 2; i >= 0; i--)
		{
			x(i, 0) = u(i, 0);
			for (j = i + 1; j < N; j++)
				x(i, 0) = x(i, 0) - M(i, j) * x(j);
			x(i) = x(i) / M(i, i);
		}
		return x;
	}

	Matrix<double> solve() {
		auto F = getF();
		auto J = getJ(F);

		unsigned int N = F.nRows();

		map<string, double> U;
		Matrix<double> tmpF(N, 1);
		Matrix<double> tmpJ(N, N);
		for (unsigned int i = 0; i < N; i++) {
			U.insert({ "u" + to_string(i),(double)rand() / RAND_MAX * 5 });
		}

		for (int _ = 0; _ < 100; _++) {
			for (unsigned int i = 0; i < N; i++) {
				tmpF(i, 0) = -F(i, 0).evaluate(U);
				for (unsigned int j = 0; j < N; j++) {
					tmpJ(i, j) = J(i, j).evaluate(U);
				}
			}

			Matrix<double> tmpu = GaussElimination(tmpJ, tmpF);

			for (unsigned int i = 0; i < N; i++) {
				std::map<string, double>::iterator it = U.find("u" + to_string(i));
				if (it != U.end())
					it->second += tmpu(i);
			}
		}

		Matrix<double> ret(N, 1);
		for (unsigned int i = 0; i < N; i++) {
			ret(i) = U.find("u" + to_string(i))->second;
		}
		return ret;
	}
};