#pragma once

#include"AD.h"
#include<functional>


class Solver {
public:
	AD::ptr cUxx, cUyy, cUx, cUy, cc, fXL, fXU, fYL, fYU;
	double bXL = 0, bYL = 0, bXU = 1, bYU = 1;
	unsigned int resx = 3;
	unsigned int resy = 3;

private:
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
	static Vec<double> GaussElimination(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.nRows();
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
					M.Populate(k, j, M(i, j));
					M.Populate(i, j, temp);
				}
				temp = u(k);
				u(k) = u(i);
				u(i) = temp;
			}
			for (j = i + 1; j < N; j++)
			{
				temp = M(j, i) / M(i, i);
				for (k = 0; k < N; k++)
					M.Populate(j, k, M(j, k) - temp * M(i, k));
				u(j) = u(j) - temp * u(i);
			}
		}
		x(N - 1) = u(N - 1) / M(N - 1, N - 1);
		for (i = N - 2; i >= 0; i--)
		{
			x(i) = u(i);
			for (j = i + 1; j < N; j++)
				x(i) = x(i) - M(i, j) * x(j);
			x(i) = x(i) / M(i, i);
		}
		return x;
	}

	Vec<double> LU_Decomposition(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.nRows();
		Vec<double> x(N), z(N);
		Matrix<double> l(N, N), u1(N, N);
		int i, j, k;

		for (i = 0; i < N; i++)
			l.Populate(i, 0, M(i, 0));
		for (j = 1; j < N; j++)
			u1.Populate(0, j, M(0, j) / l(0, 0));
		for (i = 0; i < N; i++)
			u1.Populate(i, i, 1);
		for (i = 1; i < N; i++)
			for (j = 1; j < N; j++)
				if (i >= j)
				{
					l.Populate(i, j, M(i, j));
					for (k = 0; k <= j - 1; k++) {

						l.Populate(i, j, l(i, j) - l(i, k) * u1(k, j));
					}
				}
				else
				{
					u1.Populate(i, j, M(i, j));
					for (k = 0; k <= i - 1; k++)
						u1.Populate(i, j, u1(i, j) - l(i, k) * u1(k, j));
					u1.Populate(i, j, u1(i, j) / l(i, i));
				}

		z[0] = u[0] / l(0, 0);
		for (i = 1; i < N; i++)
		{
			z[i] = u[i];
			for (j = 0; j <= i - 1; j++)
				z[i] -= l(i, j) * z[j];
			z[i] /= l(i, i);
		}
		x[N - 1] = z[N - 1];
		for (i = N - 2; i >= 0; i--)
		{
			x[i] = z[i];
			for (j = i + 1; j < N; j++)
				x[i] -= u1(i, j) * x[j];
		}
		return x;
	}

	Vec<double> TriDiagonal(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.size();
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
			D[i] = u[i];
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
		x[N - 1] = D_star[N - 1];
		for (i = N - 2; i >= 0; i--)
		{
			x[i] = D_star[i] - C_star[i] * x[i + 1];
		}
		return x;
	}

	Vec<double> Gauss_Jacobi(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.size();
		Vec<double> x(N), xn(N);
		int i, j, flag;
		double sum, eps = 0.0001;

		for (i = 0; i < N; i++)
			x[i] = 0;
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
				sum = u[i];
				for (j = 0; j < N; j++)
					if (j != i)
						sum -= M(i, j) * x[j];
				xn[i] = sum / M(i, i);
			}
			flag = 0;
			for (i = 0; i < N; i++)
				if (fabs(x[i] - xn[i]) > eps)
					flag = 1;
			if (flag == 1)
				for (i = 0; i < N; i++)
					x[i] = xn[i];
		} while (flag == 1);
		return xn;
	}

	Vec<double> Gauss_Seidal(Matrix<double>& M, Vec<double>& u)
	{
		int N = u.size();
		int test, i, j;
		double sum, eps = 0.0001;

		Vec<double> x(N), xn(N);
		for (i = 0; i < N; i++)
		{
			x[i] = 0;
		}

		for (i = 0; i < N; i++)
			x[i] = 0;


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
				sum = u[i];
				for (j = 0; j < N; j++)
				{
					if (j < i)
						sum -= M(i, j) * xn[j];
					else if (j > i)
						sum -= M(i, j) * x[j];
				}
				xn[i] = sum / M(i, i);
			}
			test = 0;
			for (i = 0; i < N; i++)
				if (fabs(x[i] - xn[i]) > eps)
					test = 1;
			if (test == 1)
				for (i = 0; i < N; i++)
					x[i] = xn[i];
		} while (test == 1);
		return xn;
	}

	Vec<double> SOR(Matrix<double>& M, Vec<double>& u) {
		Vec<double> x(u.size()), xn(u.size());
		unsigned int i, j, flag;
		double sum, eps = 0.001, w = 1;

		for (i = 0; i < u.size(); i++)
			x[i] = 0;
		do
		{
			for (i = 0; i < u.size(); i++)
			{
				sum = u[i] * w + M(i, i) * x[i];
				for (j = 0; j < u.size(); j++)
				{
					if (j < i)
						sum -= M(i, j) * xn[j] * w;
					else if (j >= i)
						sum -= M(i, j) * x[j] * w;
					xn[i] = sum / M(i, i);
				}
			}
			flag = 0;
			for (i = 0; i < u.size(); i++)
				if (fabs(x[i] - xn[i]) > eps)
					flag = 1;
			if (flag == 1)
				for (i = 0; i < u.size(); i++)
					x[i] = xn[i];
		} while (flag == 1);
		return xn;
	}

	Vec<double> SolveNewton(function<Vec<double>(Matrix<double>&, Vec<double>&)> linearSolver) {
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

	Vec<double> SolveBroyden(function<Vec<double>(Matrix<double>&, Vec<double>&)> linearSolver) {
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

			//REPORT
			tmpB = tmpB + (tmpF * tmpu.Transpose()) / mag;
		}

		Vec<double> ret(N);
		for (unsigned int i = 0; i < N; i++) {
			ret(i) = U.find("u" + to_string(i))->second;
		}
		return ret;
	}
};