#include<iostream>
#include"AD.h"
#include"Solver.h"
#include<fstream>

// struct Data {
// 	AD::ptr* ptr;
// 	string s = "0";
// 	Data(AD::ptr* a) {
// 		ptr = a;
// 	}
// };

int main() {

	std::string function="1/-sec(x)";
	AD::ptr parsed=AD::parse(function);
	std::cout<<*parsed;
	AD::ptr der=parsed->derivative("x");
	std::cout<<der->evaluate({{"x",20}});

	//std::cout<<der->evaluate({{"x",2}});

	//for (auto x : AD::infixToPostfix(function))
	//{
	//	cout<<x<<" ";
	//}

	// //Solver
	// //llx->lower limit x
	// //ulx->upper limit x
	// //lly->lower limit y
	// //ulx->upper limit y
	// double llx;
	// double ulx;
	// double lly;
	// double uly;
	// cin >> llx >> ulx >> lly >> uly;
	// Solver s(2, 2, llx, ulx, lly, uly);
	// //cUxx-> coefficient of d2z/dx2
	// //cUyy-> coefficient of d2y/dy2
	// //cUx-> coefficient of dz/dx
	// //cUy-> coefficient of dz/dy
	// //cc-> value of the constant term
	// //in cUxx*d2z/dx2+cUyy*d2y/dy2+cUx*dz/dx+cUy*dz/dy+cc=0
	// //fXL-> Lower limit function on X
	// //fXU-> Upper limit function on X
	// //fYL-> Lower limit function on Y
	// //fYU-> Upper limit function on X

	// vector<pair<string, Data>> sinps = {
	// {"cUxx",&s.cUxx},
	// {"cUyy",&s.cUyy},
	// {"cUx",&s.cUx},
	// {"cUy",&s.cUy},
	// {"cc",&s.cc},
	// {"fXL",&s.fXL},
	// {"fXU",&s.fXU},
	// {"fYL",&s.fYL},
	// {"fYU",&s.fYU},
	// };

	// for (auto& a : sinps) {
	// 	cout << a.first << ": ";
	// 	cin >> a.second.s;
	// 	cout << "\n";
	// 	*a.second.ptr = AD::parse(a.second.s);
	// }

	// auto t = s.SolveNewton(linearSolvers::GaussElimination);

	// //ofstream out;
	// //out.open("out.txt");
	// for (unsigned int i = 0; i < t.dim(); i++) {
	// 	cout << t(i) << "\n";
	// }


	return 0;

}


