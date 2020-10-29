#include<iostream>
#include"AD.h"

int main()
{
	AD parsed = AD::parse("x^x");
	cout << parsed << "\n";
	cout << (parsed.derivative("x"))->evaluate({{"x",2}});
	return 0;
}