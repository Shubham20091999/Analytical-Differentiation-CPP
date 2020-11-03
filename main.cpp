#include<iostream>
#include"AD.h"

int main()
{
	string s = "cos(sin(tan(y1))++++++-+++++20)";
	AD parsed = AD::parse(s);
	for (auto a : AD::infixToPostfix(s))
	{
		cout << a << " ";
	}
	cout << '\n';
	cout << parsed;
	cout << (parsed.evaluate({ {"y1",2} }));
	return 0;
}