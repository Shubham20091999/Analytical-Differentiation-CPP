#include<iostream>
#include"AD.h"

int main()
{
	AD parsed = AD::parse("sin(-x)");
	cout << parsed;
	cout<<parsed.evaluate({{"x",1.570796327}});
	
	return 0;
}