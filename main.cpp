#include<iostream>
#include"AD.h"

int main()
{
	AD _2(2);
	AD _3(3);
	AD _4(4);
	AD x("x");
	AD sin("sin", &x);
	AD _add("+", &_2, &_3);
	AD _mul("*", &_add, &sin);
	cout<<_mul;
	return 0;
}