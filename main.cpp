#include<iostream>
#include"AD.h"

int main()
{
	AD parsed = AD::parse("sjavkdsbv*(2^2+y1)");
	cout << parsed << "\n";
	cout << (parsed.evaluate({{"y1",2},{"sjavkdsbv",0}}));
	return 0;
}