#pragma once
#include <iostream>

class VariableNotFound {
public:
	VariableNotFound(std::string s)
	{
		std::cout << "Variable not found: " << s << '\n';
	}
};

class FunctionNotFound {
public:
	FunctionNotFound(std::string s)
	{
		std::cout << "Function not found: " << s << '\n';
	}
};

class FunctionDerivativeNotFound {
public:
	FunctionDerivativeNotFound(std::string s)
	{
		std::cout << "Derivative not found for Function: " << s << '\n';
	}
};

class OperatorNotFound {
public:
	OperatorNotFound(std::string s)
	{
		std::cout << "Operator not found: " << s << '\n';
	}
};

class OperatorDerivativeNotFound {
public:
	OperatorDerivativeNotFound(std::string s)
	{
		std::cout << "Derivative not found for Operator: " << s << '\n';
	}
};

class IndexOutOfRange {
public:
	IndexOutOfRange(unsigned int i, unsigned int m)
	{
		std::cout << "ERROR: Index error( " << i << " should be less than " << m << ")\n";
	}
};

class GenException {
public:
	GenException(const std::string& s)
	{
		std::cout << s << "\n";
	}
};