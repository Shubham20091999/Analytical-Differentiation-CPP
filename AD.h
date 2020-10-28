#pragma once
#include <variant>
#include <string>
#include<map>
#include<list>
#include<stack>

using namespace std;


class  AD
{
	enum class type {
		function,
		operation,
		known,
		unknown,
	};

	const map<string, double (*)(double)> m = { {"sin",sin} };

	std::variant<std::string, double> value;
	AD* left = nullptr;
	AD* right = nullptr;
	type typ;

public:
	AD(double x)
	{
		value = x;
		typ = type::known;
	}

	AD(string x)
	{
		value = x;
		typ = type::unknown;
	}

	AD(string x, AD* r)
	{
		value = x;
		right = r;
		typ = type::function;
	}

	AD(string x, AD* l, AD* r)
	{
		value = x;
		left = l;
		right = r;
		typ = type::operation;
	}

	friend ostream& operator<<(ostream& out, const  AD& exp)
	{
		if (exp.typ == type::known)
			out << std::get<double>(exp.value);
		else if (exp.typ == type::unknown)
			out << std::get<string>(exp.value);
		else if (exp.typ == type::function)
			out << std::get<string>(exp.value) << '(' << *(exp.right) << ')';
		else if (exp.typ == type::operation)
			out << '(' << *exp.left << std::get<string>(exp.value) << *exp.right << ')';
		return out;
	}

	static AD parse(const string& exp)
	{

	}

private:
	static int precedence(char c)
	{
		if (c == '^')
			return 3;
		else if (c == '*' || c == '/')
			return 2;
		else if (c == '+' || c == '-')
			return 1;
		else if (c == '(' || c == ')')
			return 0;
		return -1;
	}

	static list<string> infixToPostfix(const string& exp)
	{
		list<variant<AD, string>> a;
		stack<string> stack;

		for (int i = 0; i < exp.size(); i++)
		{
			if (precedence(exp[i]) >= 0)
			{
				
			}
		}
	}

private:
	static  list<string> parse_toList(const string& exp)
	{

		for (int i = 0; i < )
	}
};