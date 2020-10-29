#pragma once
#include <variant>
#include <string>
#include<map>
#include<list>
#include<stack>
//#include<sstream>
#include<functional>
#include"Functions.h"

using namespace std;

class  AD
{
public:
	enum class type {
		function,
		operation,
		known,
		unknown,
	};

	static map < string, double (*)(double) > functions;
	static map<string, double (*)(double, double) > operations;

	static map<string, std::shared_ptr<AD>(*)(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)> doperations;

	std::variant<std::string, double> value;
	shared_ptr<AD> left = nullptr;
	shared_ptr<AD> right = nullptr;
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

	AD(string x, AD* r) :
		value(x), right(r), typ(type::function)
	{
	}

	AD(string x, AD* l, AD* r) :
		value(x), left(l), right(r), typ(type::operation)
	{
	}

	AD(string x, std::shared_ptr<AD> r) :
		value(x), right(r), typ(type::function)
	{

	}

	AD(string x, std::shared_ptr<AD> l, std::shared_ptr<AD> r) :
		value(x), left(l), right(r), typ(type::operation)
	{
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
		{
			//std::stringstream l;
			//l << *exp.left;
			//std::ostringstream r;
			//r << *exp.right;
			//out << '(' << l.str() << std::get<string>(exp.value) << r.str() << ')';
			out << '(' << *exp.left << std::get<string>(exp.value) << *exp.right << ')';
		}
		return out;
	}


	double evaluate(const map<string, double>& vars = {})
	{

		if (typ == type::known)
		{
			return std::get<double>(value);
		}
		else if (typ == type::unknown)
		{
			std::map<string, double>::const_iterator it = vars.find(std::get<string>(value));
			if (it == vars.end())
			{
				throw;
			}
			return it->second;
		}
		else if (typ == type::function)
		{
			map<string, double (*)(double)>::const_iterator it = AD::functions.find(std::get<string>(value));
			if (it == AD::functions.end())
			{
				throw;
			}
			return (it->second)(right->evaluate(vars));
		}
		else/* (typ == type::operation)*/
		{
			map<string, double (*)(double, double)>::const_iterator it = AD::operations.find(std::get<string>(value));
			if (it == AD::operations.end()) {
				throw;
			}
			return (it->second)(left->evaluate(vars), right->evaluate(vars));
		}
	}

	static AD parse(const string& exp)
	{
		list<string> lst = infixToPostfix(exp);

		stack<AD*> stack;

		for (auto e = lst.begin(); e != lst.end(); e++)
		{
			int pre = precedence(*e);
			if (pre > 1 and pre < 5)
			{
				AD* b = stack.top();
				stack.pop();
				AD* a = stack.top();
				stack.pop();
				stack.push(new AD(*e, a, b));
			}
			else if (pre == 5)
			{
				AD* x = stack.top();
				stack.pop();
				stack.push(new AD(*e, x));
			}
			else {
				if (AD::isNumber(*e))
					stack.push(new AD(std::stod(*e)));
				else
					stack.push(new AD(*e));
			}
		}
		return *stack.top();
	}

	static int precedence(char c)
	{
		if (c == '^')
			return 4;
		else if (c == '*' || c == '/')
			return 3;
		else if (c == '+' || c == '-')
			return 2;
		else if (c == '(')
			return 1;
		else if (c == ')')
			return 0;
		return -1;
	}

	static int precedence(const string& c)
	{
		if (AD::functions.find(c) != AD::functions.end())
			return 5;
		if (c == "^")
			return 4;
		if (c == "*" || c == "/")
			return 3;
		if (c == "+" || c == "-")
			return 2;
		if (c == "(")
			return 1;
		if (c == ")")
			return 0;
		return -1;
	}

	static list<string> infixToPostfix(const string& exp)
	{
		std::list<string> lst = AD::parse_toList(exp);
		list<string> a;
		stack<string> stack;

		for (auto e = lst.begin(); e != lst.end(); e++)
		{
			if (*e == "(")
			{
				stack.push(*e);
			}
			else if (*e == ")")
			{
				while (!stack.empty() && stack.top() != "(")
				{
					a.push_back(stack.top());
					stack.pop();
				}
				if (stack.top() == "(")
					stack.pop();
			}
			else if (precedence(*e) >= 2) {
				while (!stack.empty() && precedence(*e) <= precedence(stack.top()))
				{
					a.push_back(stack.top());
					stack.pop();
				}
				stack.push(*e);
			}
			else {
				a.push_back(*e);
			}
		}
		while (!stack.empty())
		{
			a.push_back(stack.top());
			stack.pop();
		}
		return a;
	}

	static bool isNumber(const string& s)
	{
		size_t i = 0;
		if (s[0] == '-')
			i++;
		for (; i < s.size(); i++)
			if (isdigit(s[i]) == false and s[i] != '.')
				return false;
		return true;
	}

	static  std::list<string> parse_toList(const string& exp)
	{
		string container = "";
		std::list<string> ans;
		for (auto e = exp.begin(); e != exp.end(); e++)
		{
			if (*e == ' ')
				continue;
			if (AD::precedence(*e) >= 0)
			{
				if (container.size() == 0 && *e == '-')
				{
					ans.push_back("-1");
					ans.push_back("*");
				}
				else if (container.size() == 0 && *e == '+')
					continue;
				else
				{
					if (container.size() != 0)
						ans.push_back(container);
					ans.push_back(string(1, *e));
				}
				container.clear();
			}
			else {
				container.push_back(*e);
			}
		}
		if (container.size() != 0)
			ans.push_back(container);

		return ans;
	}

	static double add(double a, double b) {
		return a + b;
	}

	static double pwr(double a, double b)
	{
		return pow(a, b);
	}

	static double sub(double a, double b)
	{
		return a - b;
	}

	static double mul(double a, double b)
	{
		return a * b;
	}

	static double div(double a, double b)
	{
		return a / b;
	}


	bool operator ==(double b)
	{
		if (typ == type::known)
		{
			if (std::get<double>(value) == b)
			{
				return true;
			}
		}
		return false;
	}

	int valType_num()
	{
		if (typ == type::known)
			return 0;
		if (typ == type::unknown)
			return 1;
		if (typ == type::function)
			return 2;
		if (typ == type::operation)
			return 3;
	}

	friend std::shared_ptr<AD> operator *(std::shared_ptr<AD> a, std::shared_ptr<AD> b)
	{
		if (a->valType_num() == 0 and b->valType_num() == 0)
			return std::shared_ptr<AD>(new AD(std::get<double>(a->value) * std::get<double>(b->value)));
		if (*a == 0 or *b == 0)
			return std::shared_ptr<AD>(new AD(0));
		if (*a == 1)
			return b;
		if (*b == 1)
			return a;
		return std::shared_ptr<AD>(new AD("*", a, b));
	}

	friend std::shared_ptr<AD> operator +(std::shared_ptr<AD> a, std::shared_ptr<AD> b)
	{
		if (a->valType_num() == 0 and b->valType_num() == 0)
			return std::shared_ptr<AD>(new AD(std::get<double>(a->value) + std::get<double>(b->value)));

		if (*a == 0)
			return b;
		if (*b == 0)
			return a;
		return std::shared_ptr<AD>(new AD("+", a, b));
	}

	friend std::shared_ptr<AD> operator -(std::shared_ptr<AD> a, std::shared_ptr<AD> b)
	{
		if (a->valType_num() == 0 and b->valType_num() == 0)
			return std::shared_ptr<AD>(new AD(std::get<double>(a->value) - std::get<double>(b->value)));
		if (*a == 0)
			return std::shared_ptr<AD>(new AD(-1)) * b;
		if (*b == 0)
			return a;
		return std::shared_ptr<AD>(new AD("-", a, b));
	}

	friend std::shared_ptr<AD> operator /(std::shared_ptr<AD> a, std::shared_ptr<AD> b)
	{
		if (a->valType_num() == 0 and b->valType_num() == 0)
			return std::shared_ptr<AD>(new AD(std::get<double>(a->value) / std::get<double>(b->value)));
		if (*a == 0)
			return  std::shared_ptr<AD>(new AD(0));
		return std::shared_ptr<AD>(new AD("/", a, b));
	}

	friend std::shared_ptr<AD> operator ^(std::shared_ptr<AD> a, std::shared_ptr<AD> b)
	{
		if (a->valType_num() == 0 and b->valType_num() == 0)
			return std::shared_ptr<AD>(new AD(pow(std::get<double>(a->value), std::get<double>(b->value))));
		if (*a == 0)
			return std::shared_ptr<AD>(new AD(0));
		if (*b == 0)
			return std::shared_ptr<AD>(new AD(1));
		if (*a == 1)
			return std::shared_ptr<AD>(new AD(1));
		if (*b == 1)
			return a;
		return std::shared_ptr<AD>(new AD("^", a, b));
	}


	static std::shared_ptr<AD> dmul(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)
	{
		return a * b->derivative(x) + b * a->derivative(x);
	}

	static std::shared_ptr<AD> dadd(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)
	{
		return a->derivative(x) + b->derivative(x);
	}

	static std::shared_ptr<AD> dsub(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)
	{
		return a->derivative(x) - b->derivative(x);
	}

	static std::shared_ptr<AD> ddiv(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)
	{
		return (a->derivative(x) * b - b->derivative(x) * a) / (b ^ std::shared_ptr<AD>(new AD(2)));
	}

	static std::shared_ptr<AD> dpow(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)
	{
		return (a ^ b) * (a->derivative(x) * std::shared_ptr<AD>(new AD("ln", b)) + a * b->derivative(x) / b);
	}


	std::shared_ptr<AD> dtan(std::shared_ptr<AD> a)
	{
		return std::shared_ptr<AD>(new AD("^", new AD("cos", a), new AD(-2)));
	}


	std::shared_ptr<AD> derivative(const string& x)
	{
		if (typ == type::known)
		{
			return std::shared_ptr<AD>(new AD(0));
		}
		if (typ == type::unknown)
		{
			if (std::get<string>(value) == x)
				return std::shared_ptr<AD>(new AD(1));
			else
				return std::shared_ptr<AD>(new AD(0));
		}
		if (typ == type::function)
		{
			if (std::get<string>(value) == "tan")
			{
				return AD::dtan(right) * right->derivative(x);
			}
		}
		if (typ == type::operation)
		{
			std::map<string, std::shared_ptr<AD>(*)(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)>::const_iterator it = doperations.find(std::get<string>(value));
			if (it == AD::doperations.end())
				throw;
			return (it->second)(left, right, x);
		}
	}
};

map<string, double (*)(double) > AD::functions = { {"sin",std::sin} ,{"cos",std::cos},{"tan",std::tan},{"log10",std::log10},{"exp",std::exp},{"ln",std::log} };

map<string, std::shared_ptr<AD>(*)(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)> AD::doperations =
{ {"+",AD::dadd},{"*",AD::dmul},{"-",AD::dsub},{"/",AD::ddiv},{"^",AD::dpow} };


map<string, double (*)(double, double) > AD::operations = { {"^", AD::pwr} , { "+", AD::add } ,{"-",AD::sub},{"*",AD::mul},{"/",AD::div} };