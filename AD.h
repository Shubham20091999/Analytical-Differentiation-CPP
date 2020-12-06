#pragma once
#include <variant>
#include <string>
#include<map>
#include<list>
#include<stack>
#include<functional>
#include"Exceptions.h"
#include"Matrix.h"

using namespace std;

namespace extrafncs {
	static double add(double a, double b) {
		return a + b;
	}

	static double pow(double a, double b)
	{
		return std::pow(a, b);
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

	static double sec(double x)
	{
		return 1 / std::cos(x);
	}

	static double cosec(double x)
	{
		return 1 / std::sin(x);
	}

	static double cot(double x)
	{
		return 1 / std::tan(x);
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

	static int precedence(const string& c);
}

namespace nums {
	const double E = 2.71828182845904523536;
	const double PI = 3.141592653589793238463;
	const double LN10 = 0.434294481903251827651;
}

class  AD
{
public:
	typedef  std::shared_ptr<AD> ptr;

	enum class type {
		function,
		operation,
		known,
		unknown,
	};

	static map < string, double (*)(double) > functions;
	static map<string, double (*)(double, double) > operations;

	static map<string, ptr(*)(ptr a, ptr b, const string& x)> doperations;
	static map<string, std::function<ptr(shared_ptr<AD>)>> dfunctions;

	std::variant<std::string, double> value;
	shared_ptr<AD> left = nullptr;
	shared_ptr<AD> right = nullptr;
	type typ;

	AD() {
		value = 0.0;
		typ = AD::type::known;
	}


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

	AD(string x, ptr r) :
		value(x), right(r), typ(type::function)
	{

	}

	AD(string x, ptr l, ptr r) :
		value(x), left(l), right(r), typ(type::operation)
	{
	}

public:
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
				throw VariableNotFound(std::get<string>(value));
			}
			return it->second;
		}
		else if (typ == type::function)
		{
			map<string, double (*)(double)>::const_iterator it = AD::functions.find(std::get<string>(value));
			if (it == AD::functions.end())
			{
				throw FunctionNotFound(std::get<string>(value));
			}
			return (it->second)(right->evaluate(vars));
		}
		else/* (typ == type::operation)*/
		{
			map<string, double (*)(double, double)>::const_iterator it = AD::operations.find(std::get<string>(value));
			if (it == AD::operations.end()) {
				throw OperatorNotFound(std::get<string>(value));
			}
			return (it->second)(left->evaluate(vars), right->evaluate(vars));
		}
	}

	static AD parse(const string& exp)
	{
		list<string> lst = infixToPostfix(exp);

		stack<ptr> stack;

		for (auto e = lst.begin(); e != lst.end(); e++)
		{
			int pre = extrafncs::precedence(*e);
			if (pre > 1 and pre < 5)
			{
				ptr b = stack.top();
				stack.pop();
				ptr a = stack.top();
				stack.pop();
				stack.push(ptr(new AD(*e, a, b)));
			}
			else if (pre == 5)
			{
				ptr x = stack.top();
				stack.pop();
				stack.push(ptr(new AD(*e, x)));
			}
			else {
				if (extrafncs::isNumber(*e))
					stack.push(ptr(new AD(std::stod(*e))));
				else
					stack.push(ptr(new AD(*e)));
			}
		}
		return *stack.top();
	}

	ptr derivative(const string& x)
	{
		if (typ == type::known)
		{
			return getNum(0);
		}
		else if (typ == type::unknown)
		{
			if (std::get<string>(value) == x)
				return getNum(1);
			else
				return getNum(0);
		}
		else if (typ == type::function)
		{
			map<string, std::function<ptr(shared_ptr<AD>)>>::const_iterator it = AD::dfunctions.find(std::get<string>(value));
			if (it == AD::dfunctions.end())
				throw FunctionDerivativeNotFound(std::get<string>(value));

			return (it->second)(right) * right->derivative(x);
		}
		else /*if (typ == type::operation)*/
		{
			std::map<string, ptr(*)(ptr a, ptr b, const string& x)>::const_iterator it = doperations.find(std::get<string>(value));
			if (it == AD::doperations.end())
				throw OperatorDerivativeNotFound(std::get<string>(value));

			return (it->second)(left, right, x);
		}
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
			else if (extrafncs::precedence(*e) >= 2) {
				while (!stack.empty() && extrafncs::precedence(*e) <= extrafncs::precedence(stack.top()))
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

private:
	static  std::list<string> parse_toList(const string& exp)
	{
		string container = "";
		std::list<string> ans;
		for (auto e = exp.begin(); e != exp.end(); e++)
		{
			if (*e == ' ')
				continue;
			if (extrafncs::precedence(*e) >= 0)
			{
				if (container.size() == 0 && *e == '-' && (e == exp.begin() || *(e - 1) != ')'))
				{
					ans.push_back("-1");
					ans.push_back("*");
				}
				else if (container.size() == 0 && *e == '+' && (e == exp.begin() || *(e - 1) != ')'))
				{
					continue;
				}
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

	bool isConst()
	{
		return typ == type::known;
	}

	//-------------------------------------------

	friend ptr operator *(ptr a, ptr b)
	{
		if (a->isConst() and b->isConst())
			return ptr(new AD(std::get<double>(a->value) * std::get<double>(b->value)));
		if (*a == 0 or *b == 0)
			return getNum(0);
		if (*a == 1)
			return b;
		if (*b == 1)
			return a;
		return ptr(new AD("*", a, b));
	}

	friend ptr operator +(ptr a, ptr b)
	{
		if (a->isConst() and b->isConst())
			return ptr(new AD(std::get<double>(a->value) + std::get<double>(b->value)));

		if (*a == 0)
			return b;
		if (*b == 0)
			return a;
		return ptr(new AD("+", a, b));
	}

	friend ptr operator -(ptr a, ptr b)
	{
		if (a->isConst() and b->isConst())
			return ptr(new AD(std::get<double>(a->value) - std::get<double>(b->value)));
		if (*a == 0)
			return getNum(-1) * b;
		if (*b == 0)
			return a;
		return ptr(new AD("-", a, b));
	}

	friend ptr operator /(ptr a, ptr b)
	{
		if (a->isConst() and b->isConst())
			return ptr(new AD(std::get<double>(a->value) / std::get<double>(b->value)));
		if (*a == 0 and !(*b == 0))
		{
			return  getNum(0);
		}
		return ptr(new AD("/", a, b));
	}

	friend ptr operator ^(ptr a, ptr b)
	{
		if (a->isConst() and b->isConst())
			return ptr(new AD(pow(std::get<double>(a->value), std::get<double>(b->value))));
		if (*a == 0 and !(*b == 0))
			return getNum(0);
		if (*b == 0 and !(*a == 0))
			return getNum(1);
		if (*a == 1)
			return getNum(1);
		if (*b == 1)
			return a;
		return ptr(new AD("^", a, b));
	}

	//Functions
	friend ptr cos(ptr a)
	{
		return ptr(new AD("cos", a));
	}

	friend ptr sin(ptr a)
	{
		return ptr(new AD("sin", a));
	}

	friend ptr tan(ptr a)
	{
		return ptr(new AD("tan", a));
	}

	friend ptr sec(ptr a)
	{
		return ptr(new AD("sec", a));
	}

	friend ptr cosec(ptr a)
	{
		return ptr(new AD("cosec", a));
	}

	friend ptr cot(ptr a)
	{
		return ptr(new AD("cot", a));
	}

	friend ptr ln(ptr a)
	{
		return ptr(new AD("ln", a));
	}

	friend ptr log10(ptr a)
	{
		return ptr(new AD("log10", a));
	}

	friend ptr sinh(ptr a)
	{
		return ptr(new AD("sinh", a));
	}

	friend ptr cosh(ptr a)
	{
		return ptr(new AD("cosh", a));
	}

	friend ptr tanh(ptr a)
	{
		return ptr(new AD("tanh", a));
	}


	//d<Functions>
	static ptr dtan(ptr a)
	{
		return cos(a) ^ getNum(-2);
	}

	static ptr dcos(ptr a)
	{
		return getNum(-1) * sin(a);
	}

	static ptr dsin(ptr a)
	{
		return cos(a);
	}

	static ptr dsec(ptr a)
	{
		return cos(a) * tan(a);
	}

	static ptr dcosec(ptr a)
	{
		return getNum(-1) * cosec(a) * cot(a);
	}

	static ptr dcot(ptr a)
	{
		return getNum(-1) * (sin(a) ^ getNum(-2));
	}

	static ptr dlog10(ptr a)
	{
		return getNum(1 / nums::LN10) / a;
	}

	static ptr dln(ptr a)
	{
		return getNum(1) / a;
	}

	static ptr dcosh(ptr a)
	{
		return sinh(a);
	}

	static ptr dsinh(ptr a)
	{
		return cosh(a);
	}

	static ptr dtanh(ptr a)
	{
		return getNum(1) / (cosh(a) ^ getNum(2));
	}

	static ptr darcsin(ptr a)
	{
		return getNum(1) / ((getNum(1) - (a ^ getNum(2))) ^ getNum(0.5));
	}

	static ptr darccos(ptr a)
	{
		return getNum(-1) / ((getNum(1) - (a ^ getNum(2))) ^ getNum(0.5));
	}

	static ptr darctan(ptr a)
	{
		return getNum(1) / (getNum(1) + (a ^ getNum(2)));
	}

private:
	//d<operators>
	static ptr dmul(ptr a, ptr b, const string& x)
	{
		return a * b->derivative(x) + b * a->derivative(x);
	}

	static ptr dadd(ptr a, ptr b, const string& x)
	{
		return a->derivative(x) + b->derivative(x);
	}

	static ptr dsub(ptr a, ptr b, const string& x)
	{
		return a->derivative(x) - b->derivative(x);
	}

	static ptr ddiv(ptr a, ptr b, const string& x)
	{
		return (a->derivative(x) * b - b->derivative(x) * a) / (b ^ getNum(2));
	}

	static ptr dpow(ptr f, ptr g, const string& x)
	{
		if (g->isConst() and not(g == 0))
		{
			return (f ^ getNum(std::get<double>(g->value) - 1)) * getNum(std::get<double>(g->value)) * f->derivative(x);
		}
		return (f ^ g) * (g->derivative(x) * ptr(new AD("ln", f)) + g * f->derivative(x) / f);
	}


public:
	static ptr getNum(double a)
	{
		return std::make_shared<AD>(AD(a));
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

	void replaceUnknown(string p, string n) {
		replaceUnknown(p, n, *this);
	}

	static void replaceUnknown(string p, string n, AD& exp) {
		if (exp.typ == type::unknown)
		{
			if (std::get<string>(exp.value) == p) {
				exp.value = n;
			}
		}
		else if (exp.typ == type::function)
			replaceUnknown(p, n, *exp.right);
		else if (exp.typ == type::operation)
		{
			replaceUnknown(p, n, *exp.right);
			replaceUnknown(p, n, *exp.left);
		}
	}

	AD::ptr copy() {
		if (typ == type::known) {
			return std::make_shared<AD>(AD(std::get<double>(value)));
		}
		else if (typ == type::unknown) {
			return std::make_shared<AD>(AD(std::get<string>(value)));
		}
		else if (typ == type::function) {
			return std::make_shared<AD>(AD(std::get<string>(value), right->copy()));
		}
		else { //if(typ == type::operation) {
			return std::make_shared<AD>(AD(std::get<string>(value), left->copy(), right->copy()));
		}
	}

	static void putVal(const map<string, double>& vals, AD& exp) {
		if (exp.typ == type::unknown)
		{
			std::map<string, double>::const_iterator it = vals.find(std::get<string>(exp.value));
			if (it != vals.end()) {
				exp.value = it->second;
				exp.typ = type::known;
			}
		}
		else if (exp.typ == type::function)
			putVal(vals, *exp.right);
		else if (exp.typ == type::operation)
		{
			putVal(vals, *exp.right);
			putVal(vals, *exp.left);
		}
	}

	void putVal(const map<string, double>& vals) {
		AD::putVal(vals, *this);
	}

};

map<string, double (*)(double) > AD::functions = {
{"sin",std::sin} ,
{"cos",std::cos},
{"tan",std::tan},
{"log10",std::log10},
{"exp",std::exp},
{"ln",std::log},
{"sec",extrafncs::sec},
{"cosec",extrafncs::cosec},
{"cot",extrafncs::cot},
{"sinh",std::sinh},
{"cosh",std::cosh},
{"tanh",std::tanh},
{"arcsin",std::asin},
{"arccos",std::acos},
{"arctan",std::atan} };

map<string, std::shared_ptr<AD>(*)(std::shared_ptr<AD> a, std::shared_ptr<AD> b, const string& x)> AD::doperations = {
{"+",AD::dadd},
{"*",AD::dmul},
{"-",AD::dsub},
{"/",AD::ddiv},
{"^",AD::dpow} };

map<string, double (*)(double, double) > AD::operations = {
{"^", extrafncs::pow} ,
{ "+", extrafncs::add} ,
{"-",extrafncs::sub},
{"*",extrafncs::mul},
{"/",extrafncs::div} };

map<string, std::function< std::shared_ptr<AD>(shared_ptr<AD>)>> AD::dfunctions = {
{"cos",AD::dcos},
{"sin",AD::dsin},
{"tan",AD::dtan},
{"sec",AD::dsec},
{"cosec",AD::dcosec},
{"cot",AD::dcot},
{"ln",AD::dln},
{"sinh",AD::dsinh},
{"cosh",AD::dcosh},
{"tanh",AD::dtanh},
{"arcsin",AD::darcsin},
{"arccos",AD::darccos},
{"arctan",AD::darctan},
{"log10",AD::dlog10} };

static int extrafncs::precedence(const string& c)
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