#pragma once
#include <string>
#include<map>
#include<list>
#include<stack>
#include<functional>
#include"Exceptions.h"
#include"Matrix.h"

using namespace std;

class STRDOUBLE {
	std::string str = "";
	double num = 0;
	bool isItStr = false;

public:
	STRDOUBLE() {

	}
	STRDOUBLE(const std::string& str) :str(str) {
		isItStr = true;
	}

	STRDOUBLE(const double& num) :num(num) {
		isItStr = false;
	}

	STRDOUBLE(const STRDOUBLE& a) {
		str = a.str;
		num = a.num;
		isItStr = a.isItStr;
	}

	bool isStr() const {
		return isItStr;
	}

	bool isNum() const {
		return !isItStr;
	}

	double Num() const {
		return num;
	}

	std::string Str() const {
		return str;
	}

	STRDOUBLE operator=(const std::string& other) {
		str = other;
		isItStr = true;
		num = 0;
		return *this;
	}

	STRDOUBLE operator=(const double& other) {
		num = other;
		isItStr = false;
		str = "";
		return *this;
	}

	STRDOUBLE operator=(const STRDOUBLE& v) {
		str = v.str;
		num = v.num;
		isItStr = v.isItStr;
		return *this;
	}

	bool operator==(const STRDOUBLE& o) const {
		if (o.isStr() and isStr()) {
			return o.str == str;
		}
		if (o.isNum() and isNum()) {
			return o.num == num;
		}
		return false;
	}

	bool operator<(const STRDOUBLE& o) {
		if (o.isStr() and isStr()) {
			return o.str < str;
		}
		if (o.isNum() and isNum()) {
			return o.num < num;
		}
		return true;
	}

	friend ostream& operator<<(ostream& o, const STRDOUBLE& a) {
		if (a.isItStr) o << a.str;
		else o << a.num;
		return o;
	}
};

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

	STRDOUBLE value;
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
			out << exp.value;
		else if (exp.typ == type::unknown)
			out << exp.value;
		else if (exp.typ == type::function)
			out << (exp.value) << '(' << *(exp.right) << ')';
		else if (exp.typ == type::operation)
			out << '(' << *exp.left << exp.value << *exp.right << ')';
		return out;
	}

	double evaluate(const map<string, double>& vars = {})
	{

		if (typ == type::known)
		{
			return value.Num();
		}
		else if (typ == type::unknown)
		{
			std::map<string, double>::const_iterator it = vars.find(value.Str());
			if (it == vars.end())
			{
				throw VariableNotFound(value.Str());
			}
			return it->second;
		}
		else if (typ == type::function)
		{
			map<string, double (*)(double)>::const_iterator it = AD::functions.find(value.Str());
			if (it == AD::functions.end())
			{
				throw FunctionNotFound(value.Str());
			}
			return (it->second)(right->evaluate(vars));
		}
		else/* (typ == type::operation)*/
		{
			map<string, double (*)(double, double)>::const_iterator it = AD::operations.find(value.Str());
			if (it == AD::operations.end()) {
				throw OperatorNotFound(value.Str());
			}
			return (it->second)(left->evaluate(vars), right->evaluate(vars));
		}
	}

	static AD::ptr parse(const string& exp)
	{
		list<string> lst = infixToPostfix(exp);

		stack<ptr> stack;

		for (auto e = lst.begin(); e != lst.end(); e++)
		{
			int pre = extrafncs::precedence(*e);
			if (pre > 1 and pre < 6)
			{
				ptr b = stack.top();
				stack.pop();
				ptr a = stack.top();
				stack.pop();
				stack.push(ptr(new AD(*e == "" ? "*" : *e, a, b)));
			}
			else if (pre >= 1000)
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
		return stack.top();
	}

	ptr derivative(const string& x) const
	{
		if (typ == type::known)
		{
			return getNum(0);
		}
		else if (typ == type::unknown)
		{
			if (value.Str() == x)
				return getNum(1);
			else
				return getNum(0);
		}
		else if (typ == type::function)
		{
			map<string, std::function<ptr(shared_ptr<AD>)>>::const_iterator it = AD::dfunctions.find(value.Str());
			if (it == AD::dfunctions.end())
				throw FunctionDerivativeNotFound(value.Str());

			return (it->second)(right) * right->derivative(x);
		}
		else /*if (typ == type::operation)*/
		{
			std::map<string, ptr(*)(ptr a, ptr b, const string& x)>::const_iterator it = doperations.find(value.Str());
			if (it == AD::doperations.end())
				throw OperatorDerivativeNotFound(value.Str());

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

	static AD::ptr getptr(AD&& a) {
		return make_shared<AD>(a);
	}

	//private:
	static  std::list<string> parse_toList(const string& exp)
	{
		//temporary container for various variables, numbers or functions
		string container = "";

		//what we have to return
		std::list<string> ans;
		for (auto e = exp.begin(); e != exp.end(); e++)
		{
			if (*e == ' ')
				continue;
			//precedence>=0 means its a bracket or an operation like *,/,etc
			if (extrafncs::precedence(*e) >= 0)
			{
				//if container is empty and *e is '-' we will know that, that '-' is used as unary operator
				if (container.size() == 0 && *e == '-' && (e == exp.begin() || *(e - 1) != ')'))
				{
					//we multiply by -1
					ans.push_back("-1");
					ans.push_back("");
				}
				//if above condition and if its a '+' then just ignore
				else if (container.size() == 0 && *e == '+' && (e == exp.begin() || *(e - 1) != ')'))
				{
					continue;
				}
				else
				{
					//add to the end of the list
					if (container.size() != 0)
						ans.push_back(container);
					ans.push_back(string(1, *e));
				}
				//clear container for next iterations
				container.clear();
			}
			else {
				container.push_back(*e);
			}
		}
		//if container is not empty then push the remaining at the end of the ans
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
			return ptr(new AD((a->value).Num() * (b->value).Num()));
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
			return ptr(new AD((a->value).Num() + (b->value).Num()));

		if (*a == 0)
			return b;
		if (*b == 0)
			return a;
		return ptr(new AD("+", a, b));
	}

	friend ptr operator -(ptr a, ptr b)
	{
		if (a->isConst() and b->isConst())
			return ptr(new AD((a->value).Num() - (b->value).Num()));
		if (*a == 0)
			return getNum(-1) * b;
		if (*b == 0)
			return a;
		return ptr(new AD("-", a, b));
	}

	friend ptr operator /(ptr a, ptr b)
	{
		if (a->isConst() and b->isConst())
			return ptr(new AD((a->value).Num() / (b->value).Num()));
		if (*a == 0 and !(*b == 0))
		{
			return  getNum(0);
		}
		return ptr(new AD("/", a, b));
	}

	friend ptr operator ^(ptr a, ptr b)
	{
		if (a->isConst() and b->isConst())
			return ptr(new AD(pow((a->value).Num(), (b->value).Num())));
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
		return sec(a) * tan(a);
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
			return (f ^ getNum((g->value).Num() - 1)) * getNum((g->value).Num()) * f->derivative(x);
		}
		return (f ^ g) * (g->derivative(x) * ptr(new AD("ln", f)) + g * f->derivative(x) / f);
	}


public:
	static ptr getNum(double a)
	{
		return std::make_shared<AD>(AD(a));
	}

	bool operator ==(double b) const
	{
		if (typ == type::known)
		{
			if ((value).Num() == b)
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
			if (exp.value.Str() == p) {
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
			return std::make_shared<AD>(AD(value.Num()));
		}
		else if (typ == type::unknown) {
			return std::make_shared<AD>(AD(value.Str()));
		}
		else if (typ == type::function) {
			return std::make_shared<AD>(AD(value.Str(), right->copy()));
		}
		else { //if(typ == type::operation) {
			return std::make_shared<AD>(AD(value.Str(), left->copy(), right->copy()));
		}
	}

	static void putVal(const map<string, double>& vec, AD& exp) {
		if (exp.typ == type::unknown)
		{
			std::map<string, double>::const_iterator it = vec.find((exp.value).Str());
			if (it != vec.end()) {
				exp.value = it->second;
				exp.typ = type::known;
			}
		}
		else if (exp.typ == type::function)
			putVal(vec, *exp.right);
		else if (exp.typ == type::operation)
		{
			putVal(vec, *exp.right);
			putVal(vec, *exp.left);
		}
	}

	void putVal(const map<string, double>& vec) {
		AD::putVal(vec, *this);
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
		return 1000;
	if (c == "")
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