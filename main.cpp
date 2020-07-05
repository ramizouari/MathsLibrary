#include "integer.h"
#include "guassian_integer.h"
#include "polynomial.h"
#include "ring.h"
#include <iostream>

using namespace std;

int main()
{
	guassian_integer a(2,3);
	polynomial<integer> p({ 1,1,1 }), q({ 1,1 });
	cout << spow(q,10);
	return false;
}