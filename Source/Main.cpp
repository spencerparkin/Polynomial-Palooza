#include "ComplexNumber.h"
#include "Polynomial.h"
#include <iostream>

int main(int argc, char** argv)
{
	Polynomial polynomialA, polynomialB;

	polynomialA.coefficientArray.push_back(ComplexNumber(-3.0, 0.0));
	polynomialA.coefficientArray.push_back(ComplexNumber(4.0, 0.0));
	polynomialA.coefficientArray.push_back(ComplexNumber(-7.0, 0.0));

	polynomialB.coefficientArray.push_back(ComplexNumber(2.0, 0.0));
	polynomialB.coefficientArray.push_back(ComplexNumber(-9.0, 0.0));
	polynomialB.coefficientArray.push_back(ComplexNumber(-2.0, 0.0));
	polynomialB.coefficientArray.push_back(ComplexNumber(12.0, 0.0));

	std::cout << "Polynomial A: " << std::string(polynomialA) << std::endl;
	std::cout << "Polynomial B: " << std::string(polynomialB) << std::endl;

	Polynomial product = polynomialA * polynomialB;

	std::cout << "Product: " << std::string(product) << std::endl;

	return 0;
}