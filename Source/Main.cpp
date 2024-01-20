#include "ComplexNumber.h"
#include "Polynomial.h"
#include "DiscreteFourierTransform.h"
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

	std::cout << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << std::endl;

	std::string error;

	DFT dftA, dftB;
	dftA.FromPolynomial(polynomialA, error);
	dftB.FromPolynomial(polynomialB, error);

	std::cout << "DFT A: " << std::string(dftA) << std::endl;
	std::cout << "DFT B: " << std::string(dftB) << std::endl;

	std::cout << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << std::endl;

	FFT fftA, fftB;
	fftA.FromPolynomial(polynomialA, error);
	fftB.FromPolynomial(polynomialB, error);

	std::cout << "FFT A: " << std::string(fftA) << std::endl;
	std::cout << "FFT B: " << std::string(fftB) << std::endl;

	return 0;
}