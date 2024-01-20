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

	std::cout << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << std::endl;

	Polynomial product = polynomialA * polynomialB;

	std::cout << "Product: " << std::string(product) << std::endl;

	std::cout << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << std::endl;

	polynomialA.SetDegreeBound(4);
	polynomialB.SetDegreeBound(4);

	Polynomial dftA, dftB;

	dftA.DFT(polynomialA, false);
	dftB.DFT(polynomialB, false);

	std::cout << "DFT A: " << std::string(dftA) << std::endl;
	std::cout << "DFT B: " << std::string(dftB) << std::endl;

	std::cout << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << std::endl;

	Polynomial invDftA, invDftB;

	invDftA.DFT(dftA, true);
	invDftB.DFT(dftB, true);

	std::cout << "inv-DFT A: " << std::string(invDftA) << std::endl;
	std::cout << "inv-DFT B: " << std::string(invDftB) << std::endl;

	std::cout << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << std::endl;

	Polynomial fftA, fftB;

	fftA.FFT(polynomialA, false);
	fftB.FFT(polynomialB, false);

	std::cout << "FFT A: " << std::string(fftA) << std::endl;
	std::cout << "FFT B: " << std::string(fftB) << std::endl;

	std::cout << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << std::endl;

	Polynomial invFftA, invFftB;

	invFftA.FFT(fftA, true);
	invFftB.FFT(fftB, true);

	std::cout << "inv-DFT A: " << std::string(invFftA) << std::endl;
	std::cout << "inv-DFT B: " << std::string(invFftB) << std::endl;

	std::cout << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	std::cout << std::endl;

	Polynomial fastProduct;

	fastProduct.FastMultiply(polynomialA, polynomialB);

	std::cout << "Fast product: " << std::string(fastProduct) << std::endl;

	return 0;
}