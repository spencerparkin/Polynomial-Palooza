#include "DiscreteFourierTransform.h"
#include "Polynomial.h"
#include <math.h>

//----------------------------- DFT -----------------------------

DFT::DFT()
{
}

/*virtual*/ DFT::~DFT()
{
}

DFT::operator std::string() const
{
	std::string pointArrayStr;

	for (uint32_t i = 0; i < this->pointArray.size(); i++)
	{
		const ComplexNumber& point = this->pointArray[i];

		if (pointArrayStr.length() > 0)
			pointArrayStr += ", ";

		pointArrayStr += std::string(point);
	}

	return pointArrayStr;
}

/*virtual*/ bool DFT::FromPolynomial(const Polynomial& polynomial, std::string& error)
{
	uint32_t degreeBound = polynomial.GetDegreeBound();
	this->pointArray.clear();

	for (uint32_t i = 0; i < degreeBound; i++)
	{
		ComplexNumber rootOfUnity;
		rootOfUnity.ExpI(2.0 * M_PI * double(i) / double(degreeBound));
		this->pointArray.push_back(polynomial.Evaluate(rootOfUnity));
	}

	return true;
}

/*virtual*/ bool DFT::ToPolynomial(Polynomial& polynomial, std::string& error) const
{
	// Crap.  Here we need to invert an NxN Vandermonde matrix.  :/
	error = "Not yet implimented.";
	return false;
}

//----------------------------- FFT -----------------------------

FFT::FFT()
{
}

/*virtual*/ FFT::~FFT()
{
}

/*virtual*/ bool FFT::FromPolynomial(const Polynomial& polynomial, std::string& error)
{
	this->pointArray.clear();

	uint32_t degreeBound = polynomial.GetDegreeBound();
	if (degreeBound == 1)
	{
		this->pointArray.push_back(polynomial[0]);
		return true;
	}

	if (degreeBound % 2 != 0)
	{
		error = "Degree-bound must be even if not one.  Initially, it must be a power of 2.";
		return false;
	}

	Polynomial evenPoly, oddPoly;

	for (uint32_t i = 0; i < degreeBound; i++)
	{
		const ComplexNumber& coefficient = polynomial[i];
		if (i % 2 == 0)
			evenPoly.coefficientArray.push_back(coefficient);
		else
			oddPoly.coefficientArray.push_back(coefficient);
	}

	FFT evenFft, oddFft;

	if (!evenFft.FromPolynomial(evenPoly, error))
		return false;

	if (!oddFft.FromPolynomial(oddPoly, error))
		return false;

	for (uint32_t i = 0; i < degreeBound; i++)
	{
		ComplexNumber rootOfUnity;
		rootOfUnity.ExpI(2.0 * M_PI * double(i) / double(degreeBound));
		uint32_t j = i % (degreeBound / 2);
		this->pointArray.push_back(evenFft.pointArray[j] + rootOfUnity * oddFft.pointArray[j]);
	}

	return true;
}

/*virtual*/ bool FFT::ToPolynomial(Polynomial& polynomial, std::string& error) const
{
	error = "Not yet implimented.";
	return false;
}