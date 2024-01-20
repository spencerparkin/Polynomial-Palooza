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

/*virtual*/ void DFT::FromPolynomial(const Polynomial& polynomial)
{
	uint32_t degreeBound = polynomial.GetDegreeBound();
	this->pointArray.clear();

	for (uint32_t i = 0; i < degreeBound; i++)
	{
		ComplexNumber rootOfUnity;
		rootOfUnity.ExpI(2.0 * M_PI * double(i) / double(degreeBound));
		this->pointArray.push_back(polynomial.Evaluate(rootOfUnity));
	}
}

/*virtual*/ void DFT::ToPolynomial(Polynomial& polynomial) const
{
	// Crap.  Here we need to invert an NxN Vandermonde matrix.  :/
}

//----------------------------- FFT -----------------------------

FFT::FFT()
{
}

/*virtual*/ FFT::~FFT()
{
}

/*virtual*/ void FFT::FromPolynomial(const Polynomial& polynomial)
{
}

/*virtual*/ void FFT::ToPolynomial(Polynomial& polynomial) const
{
}