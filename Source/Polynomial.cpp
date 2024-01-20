#include "Polynomial.h"
#include <algorithm>
#include <format>

Polynomial::Polynomial()
{
}

Polynomial::Polynomial(const Polynomial& polynomial)
{
	for (const ComplexNumber& coefficient : polynomial.coefficientArray)
		this->coefficientArray.push_back(coefficient);
}

/*virtual*/ Polynomial::~Polynomial()
{
}

void Polynomial::operator=(const Polynomial& polynomial)
{
	this->coefficientArray.clear();
	for (const ComplexNumber& coefficient : polynomial.coefficientArray)
		this->coefficientArray.push_back(coefficient);
}

const ComplexNumber& Polynomial::operator[](uint32_t i) const
{
	return this->coefficientArray[i];
}

ComplexNumber& Polynomial::operator[](uint32_t i)
{
	return this->coefficientArray[i];
}

Polynomial::operator std::string() const
{
	std::string polynomialStr;

	uint32_t degree = this->Degree();
	for (int32_t i = degree - 1; i >= 0; i--)
	{
		const ComplexNumber& coefficient = this->coefficientArray[i];

		if (polynomialStr.length() > 0)
			polynomialStr += " + ";

		polynomialStr += "(" + std::string(coefficient) + ")" + std::format("x^{}", i);
	}

	return polynomialStr;
}

uint32_t Polynomial::Degree() const
{
	uint32_t degree = 0;

	while (degree < this->coefficientArray.size())
	{
		const ComplexNumber& coefficient = this->coefficientArray[degree];
		if (coefficient.SquareMagnitude() == 0.0)
			break;

		degree++;
	}

	return degree;
}

uint32_t Polynomial::GetDegreeBound() const
{
	return (uint32_t)this->coefficientArray.size();
}

void Polynomial::SetDegreeBound(uint32_t degreeBound) const
{
	const_cast<Polynomial*>(this)->SetDegreeBound(degreeBound);
}

void Polynomial::SetDegreeBound(uint32_t degreeBound)
{
	while (this->coefficientArray.size() < degreeBound)
		this->coefficientArray.push_back(ComplexNumber(0.0, 0.0));

	while (this->coefficientArray.size() > degreeBound)
	{
		const ComplexNumber& coefficient = this->coefficientArray[this->coefficientArray.size() - 1];
		if (coefficient == ComplexNumber(0.0, 0.0))
			this->coefficientArray.pop_back();
		else
			break;
	}
}

ComplexNumber Polynomial::Evaluate(const ComplexNumber& complexArg) const
{
	ComplexNumber result(0.0, 0.0);
	ComplexNumber complexArgPower(1.0, 0.0);

	for (uint32_t i = 0; i < this->coefficientArray.size(); i++)
	{
		const ComplexNumber& coefficient = this->coefficientArray[i];

		result += complexArgPower * coefficient;
		complexArgPower *= complexArg;
	}

	return result;
}

/*static*/ uint32_t Polynomial::SmallestPowerOfTwoGreaterThanOrEqualTo(uint32_t givenInt)
{
	if (givenInt == 0 || (givenInt & (givenInt - 1)) == 0)
		return givenInt;

	uint32_t powerOfTwo = 1;
	while (powerOfTwo < givenInt)
		powerOfTwo <<= 1;

	return powerOfTwo;
}

void Polynomial::FastMultiply(const Polynomial& polynomialA, const Polynomial& polynomialB)
{
	uint32_t maxDegree = std::max(polynomialA.Degree(), polynomialB.Degree());
	uint32_t degreeBound = 2 * SmallestPowerOfTwoGreaterThanOrEqualTo(maxDegree);

	polynomialA.SetDegreeBound(degreeBound);
	polynomialB.SetDegreeBound(degreeBound);

	Polynomial dftA, dftB;

#define USE_FFT

#if defined USE_FFT
	dftA.FFT(polynomialA, false);
	dftB.FFT(polynomialB, false);
#else
	dftA.DFT(polynomialA, false);
	dftB.DFT(polynomialB, false);
#endif

	Polynomial dftProduct;
	dftProduct.SetDegreeBound(degreeBound);
	for (uint32_t i = 0; i < degreeBound; i++)
		dftProduct[i] = dftA[i] * dftB[i];

#if defined USE_FFT
	this->FFT(dftProduct, true);
#else
	this->DFT(dftProduct, true);
#endif
}

bool Polynomial::DFT(const Polynomial& polynomial, bool inverse)
{
	uint32_t degreeBound = polynomial.GetDegreeBound();
	this->coefficientArray.clear();

	for (uint32_t i = 0; i < degreeBound; i++)
	{
		ComplexNumber rootOfUnity;
		rootOfUnity.ExpI((inverse ? -1.0 : 1.0) * 2.0 * M_PI * double(i) / double(degreeBound));
		ComplexNumber coefficient = polynomial.Evaluate(rootOfUnity);
		if (inverse)
			coefficient /= double(degreeBound);
		this->coefficientArray.push_back(coefficient);
	}

	return true;
}

// TODO: The inverse version of this is wrong and needs to be fixed.
bool Polynomial::FFT(const Polynomial& polynomial, bool inverse, bool recursed /*= false*/)
{
	this->coefficientArray.clear();

	uint32_t degreeBound = polynomial.GetDegreeBound();
	if (degreeBound == 1)
	{
		this->coefficientArray.push_back(polynomial[0]);
		return true;
	}

	if (degreeBound % 2 != 0)
		return false;

	Polynomial evenPoly, oddPoly;

	for (uint32_t i = 0; i < degreeBound; i++)
	{
		const ComplexNumber& coefficient = polynomial[i];
		if (i % 2 == 0)
			evenPoly.coefficientArray.push_back(coefficient);
		else
			oddPoly.coefficientArray.push_back(coefficient);
	}

	Polynomial evenFFT, oddFFT;

	if (!evenFFT.FFT(evenPoly, inverse, true))
		return false;

	if (!oddFFT.FFT(oddPoly, inverse, true))
		return false;

	for (uint32_t i = 0; i < degreeBound; i++)
	{
		ComplexNumber rootOfUnity;
		rootOfUnity.ExpI((inverse ? -1.0 : 1.0) * 2.0 * M_PI * double(i) / double(degreeBound));
		uint32_t j = i % (degreeBound / 2);
		ComplexNumber coefficient = evenFFT.coefficientArray[j] + rootOfUnity * oddFFT.coefficientArray[j];
		if (!recursed && inverse)
			coefficient /= double(degreeBound);
		this->coefficientArray.push_back(coefficient);
	}

	return true;
}

void Polynomial::Trim()
{
	constexpr double eps = 1e-5;

	for (ComplexNumber& complexNumber : this->coefficientArray)
	{
		if (::fabs(complexNumber.realPart) < eps)
			complexNumber.realPart = 0.0;

		if (::fabs(complexNumber.imagPart) < eps)
			complexNumber.imagPart = 0.0;
	}

	while (this->coefficientArray.size() > 1)
	{
		ComplexNumber& complexNumber = this->coefficientArray[this->coefficientArray.size() - 1];
		double squareMag = complexNumber.SquareMagnitude();
		if (squareMag < eps)
			this->coefficientArray.pop_back();
		else
			break;
	}
}

Polynomial operator+(const Polynomial& polynomialA, const Polynomial& polynomialB)
{
	Polynomial sum;

	sum.SetDegreeBound(std::max(polynomialA.Degree(), polynomialB.Degree()));
	polynomialA.SetDegreeBound(sum.GetDegreeBound());
	polynomialB.SetDegreeBound(sum.GetDegreeBound());

	for (uint32_t i = 0; i < sum.GetDegreeBound(); i++)
		sum[i] = polynomialA[i] + polynomialB[i];

	return sum;
}

Polynomial operator-(const Polynomial& polynomialA, const Polynomial& polynomialB)
{
	Polynomial diff;

	diff.SetDegreeBound(std::max(polynomialA.Degree(), polynomialB.Degree()));
	polynomialA.SetDegreeBound(diff.GetDegreeBound());
	polynomialB.SetDegreeBound(diff.GetDegreeBound());

	for (uint32_t i = 0; i < diff.GetDegreeBound(); i++)
		diff[i] = polynomialA[i] - polynomialB[i];

	return diff;
}

Polynomial operator*(const Polynomial& polynomialA, const Polynomial& polynomialB)
{
	Polynomial product;
	product.SetDegreeBound(polynomialA.Degree() + polynomialB.Degree());

	// This is the naive method of polynomial multiplication.  A goal of this
	// program is to see if we can do this in O(N ln N) time using an FFT.
	// As a check on such an optimization, we provide this method also.
	for (uint32_t i = 0; i < polynomialA.GetDegreeBound(); i++)
	{
		const ComplexNumber& coefficientA = polynomialA[i];

		for (uint32_t j = 0; j < polynomialB.GetDegreeBound(); j++)
		{
			const ComplexNumber& coefficientB = polynomialB[j];

			product.coefficientArray[i + j] += coefficientA * coefficientB;
		}
	}

	return product;
}