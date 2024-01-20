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

	uint32_t maxDegree = std::max(polynomialA.Degree(), polynomialB.Degree());
	product.SetDegreeBound(2 * maxDegree);
	polynomialA.SetDegreeBound(maxDegree);
	polynomialB.SetDegreeBound(maxDegree);

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