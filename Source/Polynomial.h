#pragma once

#include "ComplexNumber.h"
#include <vector>
#include <string>

class Polynomial
{
public:
	Polynomial();
	Polynomial(const Polynomial& polynomial);
	virtual ~Polynomial();

	void operator=(const Polynomial& polynomial);
	const ComplexNumber& operator[](uint32_t i) const;
	ComplexNumber& operator[](uint32_t i);
	operator std::string() const;

	uint32_t Degree() const;
	uint32_t GetDegreeBound() const;
	void SetDegreeBound(uint32_t degreeBound) const;
	void SetDegreeBound(uint32_t degreeBound);

	bool DFT(const Polynomial& polynomial, bool inverse);
	bool FFT(const Polynomial& polynomial, bool inverse, bool recursed = false);

	void FastMultiply(const Polynomial& polynomialA, const Polynomial& polynomialB);

	static uint32_t SmallestPowerOfTwoGreaterThanOrEqualTo(uint32_t givenInt);

	ComplexNumber Evaluate(const ComplexNumber& complexArg) const;

	void Trim();

	std::vector<ComplexNumber> coefficientArray;
};

Polynomial operator+(const Polynomial& polynomialA, const Polynomial& polynomialB);
Polynomial operator-(const Polynomial& polynomialA, const Polynomial& polynomialB);
Polynomial operator*(const Polynomial& polynomialA, const Polynomial& polynomialB);