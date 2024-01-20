#pragma once

#include "ComplexNumber.h"
#include <vector>

class Polynomial
{
public:
	Polynomial();
	Polynomial(const Polynomial& polynomial);
	virtual ~Polynomial();

	void operator=(const Polynomial& polynomial);
	const ComplexNumber& operator[](uint32_t i) const;
	ComplexNumber& operator[](uint32_t i);

	uint32_t Degree() const;
	uint32_t GetDegreeBound() const;
	void SetDegreeBound(uint32_t degreeBound) const;
	void SetDegreeBound(uint32_t degreeBound);

	ComplexNumber Evaluate(const ComplexNumber& complexArg) const;

	std::vector<ComplexNumber> coefficientArray;
};

Polynomial operator+(const Polynomial& polynomialA, const Polynomial& polynomialB);
Polynomial operator-(const Polynomial& polynomialA, const Polynomial& polynomialB);
Polynomial operator*(const Polynomial& polynomialA, const Polynomial& polynomialB);