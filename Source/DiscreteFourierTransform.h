#pragma once

#include "ComplexNumber.h"
#include <vector>

class Polynomial;

class DFT
{
public:
	DFT();
	virtual ~DFT();

	operator std::string() const;

	virtual void FromPolynomial(const Polynomial& polynomial);
	virtual void ToPolynomial(Polynomial& polynomial) const;

	// If there are N points here, then they were found by evaluating
	// the polynomial at the N complex Nth roots of unity.
	std::vector<ComplexNumber> pointArray;
};

class FFT : public DFT
{
public:
	FFT();
	virtual ~FFT();

	virtual void FromPolynomial(const Polynomial& polynomial) override;
	virtual void ToPolynomial(Polynomial& polynomial) const override;
};