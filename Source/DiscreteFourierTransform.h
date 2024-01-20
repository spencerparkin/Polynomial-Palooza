#pragma once

#include "ComplexNumber.h"
#include <vector>
#include <string>

class Polynomial;

class DFT
{
public:
	DFT();
	virtual ~DFT();

	operator std::string() const;

	virtual bool FromPolynomial(const Polynomial& polynomial, std::string& error);
	virtual bool ToPolynomial(Polynomial& polynomial, std::string& error) const;

	// If there are N points here, then they were found by evaluating
	// the polynomial at the N complex Nth roots of unity.
	std::vector<ComplexNumber> pointArray;
};

class FFT : public DFT
{
public:
	FFT();
	virtual ~FFT();

	virtual bool FromPolynomial(const Polynomial& polynomial, std::string& error) override;
	virtual bool ToPolynomial(Polynomial& polynomial, std::string& error) const override;
};