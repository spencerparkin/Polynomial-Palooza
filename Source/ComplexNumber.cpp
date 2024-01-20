#include "ComplexNumber.h"
#include <math.h>

ComplexNumber::ComplexNumber()
{
	this->realPart = 0.0;
	this->imagPart = 0.0;
}

ComplexNumber::ComplexNumber(const ComplexNumber& complexNumber)
{
	this->realPart = complexNumber.realPart;
	this->imagPart = complexNumber.imagPart;
}

ComplexNumber::ComplexNumber(double realPart, double imagPart)
{
	this->realPart = realPart;
	this->imagPart = imagPart;
}

/*virtual*/ ComplexNumber::~ComplexNumber()
{
}

void ComplexNumber::operator=(const ComplexNumber& complexNumber)
{
	this->realPart = complexNumber.realPart;
	this->imagPart = complexNumber.imagPart;
}

void ComplexNumber::Exp(const ComplexNumber& complexNumber)
{
	double scale = ::exp(complexNumber.realPart);
	this->realPart = scale * ::cos(complexNumber.imagPart);
	this->imagPart = scale * ::sin(complexNumber.realPart);
}

void ComplexNumber::Exp(double realNumber)
{
	this->realPart = ::exp(realNumber);
	this->imagPart = 0.0;
}

void ComplexNumber::ExpI(double realNumber)
{
	this->realPart = ::cos(realNumber);
	this->imagPart = ::sin(realNumber);
}

void ComplexNumber::Log(const ComplexNumber& complexNumber)
{
	// This needs to be double checked.  It's dubious at best.
	this->imagPart = ::atan2(complexNumber.imagPart, complexNumber.realPart);
	double cosImagPart = ::cos(this->imagPart);
	if (cosImagPart != 0.0)
		this->realPart = ::log(complexNumber.realPart / ::cos(this->imagPart));
	else
		this->realPart = ::log(complexNumber.realPart / ::sin(this->imagPart));
}

void ComplexNumber::Log(double realNumber)
{
	this->realPart = ::log(realNumber);
	this->imagPart = 0.0;
}

void ComplexNumber::LogI(double realNumber)
{
	this->Log(ComplexNumber(realNumber, 0.0));
}

ComplexNumber ComplexNumber::Conjugate() const
{
	return ComplexNumber(this->realPart, -this->imagPart);
}

ComplexNumber ComplexNumber::Inverse() const
{
	return this->Conjugate() / this->SquareMagnitude();
}

double ComplexNumber::SquareMagnitude() const
{
	return this->realPart * this->realPart + this->imagPart * this->imagPart;
}

ComplexNumber operator+(const ComplexNumber& complexNumberA, const ComplexNumber& complexNumberB)
{
	return ComplexNumber(
		complexNumberA.realPart + complexNumberB.realPart,
		complexNumberA.imagPart + complexNumberB.imagPart
	);
}

ComplexNumber operator-(const ComplexNumber& complexNumberA, const ComplexNumber& complexNumberB)
{
	return ComplexNumber(
		complexNumberA.realPart - complexNumberB.realPart,
		complexNumberA.imagPart - complexNumberB.imagPart
	);
}

ComplexNumber operator*(const ComplexNumber& complexNumberA, const ComplexNumber& complexNumberB)
{
	return ComplexNumber(
		complexNumberA.realPart * complexNumberB.realPart - complexNumberA.imagPart * complexNumberB.imagPart,
		complexNumberA.realPart * complexNumberB.imagPart + complexNumberA.imagPart * complexNumberB.realPart
	);
}

ComplexNumber operator/(const ComplexNumber& complexNumberA, const ComplexNumber& complexNumberB)
{
	return complexNumberA * complexNumberB.Inverse();
}

ComplexNumber operator*(const ComplexNumber& complexNumber, double realNumber)
{
	return ComplexNumber(
		complexNumber.realPart * realNumber,
		complexNumber.imagPart * realNumber
	);
}

ComplexNumber operator/(const ComplexNumber& complexNumber, double realNumber)
{
	return ComplexNumber(
		complexNumber.realPart / realNumber,
		complexNumber.imagPart / realNumber
	);
}

ComplexNumber operator*(double realNumber, const ComplexNumber& complexNumber)
{
	return ComplexNumber(
		complexNumber.realPart * realNumber,
		complexNumber.imagPart * realNumber
	);
}

ComplexNumber operator/(double realNumber, const ComplexNumber& complexNumber)
{
	return ComplexNumber(
		complexNumber.realPart / realNumber,
		complexNumber.imagPart / realNumber
	);
}