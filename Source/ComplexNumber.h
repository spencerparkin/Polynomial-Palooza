#pragma once

class ComplexNumber
{
public:
	ComplexNumber();
	ComplexNumber(const ComplexNumber& complexNumber);
	ComplexNumber(double realPart, double imagPart);
	virtual ~ComplexNumber();

	void operator=(const ComplexNumber& complexNumber);
	void operator+=(const ComplexNumber& complexNumber);
	void operator-=(const ComplexNumber& complexNumber);
	void operator*=(const ComplexNumber& complexNumber);
	void operator/=(const ComplexNumber& complexNumber);
	bool operator==(const ComplexNumber& complexNumber) const;
	bool operator!=(const ComplexNumber& complexNumber) const;

	void Exp(const ComplexNumber& complexNumber);
	void Exp(double realNumber);
	void ExpI(double realNumber);

	void Log(const ComplexNumber& complexNumber);
	void Log(double realNumber);
	void LogI(double realNumber);

	ComplexNumber Conjugate() const;
	ComplexNumber Inverse() const;

	double SquareMagnitude() const;

	double realPart, imagPart;
};

ComplexNumber operator+(const ComplexNumber& complexNumberA, const ComplexNumber& complexNumberB);
ComplexNumber operator-(const ComplexNumber& complexNumberA, const ComplexNumber& complexNumberB);
ComplexNumber operator*(const ComplexNumber& complexNumberA, const ComplexNumber& complexNumberB);
ComplexNumber operator/(const ComplexNumber& complexNumberA, const ComplexNumber& complexNumberB);
ComplexNumber operator*(const ComplexNumber& complexNumber, double realNumber);
ComplexNumber operator/(const ComplexNumber& complexNumber, double realNumber);
ComplexNumber operator*(double realNumber, const ComplexNumber& complexNumber);
ComplexNumber operator/(double realNumber, const ComplexNumber& complexNumber);