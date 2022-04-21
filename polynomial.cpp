#include "polynomial.h"

#include <iostream>
#include <string>
//#include <sstream>
#include <iomanip>
#include <cfloat>
#include <math.h>

using std::istream;
using std::ostream;
using std::string;
using std::stringstream;
using std::fixed;
using std::setprecision;
using std::showpos;
using std::endl;
using std::cout;

Polynomial::Polynomial(size_t degree) : _degree(degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = 0.0;
	}
}
Polynomial::Polynomial(size_t degree, const float* coefficients): _degree(degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = coefficients[i];
	}
}
Polynomial::Polynomial(const Polynomial& polynomial): _degree(polynomial._degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = polynomial._coefficients[i];
	}
}
Polynomial::~Polynomial(){

    delete[] _coefficients;

	// DO THIS FIRST TO PREVENT MEMORY LEAKS!




}
const Polynomial Polynomial::Sum(const Polynomial& rhs)const {
    //cout << this->ToString() << endl;

    if (this->_degree > rhs._degree) {
        //makes the polynomial with highest degree of this
        Polynomial retVal(*this);
        for (size_t i = 0; i < rhs._degree + 1; ++i) {
            retVal._coefficients[i] +=  rhs._coefficients[i];
        }
        //cout << retVal.ToString() << endl;
        return retVal;
    }

    else {
        //makes the polynomial with highest degree of rhs
        Polynomial retVal(rhs);
        for (size_t i = 0; i < this->_degree + 1; ++i) {
            retVal._coefficients[i] +=  this->_coefficients[i];
        }
        //cout << retVal.ToString() << endl;
        return retVal;
    }
}


const Polynomial Polynomial::Subtract(const Polynomial& rhs)const{
    //cout << this->ToString() << endl;

    if (this->_degree > rhs._degree) {
        //makes the polynomial with highest degree of this
        Polynomial retVal(*this);
        for (size_t i = 0; i < rhs._degree + 1; ++i) {
            retVal._coefficients[i] -=  rhs._coefficients[i];
        }
        //cout << retVal.ToString() << endl;
        return retVal;
    }

    else {
        //makes the polynomial with highest degree of rhs
        Polynomial retVal(rhs);
        for (size_t i = 0; i < this->_degree + 1; ++i) {
            retVal._coefficients[i] -=  this->_coefficients[i];
        }
        //cout << retVal.ToString() << endl;
        return retVal;
    }
}



const Polynomial Polynomial::Minus()const{
	Polynomial retVal(*this);
	for (size_t i = 0; i < _degree + 1; i++) {
		retVal._coefficients[i] *= -1;
	}
	return retVal;
}





const Polynomial Polynomial::Multiply(const Polynomial& rhs)const {
    //cout << this->ToString() << endl;

    //makes the polynomial with highest degree of this+rhs
        Polynomial retVal(this->_degree+rhs._degree);
        for (size_t i = 0; i < this->_degree + 1; ++i) {
            for (size_t j = 0; j < rhs._degree+1; ++j) {
                retVal._coefficients[i + j] += this->_coefficients[i] * rhs._coefficients[j];
            }
        }
        //cout << retVal.ToString() << endl;
        return retVal;

}



const Polynomial Polynomial::Divide(const Polynomial& rhs)const{
    Polynomial retVal(_degree-1);

    //NOT WORKING
//    Polynomial tmpVal(_degree-1);
//    Polynomial tmpPoly(_degree-1);
//    Polynomial thisCopyPoly(*this);
//
//    for(size_t i=_degree; i>0; --i){
//        for(size_t j=rhs._degree; i>0; --i){

//divides first cofficients of both equation and stores quotient in tmp poly
//            tmpVal._coefficients[i+j] = _coefficients[i]/rhs._coefficients[j];

//adds that qoutient to the actuall queintent
//            retVal._coefficients[i+j] += tmpVal._coefficients[i+j];

//multiply the second equation with the tmp quotient
//            tmpPoly = rhs.Multiply(tmpVal);

//subtracts the tmp polynomial with this and creates a new this with the output
//            thisCopyPoly = thisCopyPoly.Subtract(tmpPoly);
//
// resets the tmp objects
//            tmpVal = 0;
//            tmpPoly=0;
//
//        }
//    }
//
//
	return retVal;


}


const Polynomial Polynomial::Derive()const{
    //cout << this->ToString() << endl;
    //makes the polynomial with highest degree of this-1
    Polynomial retVal(this->_degree-1);
    for(size_t i =1; i<this->_degree; ++i){
        retVal._coefficients[i-1] += i *this->_coefficients[i];
    }

    //cout << retVal.ToString() << endl;
	return retVal;
}




float Polynomial::Evaluate(float x)const{
    //cout << this->ToString() << endl;
    float answer = 0;

    for (size_t i = 0; i < this->_degree+1; ++i){
        answer += this->_coefficients[i] * pow(x,i);
    }
    //cout << value << endl;

	return answer;
}





float Polynomial::Integrate(float start, float end)const{
    //makes the polynomial with highest degree of this+1
    Polynomial retVal(this->_degree+1);

    //first value has to be zero for c
    retVal._coefficients[0]= 0.0;
    for(size_t i =0; i<this->_degree+1; ++i){
        // the f is for float, got help from lab assistant Parker, said it was needed to prevent integer division
        //dont need the x to power here
        retVal._coefficients[i+1] = (1.0f/(i+1))*this->_coefficients[i];
    }

    float defInt = retVal.Evaluate(end) - retVal.Evaluate(start);
    return defInt;
}





const Polynomial& Polynomial::operator=(const Polynomial& rhs){
	if (&rhs == this){
		return *this;
	}
	if (_degree != rhs._degree){
		if (_coefficients){
			delete[] _coefficients;
		}
		_degree = rhs._degree;
		_coefficients = new float[_degree + 1];
	}
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = rhs._coefficients[i];
	}
	return *this;
}



bool Polynomial::Equals(const Polynomial& rhs)const{
	if (_degree != rhs._degree){
		return false;
	}
	for (size_t i=0; i < _degree; i++){
		if (abs(_coefficients[i] - rhs._coefficients[i]) > 0.0001){
			return false;
		}
	}
	return true;
}
string Polynomial::ToString()const{
	stringstream ss;
	for (size_t i = _degree; i > 0; i--) {
		ss << showpos << fixed << setprecision(2) << _coefficients[i] << "x^" << i << " ";
	}
	ss << showpos << fixed << setprecision(2) << _coefficients[0];
	return ss.str();
}
ostream& Polynomial::Write(ostream& output)const{
	output << _degree << " ";
	for (size_t i = 0; i < _degree + 1; i++) {
		output << _coefficients[i] << " ";
	}
	return output;
}
istream& Polynomial::Read(istream& input){
	size_t degree;
	input >> degree;
	if (input.fail()){
		return input;
	}
	float* coefficients = new float[degree + 1];
	for (size_t i = 0; i < degree + 1; i++) {
		input >> coefficients[i];
		if (input.fail()){
			delete[] coefficients;
			return input;
		}
	}

	if (degree != _degree){
		if (_coefficients){
			delete[] _coefficients;
		}
		_degree = degree;
		_coefficients = coefficients;
	}else{
		for (size_t i = 0; i < _degree + 1; i++) {
			_coefficients[i] = coefficients[i];
		}
		delete[] coefficients;
	}
	return input;
}
