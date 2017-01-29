#include "vec.h"

// Empty constructor

vec::vec() {
    x = 0;
    y = 0;
}

// Constructor with variables, for 2D z = 0

vec::vec(double first, double second) {
    x = first;
    y = second;
}

// Returns component of vector

double vec::operator[](int index) const {
    if (index == 0) {
        return x;
    } else if (index == 1) {
        return y;
    } 
}

// Changes component of vector	

void vec::set(int index, double value) {
    if (index == 0) {
        x = value;
    } else if (index == 1) {
        y = value;
    } 
}

// Sum of two vectors (vecRHS = right vector)

vec& vec::operator+=(vec vecRHS) {
    x += vecRHS.x;
    y += vecRHS.y;
}

// Difference of two vectors (vecRHS = right vector)		

vec& vec::operator-=(vec vecRHS) {
    x -= vecRHS.x;
    y -= vecRHS.y;
}

// Scalar product of a vector

vec& vec::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
}

// Calculates norm

double vec::norm() const {
    return sqrt(x * x + y * y);
}

// Calculates square of norm

double vec::norm2() const {
    return x * x + y * y;
}

// Reset vector to zero

void vec::reset() {
    x = 0;
    y = 0;
}

// Prints out values of vectors

void vec::print() const {
    cout << '(' << x << ',' << y <<  ')' << '\n';
}

// Stream to file

void vec::printFile(ofstream file) const {
    file << x << '\t' << y;
}

vec vec::operator+(vec vector) const {
	return vec(x + vector.x, y + vector.y);
}

vec vec::operator-(vec vector) const {
    return vec(x - vector.x, y - vector.y);
}

vec vec::operator*(double scalar) const {
    return vec(x*scalar, y*scalar);
}





