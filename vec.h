
#ifndef VEC_H
#define	VEC_H

/*
#ifdef	__cplusplus
extern "C" {
#endif




#ifdef	__cplusplus
}
#endif
*/


#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

class vec {
private:
    double x, y;

public:
    // Empty constructor
    vec();
    // Constructor with variables, for 2D z = 0
    vec(double first, double second);
    // Returns component of vector
    double operator[](int index) const;
    // Set x, y or z with index 0, 1, 2
    void set(int index, double value);
    // Overload += operator
    vec& operator+=(vec vecRHS);
    // Overload -= operator
    vec& operator-=(vec vecRHS);
    // Scalar product of a vector
    vec& operator*=(double scalar);
    // Calculates norm
    double norm() const;
    // Calculates square of norm
    double norm2() const;
    // Sets vector back to zero
    void reset();
    // Prints out values of vectors
    void print() const;
    // Stream to file
    void printFile(ofstream file) const;

    vec operator+(vec vector) const;
    vec operator-(vec vector) const;
    vec operator*(double scalar) const;
};

#endif	/* VEC_H */

