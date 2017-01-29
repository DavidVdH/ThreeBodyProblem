#ifndef RUNGEKUTTA_H
#define	RUNGEKUTTA_H

#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "vec.h"
#include "vectorOperatoren.cpp"
#include "globals.h"

// Acceleration for all particles in a given configuration (positions)
vector<vec> acceleration(const vector<vec>& positions, const vector<double>& mass);
// Calculates the tidal force component F_ij (with i = self, j = other)
vec tidalForce(const vector<vec>& positions, const vector<double>& mass, int self, int other);
// Calculates dw/d(tau)
vec wPrime(const vector<vec>& positions, const vec& u, const vector<double>& mass, const double& epsilon, const vector<int>& CEPos);
// Calculates d(epsilon)/d(tau)
double ePrime(const vector<vec>& positions, const vec& u, const vec& w, const vector<double>& mass, const vector<int>& CEPos);
// Calculates the acceleration of CM and all particles that are not in close encounter
vector<vec> cmAcceleration(vector<vec> positions, vector<double> mass, const vector<int>& CEPos);
// Inline function to check the sign of a number
inline int sign(double number) {

    if (number >= 0) {
        return 1;
    } else {
        return -1;
    }
}
// Transforms relative position and speed to regularised coordinates u and w
vector<vec> transform(vec r, vec v);
// Transforms regularised u and w to relative positions and speed
vector<vec> transformBack(vec u, vec w);
// calculates the updated positions (3 particles) from positions (3 particles, center of mass) and u
vector<vec> calcTempPos(vector<vec> positions, vector<double> mass, vec u, vector<int> CEPos, bool noCm);
// calculates the updated velocities (3 particles) from velocities (3 particles, center of mass), u and w
vector<vec> calcTempVel(vector<vec> velocities, vector<double> mass, vec u, vec w, vector<int> CEPos);

// 4th order Runge Kutta with regularized coordinates. Positions, velocities and mass are the properties of the physical system. No
// regularisation or transformation to a Center of Mass system has taken place yet. u, w, epsilon and t are the 4 regularised coordinates.
void regRK4(vector<vec>& positions, vector<vec>& velocities, vec& u, vec& w, double& epsilon, double& t, const vector<double>& mass,
        const vector<int>& CEPos);

// RK4 with physical coordinates. Is used when no close encounter occurs.
void RK4(vector<vec>& positions, vector<vec>& velocities, vector<double> mass, double h);


#endif	/* RUNGEKUTTA_H */

