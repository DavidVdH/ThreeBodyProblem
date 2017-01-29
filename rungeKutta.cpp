#include "rungeKutta.h"

using namespace std;

// Acceleration for all particles in a given configuration (positions)

vector<vec> acceleration(const vector<vec>& positions, const vector<double>& mass) {

    vector<vec> acceleration;

    for (int i = 0; i < positions.size(); i++) {
        vec a_i(0, 0);

        for (int j = 0; j < positions.size(); j++) {
            if (j != i) {
                vec diff = positions[i] - positions[j];
                a_i += diff * ((-1) * G * mass[j] / (diff.norm2() * diff.norm()));
            }
        }

        acceleration.push_back(a_i);
    }

    return acceleration;
}


// Calculates the tidal force component F_ij (with i = self, j = other)

vec tidalForce(const vector<vec>& positions, const vector<double>& mass, int self, int other) {

    vec tidalForce(0, 0);

    for (int i = 0; i < positions.size(); i++) {
        if (i != self && i != other) {
            vec diff = positions[self] - positions[i];
            tidalForce += diff * ((-1) * G * mass[i] / (diff.norm2() * diff.norm()));
        }
    }

    return tidalForce;
}


// Calculates dw/d(tau)

vec wPrime(const vector<vec>& positions, const vec& u, const vector<double>& mass, const double& epsilon, const vector<int>& CEPos) {

    int obj1 = CEPos[0], obj2 = CEPos[1];

    vec F = tidalForce(positions, mass, obj1, obj2) - tidalForce(positions, mass, obj2, obj1);

    vec LF(u[0] * F[0] + u[1] * F[1], -u[1] * F[0] + u[0] * F[1]);

    return (u * epsilon + LF * u.norm2()) * 0.5;
}


// Calculates d(epsilon)/d(tau)

double ePrime(const vector<vec>& positions, const vec& u, const vec& w, const vector<double>& mass, const vector<int>& CEPos) {

    int obj1 = CEPos[0], obj2 = CEPos[1];

    vec F = tidalForce(positions, mass, obj1, obj2) - tidalForce(positions, mass, obj2, obj1);

    vec LF(u[0] * F[0] + u[1] * F[1], -u[1] * F[0] + u[0] * F[1]);

    return 2 * (w[0] * LF[0] + w[1] * LF[1]);
}


// Calculates the acceleration of CM and all particles that are not in close encounter

vector<vec> cmAcceleration(vector<vec> positions, vector<double> mass, const vector<int>& CEPos) {

    int obj1 = CEPos[0], obj2 = CEPos[1];

    vector<vec> accel;

    for (int i = 0; i < positions.size() - 1; i++) {
        vec a_i(0, 0);

        for (int j = 0; j < positions.size(); j++) {
            if (i != obj1 && i != obj2 && j != obj1 && j != obj2 && j != i) {
                vec diff = positions[i] - positions[j];
                a_i += diff * ((-1) * G * mass[j] / (diff.norm2() * diff.norm()));
            }
        }
        accel.push_back(a_i);
    }

    mass.pop_back();
    positions.pop_back();

    vec cmAccel = (tidalForce(positions, mass, obj1, obj2) * mass[obj1] + tidalForce(positions, mass, obj2, obj1) * mass[obj2]) * (1 / (mass[obj1] + mass[obj2]));

    accel.push_back(cmAccel);

    return accel;
}


// Transforms relative position and speed to regularised coordinates u and w

vector<vec> transform(vec r, vec v) {

    vector<vec> transformed;

    vec u, w;

    u.set(0, sign(r[1]) * sqrt((r.norm() + r[0]) / 2));
    u.set(1, sqrt((r.norm() - r[0]) / 2));

    w.set(0, 0.5 * (u[0] * v[0] + u[1] * v[1]));
    w.set(1, 0.5 * (-u[1] * v[0] + u[0] * v[1]));

    transformed.push_back(u);
    transformed.push_back(w);

    return transformed;
}


// Transforms regularised u and w to relative positions and speed

vector<vec> transformBack(vec u, vec w) {

    vector<vec> transformBack;

    vec r, v;

    r.set(0, u[0] * u[0] - u[1] * u[1]);
    r.set(1, 2 * u[0] * u[1]);

    v.set(0, 2 / u.norm2() * (u[0] * w[0] - u[1] * w[1]));
    v.set(1, 2 / u.norm2() * (u[1] * w[0] + u[0] * w[1]));

    transformBack.push_back(r);
    transformBack.push_back(v);

    return transformBack;
}


// calculates the updated positions (3 particles) from positions (3 particles, center of mass) and u

vector<vec> calcTempPos(vector<vec> positions, vector<double> mass, vec u, vector<int> CEPos, bool noCm) {

    int obj1 = CEPos[0], obj2 = CEPos[1];

    vec w(0, 0);
    vec r = transformBack(u, w)[0];
    positions[obj2] = positions.back() - r * (mass[obj1] / (mass[obj1] + mass[obj2]));
    positions[obj1] = positions[obj2] + r;

    if (noCm) {
        positions.pop_back();
    }

    return positions;
}


// calculates the updated velocities (3 particles) from velocities (3 particles, center of mass), u and w

vector<vec> calcTempVel(vector<vec> velocities, vector<double> mass, vec u, vec w, vector<int> CEPos) {

    int obj1 = CEPos[0], obj2 = CEPos[1];

    vec v = transformBack(u, w)[1];
    velocities[obj2] = velocities.back() - v * (mass[obj1] / (mass[obj1] + mass[obj2]));
    velocities[obj1] = velocities[obj2] + v;

    velocities.pop_back();

    return velocities;
}


// 4th order Runge Kutta with regularized coordinates. Positions, velocities and mass are the properties of the physical system. No
// regularisation or transformation to a Center of Mass system has taken place yet. u, w, epsilon and t are the 4 regularised coordinates.

void regRK4(vector<vec>& positions, vector<vec>& velocities, vec& u, vec& w, double& epsilon, double& t, const vector<double>& mass,
        const vector<int>& CEPos) {

    // Initialising variables

    vector<double> cmMass;
    vector<vec> cmPositions, cmVelocities;

    cmMass = mass;
    cmPositions = positions;
    cmVelocities = velocities;

    vector<vec> k1_r, k2_r, k3_r, k4_r, k1_v, k2_v, k3_v, k4_v;
    vec k1_u, k2_u, k3_u, k4_u, k1_w, k2_w, k3_w, k4_w;
    double k1_e, k2_e, k3_e, k4_e, k1_t, k2_t, k3_t, k4_t;

    vector<vec> tempPos;

    int obj1 = CEPos[0], obj2 = CEPos[1];

    //initialisation of cm coordinates

    cmMass.push_back((mass[obj1] + mass[obj2]));
    cmPositions.push_back((positions[obj1] * mass[obj1] + positions[obj2] * mass[obj2]) * (1 / (mass[obj1] + mass[obj2])));
    cmVelocities.push_back((velocities[obj1] * mass[obj1] + velocities[obj2] * mass[obj2]) * (1 / (mass[obj1] + mass[obj2])));

    cmVelocities[obj1].reset();
    cmVelocities[obj2].reset();

    //1st RK step

    k1_t = u.norm2() * tau; // tau is declared in globals for testing purposes

    k1_v = k1_t * cmAcceleration(cmPositions, cmMass, CEPos);
    k1_r = k1_t * cmVelocities;

    k1_w = wPrime(positions, u, mass, epsilon, CEPos) * tau;
    k1_u = w * tau;
    k1_e = ePrime(positions, u, w, mass, CEPos) * tau;

    // calculating of temporary positions to account for the change due to k1_r and k1_u
    tempPos = calcTempPos(cmPositions + 0.5 * k1_r, cmMass, u + k1_u * 0.5, CEPos, false);

    //2nd RK step

    k2_t = (u + k1_u * 0.5).norm2() * tau;

    k2_v = k2_t * cmAcceleration(tempPos, cmMass, CEPos);
    k2_r = k2_t * (cmVelocities + 0.5 * k1_v);

    // remove cm coordinate
    tempPos.pop_back();

    k2_w = wPrime(tempPos, u + k1_u * 0.5, mass, epsilon + k1_e * 0.5, CEPos) * tau;
    k2_u = (w + k1_w * 0.5) * tau;
    k2_e = ePrime(tempPos, u + k1_u * 0.5, w + k1_w * 0.5, mass, CEPos) * tau;


    tempPos = calcTempPos(cmPositions + 0.5 * k2_r, cmMass, u + k2_u * 0.5, CEPos, false);

    //3rd RK step

    k3_t = (u + k2_u * 0.5).norm2() * tau;

    k3_v = k3_t * cmAcceleration(tempPos, cmMass, CEPos);
    k3_r = k3_t * (cmVelocities + 0.5 * k2_v);

    tempPos.pop_back();

    k3_w = wPrime(tempPos, u + k2_u * 0.5, mass, epsilon + k2_e * 0.5, CEPos) * tau;
    k3_u = (w + k2_w * 0.5) * tau;
    k3_e = ePrime(tempPos, u + k2_u * 0.5, w + k2_w * 0.5, mass, CEPos) * tau;


    tempPos = calcTempPos(cmPositions + k3_r, cmMass, u + k3_u, CEPos, false);

    //4th RK step

    k4_t = (u + k3_u).norm2() * tau;

    k4_v = k4_t * cmAcceleration(tempPos, cmMass, CEPos);
    k4_r = k4_t * (cmVelocities + k3_v);

    tempPos.pop_back();

    k4_w = wPrime(tempPos, u + k3_u, mass, epsilon + k3_e, CEPos) * tau;
    k4_u = (w + k3_w) * tau;
    k4_e = ePrime(tempPos, u + k3_u, w + k3_w, mass, CEPos) * tau;

    //Calculating new quantities

    w += (k1_w + k2_w * 2 + k3_w * 2 + k4_w) * (1. / 6);
    u += (k1_u + k2_u * 2 + k3_u * 2 + k4_u) * (1. / 6);
    epsilon += (k1_e + k2_e * 2 + k3_e * 2 + k4_e) * (1. / 6);
    t += (k1_t + k2_t * 2 + k3_t * 2 + k4_t) * (1. / 6);

    cmVelocities = cmVelocities + (1. / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
    cmPositions = cmPositions + (1. / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);

    positions = calcTempPos(cmPositions, cmMass, u, CEPos, true);
    velocities = calcTempVel(cmVelocities, cmMass, u, w, CEPos);
}


// RK4 with physical coordinates. Is used when no close encounter occurs.

void RK4(vector<vec>& positions, vector<vec>& velocities, vector<double> mass, double h) {

    // defining coefficients ki_r and ki_v
    vector<vec> k1_r, k2_r, k3_r, k4_r, k1_v, k2_v, k3_v, k4_v;

    k1_v = h * acceleration(positions, mass);
    k1_r = h * velocities;

    k2_v = h * acceleration(positions + 0.5 * k1_r, mass);
    k2_r = h * (velocities + 0.5 * k1_v);

    k3_v = h * acceleration(positions + 0.5 * k2_r, mass);
    k3_r = h * (velocities + 0.5 * k2_v);

    k4_v = h * acceleration(positions + k3_r, mass);
    k4_r = h * (velocities + k3_v);

    velocities = velocities + (1. / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
    positions = positions + (1. / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);
}










