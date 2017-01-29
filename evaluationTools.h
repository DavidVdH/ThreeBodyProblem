/* 
 * File:   evaluationTools.h
 * Author: Andreas
 *
 * Created on 4 december 2014, 20:18
 */

#ifndef EVALUATIONTOOLS_H
#define	EVALUATIONTOOLS_H

#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "vec.h"
#include "rungeKutta.h"
#include "globals.h"

// calculates the energy for a set of physical coordinates and velocities
double calcEnergy(const vector<vec>& positions, const vector<vec>& velocities, const vector<double>& mass);
// calculates angular momentum of the system
double calcAngularMomentum(const vector<vec>& positions, const vector<vec>& velocities, const vector<double>& mass);
// calculates the distance between the two closest particles. If this distance is smaller than the treshold distance, closeEncounter is set true.
// The indices of the two particles envolved in the close encounter are stored in CEPos.
double calcMinDist(const vector<vec>& positions, bool& closeEncounter, vector<int>& CEPos);
// Write away all values if the time exceeds a certain value.
// This global variable (writeTime) is then updated with a fixed global varibale (writeStep)
void writeAway(const double& time, const vector<vec>& positions, const vector<vec>& velocities, const vector<double>& mass);
// Function that integrates the time evolution of the particles without regularisation.
void calcNoReg(vector<vec>& positions, vector<vec>& velocities, vector<double>& mass, double endTime);
// Function that integrates the time evolution of the particles with regularisation.
void calcReg(vector<vec>& positions, vector<vec>& velocities, vector<double>& mass, double endTime);
// function for testing purposes
double calcRegNoWrite(vector<vec>& positions, vector<vec>& velocities, vector<double>& mass, double endTime);
// Function that is called to start kind of a user interface
void usualProblem();
// Function that is called to evaluate the program and its parameters doing a lot of runs
void burrauTesting(double totalTime, double TAUBEGIN, double TAUEND, double TAUFACTOR, double HBEGIN, double HEND, double HFACTOR);

#endif	/* EVALUATIONTOOLS_H */

