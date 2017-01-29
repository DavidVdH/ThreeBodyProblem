#include "evaluationTools.h"

// values used through for writing and reference 
const double writeStep = 0.0005;
double writeTime = 0;
double energyRef;
double angularMomentumRef;
ofstream file;


// calculates the energy for a set of physical coordinates and velocities

double calcEnergy(const vector<vec>& positions, const vector<vec>& velocities, const vector<double>& mass) {

    double energy = 0;

    for (int i = 0; i < mass.size(); i++) {
        energy += 0.5 * mass[i] * velocities[i].norm2();

        for (int j = 0; j < i; j++) {
            vec diff = positions[i] - positions[j];
            energy -= G * mass[i] * mass[j] / diff.norm();
        }
    }

    return energy;
}


// calculates angular momentum of the system

double calcAngularMomentum(const vector<vec>& positions, const vector<vec>& velocities, const vector<double>& mass) {

    // calculate center of mass of the whole system, on a particular moment.
    vec centerOfMass(0, 0);
    double totalMass;

    for (int i = 0; i < positions.size(); i++) {
        centerOfMass += positions[i] * mass[i];
        totalMass += mass[i];
    }

    centerOfMass = centerOfMass * (1 / totalMass);

    double L = 0;
    for (int i = 0; i < positions.size(); i++) {
        L += (positions[i][0] - centerOfMass[0]) * velocities[i][1] * mass[i] - (positions[i][1] - centerOfMass[1]) * velocities[i][0] * mass[i];
    }

    return L;
}


// calculates the distance between the two closest particles. If this distance is smaller than the treshold distance, closeEncounter is set true.
// The indices of the two particles envolved in the close encounter are stored in CEPos.

double calcMinDist(const vector<vec>& positions, bool& closeEncounter, vector<int>& CEPos) {

    double minDist = (positions[1] - positions[0]).norm();
    CEPos[0] = 0;
    CEPos[1] = 1;

    for (int i = 2; i < positions.size(); i++) {
        for (int j = 0; j < i; j++) {

            double dist = (positions[i] - positions[j]).norm();

            if (dist < minDist) {
                minDist = dist;
                CEPos[0] = i;
                CEPos[1] = j;
            }
        }
    }

    if (minDist < CETresh) // CEtresh is value for close encounter
    {
        closeEncounter = true;
    } else {
        closeEncounter = false;
    }

    return minDist;
}


// Write away all values if the time exceeds a certain value.
// This global variable (writeTime) is then updated with a fixed global varibale (writeStep)

void writeAway(const double& time, const vector<vec>& positions, const vector<vec>& velocities, const vector<double>& mass) {

    // only write down at certain times
    if (time > writeTime) {
        file << time << '\t';

        for (int j = 0; j < positions.size(); j++) {
            file << positions[j][0] << '\t' << positions[j][1] << '\t';

        }

        file << abs((calcEnergy(positions, velocities, mass) - energyRef) / energyRef) << '\t';

        file << abs((calcAngularMomentum(positions, velocities, mass) - angularMomentumRef) / angularMomentumRef) << endl;

        writeTime += writeStep;
    }

}


// Function that integrates the time evolution of the particles without regularisation.

void calcNoReg(vector<vec>& positions, vector<vec>& velocities, vector<double>& mass, double endTime) {
    energyRef = calcEnergy(positions, velocities, mass);

    bool closeEncounter = false; // vars to check close encounter
    vector<int> CEPos(2, -1);

    double t = 0;

    while (t < endTime) {
        double h = 0.01 * calcMinDist(positions, closeEncounter, CEPos);
        t += h;
        file << t << '\t';

        for (int j = 0; j < positions.size(); j++) {
            file << positions[j][0] << '\t' << positions[j][1] << '\t';

        }

        file << abs((calcEnergy(positions, velocities, mass) - energyRef) / energyRef) << endl;

        RK4(positions, velocities, mass, h);
    }
}


// Function that integrates the time evolution of the particles with regularisation.

void calcReg(vector<vec>& positions, vector<vec>& velocities, vector<double>& mass, double endTime) {

    energyRef = calcEnergy(positions, velocities, mass);
    angularMomentumRef = calcAngularMomentum(positions, velocities, mass);

    bool closeEncounter = false; // vars to check close encounter
    vector<int> CEPos(2, -1);

    double t = 0;
    double h = 0.01;

    while (t < endTime) {

        if (closeEncounter) {

            int obj1 = CEPos[0], obj2 = CEPos[1];

            vec r = positions[obj1] - positions[obj2], v = velocities[obj1] - velocities[obj2];

            vector<vec> regCoord = transform(r, v);
            vec u = regCoord[0], w = regCoord[1];

            double epsilon = (2 * w.norm2() - G * (mass[obj1] + mass[obj2])) / u.norm2();
            double time = 0;

            while (closeEncounter) {
                writeAway(t, positions, velocities, mass);

                // positions, velocities, u, w, epsilon & time updated in regRK4
                regRK4(positions, velocities, u, w, epsilon, time, mass, CEPos);

                h += time;
                time = 0;

                // calcMinDist checks if two particles are still in closeEncounter. Variable minDist isn't used anywhere.
                double minDist = calcMinDist(positions, closeEncounter, CEPos);

                // check if the same particles are closest (if for example the three particles are in close encounter)
                if (CEPos[0] != obj1 || CEPos[1] != obj2) {

                    // initialize again all functions. We are still in close encounter, so we stay in the while loop
                    int obj1 = CEPos[0], obj2 = CEPos[1];

                    vec r = positions[obj1] - positions[obj2], v = velocities[obj1] - velocities[obj2];

                    vector<vec> regCoord = transform(r, v);
                    vec u = regCoord[0], w = regCoord[1];

                    epsilon = (2 * w.norm2() - G * (mass[obj1] + mass[obj2])) / u.norm2();
                }

            }

        } else {
            writeAway(t, positions, velocities, mass);
            h = timeStep * pow(calcMinDist(positions, closeEncounter, CEPos), 3. / 2);
            RK4(positions, velocities, mass, h);
        }

        t += h;
    }

}


// for testing purposes. If we only want to know the final energy, there is no need to write away all values. This saves a lot of time.

double calcRegNoWrite(vector<vec>& positions, vector<vec>& velocities, vector<double>& mass, double endTime) {

    energyRef = calcEnergy(positions, velocities, mass);
    angularMomentumRef = calcAngularMomentum(positions, velocities, mass);

    bool closeEncounter = false; // vars to check close encounter
    vector<int> CEPos(2, -1);

    double t = 0;
    double h = 0.01;

    while (t < endTime) {

        if (closeEncounter) {

            int obj1 = CEPos[0], obj2 = CEPos[1];

            vec r = positions[obj1] - positions[obj2], v = velocities[obj1] - velocities[obj2];

            vector<vec> regCoord = transform(r, v);
            vec u = regCoord[0], w = regCoord[1];

            double epsilon = (2 * w.norm2() - G * (mass[obj1] + mass[obj2])) / u.norm2();
            double time = 0;

            while (closeEncounter) {

                // positions, velocities, u, w, epsilon & time updated in regRK4
                regRK4(positions, velocities, u, w, epsilon, time, mass, CEPos);

                h += time;
                time = 0;

                // calcMinDist checks if two particles are still in closeEncounter. Variable minDist isn't used anywhere.
                double minDist = calcMinDist(positions, closeEncounter, CEPos);

                // check if the same particles are closest (if for example the three particles are in close encounter)
                if (CEPos[0] != obj1 || CEPos[1] != obj2) {

                    // initialize again all functions. We are still in close encounter, so we stay in the while loop
                    int obj1 = CEPos[0], obj2 = CEPos[1];

                    vec r = positions[obj1] - positions[obj2], v = velocities[obj1] - velocities[obj2];

                    vector<vec> regCoord = transform(r, v);
                    vec u = regCoord[0], w = regCoord[1];

                    epsilon = (2 * w.norm2() - G * (mass[obj1] + mass[obj2])) / u.norm2();
                }

            }

        } else {
            h = timeStep * calcMinDist(positions, closeEncounter, CEPos);
            RK4(positions, velocities, mass, h);
        }

        t += h;
    }

    return abs((calcEnergy(positions, velocities, mass) - energyRef) / energyRef);

}

void usualProblem() {
    vector<vec> positions, velocities;
    vector<double> mass;


    char burrau;
    double time;

    cout << "Would you like to solve Burrau's problem? (ONLY y/n)" << endl;
    cin >> burrau;
    if (burrau == 'y') {

        positions.push_back(vec(1, -1));
        positions.push_back(vec(-2, -1));
        positions.push_back(vec(1, 3));


        velocities.push_back(vec(0, 0));
        velocities.push_back(vec(0, 0));
        velocities.push_back(vec(0, 0));

        mass.push_back(5);
        mass.push_back(4);
        mass.push_back(3);
    } else if (burrau == 'n') {
        int numbPart = 0;
        double x, y, vx, vy;
        double mass_i = 0;

        cout << "How many particles would you like to simulate with?" << endl;
        cin >> numbPart;

        for (int i = 0; i < numbPart; i++) {

            cout << "Insert coordinates of particle " << i + 1 << ", each component on a different line." << endl;
            cin >> x;
            cin >> y;
            positions.push_back(vec(x, y));

            cout << "Insert velocity of particle " << i + 1 << ", each component on a different line." << endl;
            cin >> vx;
            cin >> vy;
            velocities.push_back(vec(vx, vy));

            cout << "Insert mass of particle " << i + 1 << " (Be positive!)" << endl;
            cin >> mass_i;
            mass.push_back(mass_i);
        }
    }

    cout << "How long would you like to simulate? (in years, 15 is a good value)" << endl;
    cin >> time;

    cout << "Program is running..." << endl;


    file.open("sol.dat");
    file << setprecision(8);

    calcReg(positions, velocities, mass, time);

    file.close();
}

void burrauTesting(double totalTime, double TAUBEGIN, double TAUEND, double TAUFACTOR,
        double HBEGIN, double HEND, double HFACTOR) {

    // initialisation of Burrau problem

    vector<vec> positions, velocities;
    vector<double> mass;

    ofstream energyfile;

    energyfile.open("energy.dat");
    energyfile << setprecision(8);

    double energyError = 0;
    double tauValue = TAUBEGIN;
    double htime = HBEGIN;
    double calcTime = 0;

    energyfile << "#Tau \t htime \t time to calculate" << endl;
    energyfile << "#CETresh was 0.006, pow(mindist, 1), integrated system for 10 years" << endl;
    cout << "CETreshold \t energy \t htime \t time to calculate" << endl;

    cout << "run begonnen" << endl;

    while (htime <= HEND) {
        // changing variable in globals.cpp
        timeStep = htime;

        while (tauValue <= TAUEND) {
            // changing variable in globals.cpp
            tau = tauValue;

            // reinitializing positions, velocities, mass.
            
            positions.push_back(vec(1, -1));
            positions.push_back(vec(-2, -1));
            positions.push_back(vec(1, 3));

            velocities.push_back(vec(0, 0));
            velocities.push_back(vec(0, 0));
            velocities.push_back(vec(0, 0));

            mass.push_back(5);
            mass.push_back(4);
            mass.push_back(3);
            

            calcTime = clock();
            energyError = calcRegNoWrite(positions, velocities, mass, totalTime);
            calcTime = clock() - calcTime;
            
            energyfile << tauValue << '\t' << htime << '\t' << energyError << '\t' << ((float) calcTime) / CLOCKS_PER_SEC << endl;
            
            cout << tauValue << '\t' << htime << '\t' << ((float) calcTime) / CLOCKS_PER_SEC << endl;            
            tauValue *= TAUFACTOR;


            // clearing old positions
            positions.clear();
            velocities.clear();
            mass.clear();
            

        }
        
        tauValue = TAUBEGIN; 

        htime *= HFACTOR;
        
    }

    energyfile.close();

}
