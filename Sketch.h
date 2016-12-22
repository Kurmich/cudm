#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

class Sketch {
	private:
		int numPoints,numStrokes;
		int ptAdded,strAdded;
		double **coords;
		int *strokeIndices;
	public:
		Sketch(int,int);
		~Sketch();
		void addPoint(double,double);
		void openStroke();
		void printContents();
		int getNumPoints();
		int getNumStrokes();
		int *getStrokeIndices();
		double **getCoords();
		void getCentroid(double &x, double &y);
};

Sketch::Sketch(int numPoints,int numStrokes) : numPoints(numPoints), numStrokes(numStrokes), ptAdded(0), strAdded(0) {
	coords = new double*[numPoints];
	for ( int i = 0; i < numPoints; ++i) {
		coords[i] = new double[2];
	}
	
	strokeIndices = new int[numStrokes];
}

Sketch::~Sketch() {
	for ( int i = 0; i < numPoints; ++i) {
		delete [] coords[i];
	}
	
	delete [] coords;
	
	delete [] strokeIndices;
}

void Sketch::getCentroid(double &x, double &y) {
	if ( numPoints > 0) {
		double xsum = 0;
		double ysum = 0;
		
		for (int i = 0; i < numPoints; ++i) {
			xsum += coords[i][0];
			ysum += coords[i][1];
		}
		
		xsum /= numPoints;
		ysum /= numPoints;
		
		x = xsum;
		y = ysum;
	}
	else {
		x = y = 0;
	}
}

void Sketch::addPoint(double x, double y) {
	if (ptAdded < numPoints) {
		coords[ptAdded][0] = x;
		coords[ptAdded++][1] = y;
	}
	else {
		cout << "ERROR: Sketch is full!" << endl;
	}
}

void Sketch::openStroke() {
	if (strAdded < numStrokes) {
		strokeIndices[strAdded++] = ptAdded;
	}
	else {
		cout << "ERROR: Sketch is full!" << endl;
	}
}

void Sketch::printContents() {
	int upperBound;
	
	for ( int i = 0; i < numStrokes; ++i) {
		cout << i << "=>";
		
		if ( i == numStrokes - 1) {
			upperBound = numPoints;
		}
		else {
			upperBound = strokeIndices[i+1];
		}
		
		for ( int j = strokeIndices[i]; j < upperBound; ++j) {
			cout << "(" << coords[j][0] << "," << coords[j][1] << ")";
			
			if ( j < upperBound - 1) {
				cout << "-";
			} 
		}
		
		cout << endl;
	}
} 

int Sketch::getNumPoints() {
	return numPoints;
}

int Sketch::getNumStrokes() {
	return numStrokes;
}

double** Sketch::getCoords() {
	return coords;
}

int* Sketch::getStrokeIndices() {
	return strokeIndices;
}
