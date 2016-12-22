#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>

using namespace std;

class Sketch {
	private:
		int numPoints,numStrokes;
		int ptAdded,strAdded;
		double **coords;
		int *strokeIndices;
		void getCentroid(double &x, double &y);
		void getStd(double &x, double &y);
		double findMaxDistance();
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
		
		Sketch* resample(double rate);
		Sketch* normalized();
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

void Sketch::getStd(double &x, double &y) {
	if ( numPoints > 0) {
		double cx,cy;
		getCentroid(cx,cy);
		
		double xsqsum = 0,ysqsum = 0;
		for ( int i = 0; i < numPoints; ++i) {
			xsqsum += (coords[i][0] - cx)*(coords[i][0] - cx);
			ysqsum += (coords[i][1] - cy)*(coords[i][1] - cy);
		}
		
		xsqsum /= numPoints;
		ysqsum /= numPoints;
		
		x = sqrt(xsqsum);
		y = sqrt(ysqsum);
	}
	else {
		x = y = 0;
	}
}

double Sketch::findMaxDistance() {
	double x,y;
	getCentroid(x,y);
	
	double maxdist = -1;
	double curdist;
	for (int i = 0; i < numPoints; ++i) {
		curdist = (coords[i][1] - y)*(coords[i][1] - y) + (coords[i][0] - x)*(coords[i][0] - x);
		
		if (curdist > maxdist) {
			maxdist = curdist;
		}
	}
	
	return sqrt(maxdist);
}

Sketch* Sketch::normalized() {
	double cx,cy;
	double stdx,stdy;
	
	getCentroid(cx,cy);
	getStd(stdx,stdy);
	
	Sketch *newSketch = new Sketch(numPoints,numStrokes);
	
	int upperBound;
	
	for ( int i = 0; i < numStrokes; ++i) {
		newSketch->openStroke();
		
		if ( i == numStrokes - 1) {
			upperBound = numPoints;
		}
		else {
			upperBound = strokeIndices[i+1];
		}
		
		for ( int j = strokeIndices[i]; j < upperBound; ++j) {
			newSketch->addPoint((coords[j][0]-cx)/stdx,(coords[j][1]-cy)/stdy);
		}
	}
	
	return newSketch;
}

Sketch* Sketch::resample(double rate) {
	int newNumPoints = 0;
	double samplingInterval = findMaxDistance()*1.01 / rate;
	
	int upperBound;
	double distance;
	for (int i = 0; i < numStrokes; ++i) {
		if ( i == numStrokes - 1) {
			upperBound = numPoints;
		}
		else {
			upperBound = strokeIndices[i+1];
		}
		
		for (int j = strokeIndices[i]; j < upperBound - 1; ++j) {
			distance = sqrt((coords[j+1][1] - coords[j][1])*(coords[j+1][1] - coords[j][1]) + (coords[j+1][0] - coords[j][0])*(coords[j+1][0] - coords[j][0]));
			newNumPoints += (int) (floor(distance / samplingInterval))+1;
		}
		
		++newNumPoints;
	}
	
	Sketch* resampled = new Sketch(newNumPoints,numStrokes);
	
	int intPoints;
	double angleBtw,xdif,ydif;
	for (int i = 0; i < numStrokes; ++i) {
		resampled->openStroke();
		
		if ( i == numStrokes - 1) {
			upperBound = numPoints;
		}
		else {
			upperBound = strokeIndices[i+1];
		}
		
		for (int j = strokeIndices[i]; j < upperBound - 1; ++j) {
			resampled->addPoint(coords[j][0],coords[j][1]);
			xdif = coords[j+1][0] - coords[j][0];
			ydif = coords[j+1][1] - coords[j][1];
			distance = sqrt(ydif*ydif + xdif*xdif);
			
			angleBtw = atan2(ydif,xdif);
			intPoints = (int) (floor(distance / samplingInterval));
			
			for (int k = 1; k <= intPoints; ++k) {
				resampled->addPoint(coords[j][0]+samplingInterval*cos(angleBtw)*k,coords[j][1]+samplingInterval*sin(angleBtw)*k);
			}
		}
		
		resampled->addPoint(coords[upperBound-1][0],coords[upperBound-1][1]);
	}
	
	return resampled;
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
