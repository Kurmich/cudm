#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

class SketchNew {
	private:
		int numPoints,numStrokes;
		int ptAdded,strAdded;
		int **coords;
		int *strokeIndices;
	public:
		SketchNew(int,int);
		~SketchNew();
		void addPoint(int,int);
		void openStroke();
		void printContents();
};

class SketchIO {
	private: 
		string filename;
	
	public:
		SketchIO( string filename);
		SketchNew read();
		void setFileName( string filename);
};

SketchIO::SketchIO( string filename) : filename(filename) {}

void SketchIO::setFileName( string filename) {
	(*this).filename = filename;
}

SketchNew SketchIO::read() {
	ifstream input(filename.c_str());
	double x,y,time;
	int strokeID,prevID;
	int strokeCount, pointCount;
	
	if (input.is_open()) {
		if (input.good()) {
			input >> x;
			input >> y;
			input >> strokeID;
			input >> time;
			
			prevID = strokeID;
			
			++strokeCount;
			++pointCount;
			
			while (true) {
				input >> x;
				
				if (!(input.good())) break;
				
				input >> y;
				input >> strokeID;
				input >> time;
				
				if (prevID != strokeID) {
					++strokeCount;
				}
				
				++pointCount;
				
				prevID = strokeID;
			}
		}
	}
	
	input.close();
	
	ifstream reinput(filename.c_str());
	SketchNew newSketch(pointCount, strokeCount);
	
	if (reinput.is_open()) {
		if (reinput.good()) {
			reinput >> x;
			reinput >> y;
			reinput >> strokeID;
			reinput >> time;
			
			prevID = strokeID;
			
			newSketch.openStroke();
			newSketch.addPoint(x,y);
			
			while (true) {
				reinput >> x;
				
				if (!(reinput.good())) break;
				
				reinput >> y;
				reinput >> strokeID;
				reinput >> time;
				
				if (prevID != strokeID) {
					newSketch.openStroke();
				}
				
				newSketch.addPoint(x,y);
				
				prevID = strokeID;
			}
		}
	}
	
	reinput.close();
	
	return newSketch;
}

SketchNew::SketchNew(int numPoints,int numStrokes) : numPoints(numPoints), numStrokes(numStrokes), ptAdded(0), strAdded(0) {
	coords = new int*[numPoints];
	for ( int i = 0; i < numPoints; ++i) {
		coords[i] = new int[2];
	}
	
	strokeIndices = new int[numStrokes];
}

SketchNew::~SketchNew() {
	for ( int i = 0; i < numPoints; ++i) {
		delete [] coords[i];
	}
	
	delete [] coords;
	
	delete [] strokeIndices;
}

void SketchNew::addPoint(int x, int y) {
	if (ptAdded < numPoints) {
		coords[ptAdded][0] = x;
		coords[ptAdded++][1] = y;
	}
	else {
		cout << "ERROR: Sketch is full!" << endl;
	}
}

void SketchNew::openStroke() {
	if (strAdded < numStrokes) {
		strokeIndices[strAdded++] = ptAdded;
	}
	else {
		cout << "ERROR: Sketch is full!" << endl;
	}
}

void SketchNew::printContents() {
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
