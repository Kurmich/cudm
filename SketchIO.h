#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

class SketchIO {
	private: 
		string filename;
	
	public:
		SketchIO( string filename);
		void read();
		void setFileName( string filename);
};

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

SketchIO::SketchIO( string filename) : filename(filename) {}

void SketchIO::setFileName( string filename) {
	(*this).filename = filename;
}

void SketchIO::read() {
	ifstream input(filename.c_str());
	double x,y,time;
	int strokeID,prevID;
	
	if (input.is_open()) {
		if (input.good()) {
			input >> x;
			input >> y;
			input >> strokeID;
			input >> time;
			
			prevID = strokeID;
			
			cout << "Stroke " << strokeID << endl;
			cout << x << " " << y << " " << time << endl;
			
			while (true) {
				input >> x;
				
				if (!(input.good())) break;
				
				input >> y;
				input >> strokeID;
				input >> time;
				
				if (prevID != strokeID) {
					cout << "Stroke " << strokeID << endl;
				}
				
				cout << x << " " << y << " " << time << endl;
				
				prevID = strokeID;
			}
		}
	}
	
	
	input.close();
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
	for (int i = 0; i < numStrokes; ++i) {
		for () { // left
		}
	}
} 
