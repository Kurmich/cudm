#include "Sketch.h"

// a class to read a sketch from .sketch files
class SketchIO {
	private: 
		string filename;
	public:
		// constructor
		SketchIO( string filename);
		// a method loading a .sketch file into the memory
		Sketch* read();
		// a setter method for the sketch file
		void setFileName( string filename);
};

// constructor
SketchIO::SketchIO( string filename) : filename(filename) {}

// a setter method for the sketch file
void SketchIO::setFileName( string filename) {
	(*this).filename = filename;
}

// a method loading a .sketch file into the memory
Sketch* SketchIO::read() {
	ifstream input(filename.c_str());
	double x,y,time;
	int strokeID,prevID;
	int strokeCount=0, pointCount=0;
	
	if (input.is_open()) {
		if (input.good()) {
			// if the file is available
			// count the # of points and strokes first
			input >> x;
			input >> y;
			input >> strokeID;
			input >> time;
			
			prevID = strokeID;
			
			++strokeCount;
			++pointCount;
			
			while (true) {
				input >> x;
				
				if (!(input.good())) { 
					break; 
				}
				
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
	
	// then put these points and strokes into the sketch abstraction
	ifstream reinput(filename.c_str());
	Sketch* newSketch = new Sketch(pointCount, strokeCount);
	
	if (reinput.is_open()) {
		if (reinput.good()) {
			reinput >> x;
			reinput >> y;
			reinput >> strokeID;
			reinput >> time;
			
			prevID = strokeID;
			
			newSketch->openStroke();
			newSketch->addPoint(x,y);
			
			while (true) {
				reinput >> x;
				
				if (!(reinput.good())) break;
				
				reinput >> y;
				reinput >> strokeID;
				reinput >> time;
				
				if (prevID != strokeID) {
					newSketch->openStroke();
				}
				
				newSketch->addPoint(x,y);
				
				prevID = strokeID;
			}
		}
	}
	
	reinput.close();
	
	return newSketch;
}
