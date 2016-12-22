#include "Sketch.h"

class SketchIO {
	private: 
		string filename;
	public:
		SketchIO( string filename);
		Sketch* read();
		void setFileName( string filename);
};

SketchIO::SketchIO( string filename) : filename(filename) {}

void SketchIO::setFileName( string filename) {
	(*this).filename = filename;
}

Sketch* SketchIO::read() {
	ifstream input(filename.c_str());
	double x,y,time;
	int strokeID,prevID;
	int strokeCount=0, pointCount=0;
	
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
