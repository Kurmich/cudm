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

// a new representation? (remember that we will also transfer them to CUDA device memory)

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
