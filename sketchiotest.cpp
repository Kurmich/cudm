#include "SketchIO.h"

int main() {
	string filename = "ozan.sketch";
	SketchIO example(filename);
	example.read();
	
	return 0;
}
