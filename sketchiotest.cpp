#include "SketchIO.h"

int main() {
	string filename = "ozan.sketch";
	SketchIO example(filename);
	SketchNew result = example.read();
	result.printContents();
	
	return 0;
}
