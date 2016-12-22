#include "FeatureExtractor.h"

int main() {
	string filename = "ozan.sketch";
	SketchIO example(filename);
	Sketch* result = example.read();
	result->printContents();
	FeatureExtractor fe(result);
	
	int *si;
	double *ang;
	int numAngles;
	fe.coords2angles(si,ang,numAngles);
	
	cout << "angles = " << endl;
	for (int i = 0; i < numAngles; ++i) {
		cout << ang[i] << endl;
	}
	
	double x,y;
	result->getCentroid(x,y);
	cout << "centroid x = " << x << ", y = " << y << endl;
	
	return 0;
}
