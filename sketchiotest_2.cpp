#include "FeatureExtractor.h"

int main() {
	string filename = "ozan.sketch";
	SketchIO example(filename);
	Sketch* result = example.read();
	result->printContents();
	FeatureExtractor fe(result);
	
	int *angleIndices;
	double *angles;
	int numAngles;
	double minX, minY, maxX, maxY;
	int numOfStrokes = result->getNumStrokes();
	fe.coords2angles(angleIndices,angles,numAngles, maxX, maxY, minX, minY);
	
	Sketch* rsm = result->resample(50);
	cout<<"resampled"<<endl;
	rsm->printContents();
	Sketch* trmd = result->transform(minX, minY, maxX, maxY);
	cout<<"transformed"<<endl;
	trmd->printContents();
	
	return 0;
}
