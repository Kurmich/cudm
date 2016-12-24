#include "FeatureExtractor.h"

int main() {
	string filename = "ozan.sketch";
	SketchIO example(filename);
	Sketch* result = example.read();
	result->printContents();
	FeatureExtractor fe(result);
	
	fe.extract();
	/*
	cout<<"maxX = "<<maxX<<" maxY "<<maxY<<" minX = "<<minX<<" minY = "<<minY<<endl;
	
	cout << "angles = " << endl;
	for (int i = 0; i < numAngles; ++i) {
		cout << angles[i] << endl;
	}

	cout<<"angle indices"<<endl;
	for(int i = 0; i < numOfStrokes; ++i)
	{
		cout<<angleIndices[i]<<endl;
	}

	double* diff = fe.getMinAngleDistance(angles, 45, 45, numAngles);

	cout << "diffs = " << endl;
	for (int i = 0; i < numAngles; ++i) {
		cout << diff[i] << endl;
	}

	double* pixels = fe.pixelValues(diff, numAngles);

	cout << "pixels = " << endl;
	for (int i = 0; i < numAngles; ++i) {
		cout << pixels[i] << endl;
	}
	
	double x,y;
	result->getCentroid(x,y);
	cout << "centroid x = " << x << ", y = " << y << endl;
	cout << "Maxdist = " << result->findMaxDistance() << endl;
	*/
	return 0;
}
