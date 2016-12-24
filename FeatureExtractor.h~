#include "SketchIO.h"
#include <math.h>

#define PI 3.14159265

class FeatureExtractor {
	private:
		Sketch* sketch;
	public:
		FeatureExtractor(Sketch* sketch);
		void coords2angles(int *&strokeIndices, double *&angles, int &numAngles);
};

FeatureExtractor::FeatureExtractor(Sketch* sketch) : sketch(sketch) {}

void FeatureExtractor::coords2angles(int *&strokeIndices, double *&angles, int &numAngles) {
	int *sIndices = sketch->getStrokeIndices();
	double **sCoords = sketch->getCoords();
	
	strokeIndices = new int[sketch->getNumStrokes()];
	numAngles = sketch->getNumPoints() - sketch->getNumStrokes();
	angles = new double[numAngles];
	
	int lastIndex;
	double angle,diffy,diffx;
	
	for (int str = 0; str < sketch->getNumStrokes(); ++str) {
		if (str == sketch->getNumStrokes() - 1) {
			lastIndex = sketch->getNumPoints();
		}
		else {
			lastIndex = sIndices[str+1];
		}
		
		strokeIndices[str] = sIndices[str]-1-str;
		for (int pt = sIndices[str]+1; pt < lastIndex; ++pt) {
			diffy = sCoords[pt][1] - sCoords[pt-1][1];
			diffx = sCoords[pt][0] - sCoords[pt-1][0];
			
			angle = atan2(diffy,diffx);
			
			if ( angle < 0) {
				angle += 2*PI;
			}
			
			angle *= 180/PI;
			angles[pt-str-1] = angle;
		}
	}
}
