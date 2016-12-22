#include "SketchIO.h"
#include <math.h>
#include <cmath>

#define PI 3.14159265

class FeatureExtractor {
	private:
		Sketch* sketch;
	public:
		FeatureExtractor(Sketch* sketch);
		void coords2angles(int *&strokeIndices, double *&angles, int &numAngles);
};

FeatureExtractor::FeatureExtractor(Sketch* sketch) : sketch(sketch) {}

void FeatureExtractor::coords2angles(int *&angleIndices, double *&angles, int &numAngles) {
	//Get stroke coordinates and indices
	int *sIndices = sketch->getStrokeIndices();
	double **sCoords = sketch->getCoords();
	
	//Initialize 
	angleIndices = new int[sketch->getNumStrokes()];
	numAngles = sketch->getNumPoints() - sketch->getNumStrokes();
	angles = new double[numAngles];
	
	int lastIndex;
	double angle,diffy,diffx;
	int curAngleIndex;
	
	for (int str = 0; str < sketch->getNumStrokes(); ++str) {
		//Get last index of current stroke
		if (str == sketch->getNumStrokes() - 1) {
			lastIndex = sketch->getNumPoints();
		}
		else {
			lastIndex = sIndices[str+1];
		}
		
		//Assign start index of angles for this stroke
		angleIndices[str] = sIndices[str]-str;
		//Starting index of angles to fill
		curAngleIndex = angleIndices[str];
		for (int pt = sIndices[str]+1; pt < lastIndex; ++pt) {
			//Get differences both in x and y directions
			diffy = sCoords[pt][1] - sCoords[pt-1][1];
			diffx = sCoords[pt][0] - sCoords[pt-1][0];
			
			//Compute angle
			angle = atan2(diffy,diffx);
			angle = fmod( (angle + 2*PI), (2*PI));
			angle *= 180.0/PI;

			//Assign current angle
			angles[curAngleIndex++] = angle;
		}
	}
}
