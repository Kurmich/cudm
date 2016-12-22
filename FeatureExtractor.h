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
		double* getMinAngleDistance(double *angles, double curAngle, double curAngle2, int numAngles);
		double truncate(double curDiff);
		double* pixelValues(double* minDist, int numAngles);
};

FeatureExtractor::FeatureExtractor(Sketch* sketch) : sketch(sketch) {}

void FeatureExtractor::coords2angles(int *&angleIndices, double *&angles, int &numOfAngles) {
	//Get stroke coordinates and indices
	int *sIndices = sketch->getStrokeIndices();
	double **sCoords = sketch->getCoords();
	
	//Initialize 
	angleIndices = new int[sketch->getNumStrokes()];
	numOfAngles = sketch->getNumPoints() - sketch->getNumStrokes();
	angles = new double[numOfAngles];
	
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

double FeatureExtractor::truncate(double curDiff)
{
	//truncates curDiff between 0 and 180
	if(curDiff >= 180)
	{
		curDiff = curDiff - 360;
	}
	else if(curDiff <= -180)
	{
		curDiff = curDiff + 360;
	}
	//truncates curDiff between 0 and 180
	return abs(curDiff);
}


double* FeatureExtractor::getMinAngleDistance(double* angles, double curAngle, double curAngle2, int numAngles)
{
	//Initialize array of angle distances
	double* diff = new double[numAngles];
	int curDiff, curDiff2;
	for(int i = 0; i < numAngles; ++i)
	{
		//get angle distances relative to current angles
		curDiff = angles[i] - curAngle;
		curDiff2 = angles[i] - curAngle2;
		//truncate angle distance between 0 and 180
		curDiff = truncate(curDiff);
		curDiff2 = truncate(curDiff2);
		//Assign minimum angle distance
		diff[i] = min(curDiff, curDiff2);
	}

	return diff;
}


double* FeatureExtractor::pixelValues(double* minDist, int numAngles)
{	/*
        Assign pixel values are calculated as a difference between stroke angle and the reference angle
        and vary linearly between 1.0(if the two are equal) and 0.0(if they differ by more than 45 degrees)
      */
	double* pixValues = new double[numAngles];
	double curPixel;
	double angleThreshold = 45;
	bool valid;
	for(int i = 0; i < numAngles; ++i)
	{
		valid = minDist[i] <= angleThreshold;
		if(valid)
		{
			curPixel = 1 - (minDist[i] / angleThreshold);
		}
		else
		{
			curPixel = 0;
		}
		pixValues[i] = curPixel;
	}

	return pixValues;
}

