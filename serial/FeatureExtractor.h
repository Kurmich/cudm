#include "SketchIO.h"
#include <math.h>
#include <cmath>
#include <vector>

#define PI 3.14159265

class FeatureExtractor {
	private:
		Sketch* sketch;
	public:
		// constructor
		FeatureExtractor(Sketch* sketch);
		// a method that finds the angle of a line relative to the x axis
		// between every two consecutive sketch point in every stroke
		void coords2angles(int *&strokeIndices, double *&angles, int &numAngles, double &maxX, double &maxY, double &minX, double &minY);
		// given the reference angle, construct the difference between each angle in angles list
		// and the ref. angle
		double* getMinAngleDistance(double *angles, double curAngle, double curAngle2, int numAngles);
		// a helper method in getMinAngleDistance
		double truncate(double curDiff);
		// compute the pixels in 24x24 grid using this method
		double* pixelValues(double *angles, double curAngle, double curAngle2, int numAngles);
		// extractor method
		double* extract();
		// this method sets the sketch to be extracted to the sketch given in the parameter
		void setSketch(Sketch* newSketch);
		// range constructor used in drawBresenham
		vector<int> arange(int a, int b, int step);
		// cumulative computation in a specific operation
		// (used in drawBresenham)
		vector<int> cum(vector<int> q,  int x, char op);
		// bresenham algorithm used in construction of feature images
		void drawBresenham(  double x1,  double y1, double x2, double y2, double* pixels, int angleIndex, double** &image );
		// a helper method to construct a 2d array
		double** init2Darray(int size);
		// extraction of a single feature image
		double** extractFeatureImage(double** gfilter, double* pixels, int* angleIndices, int numAngles, Sketch* transformed,  int hsize, int gridSize, bool endpt);
		// mapper from the sketch points to a single feature image
		void pointsToImage(double* pixels, int* angleIndices, double** sCoords, int gridSize, int strokeStart, int strokeEnd, int angleStart, int angleEnd, double** &image);
		// downsampler for 24x24 images
		double** downsample(double** image, double gridSize);
		// gaussian filter for 24x24 feature images
		double** gaussianFilter(int hsize, double sigma);
		// image smoothing before downsampling of a feature image of size 24x24
		double** smoothim(double **image, double** fgauss, int hsize, int gridSize);
};

// constructor
FeatureExtractor::FeatureExtractor(Sketch* sketch) : sketch(sketch) {}

// this method sets the sketch to be extracted to the sketch given in the parameter
void FeatureExtractor::setSketch(Sketch* newSketch)
{
	sketch = newSketch;
}

// extractor method
double* FeatureExtractor::extract()
{
	// init. the intrinsic parameters
	double resampleInterval = 50.0;
    double sigma = 10.0;
    int hsize = 4.0;
    int gridSize = 12;
	int fimgSize = gridSize*gridSize;
    int idmFeatureSize = fimgSize*5;
    double* idmFeature = new double[idmFeatureSize];
	// resample the sketch
    Sketch* resampled = sketch->resample(resampleInterval);
	// normalize the sketch
    Sketch* normalized = resampled->normalized();
    int *angleIndices;
	double *angles;
	int numAngles;
	double minX, minY, maxX, maxY;
	int numOfStrokes = normalized->getNumStrokes();
	// we'll use the normalized sketch in feature extraction
	setSketch(normalized);
	// extract the line angles between every two consecutive sketch point in every stroke
	coords2angles(angleIndices,angles,numAngles, maxX, maxY, minX, minY);
	// using the extreme coordinates, compute the transformed sketch
	Sketch* transformed = normalized->transform(minX, minY, maxX, maxY);
	double** gfilter = gaussianFilter(hsize, sigma);

	// compute 5 different feature images
	int curAngle = 0;
	int curAngle2 = (curAngle + 180)%360;
	double* pixels1 = pixelValues(angles, curAngle, curAngle2 , numAngles);
	double** featImage1 = extractFeatureImage(gfilter,pixels1, angleIndices, numAngles, transformed,  hsize,  gridSize, false );


	curAngle = 45;
	curAngle2 = (curAngle + 180)%360;
	double* pixels2 = pixelValues(angles, curAngle, curAngle2, numAngles);
	double** featImage2 = extractFeatureImage(gfilter,pixels2, angleIndices, numAngles, transformed,  hsize,  gridSize, false );


	curAngle = 90;
	curAngle2 = (curAngle + 180)%360;
	double* pixels3 = pixelValues(angles, curAngle, curAngle2, numAngles);
	double** featImage3 = extractFeatureImage(gfilter,pixels3, angleIndices, numAngles, transformed,  hsize,  gridSize, false );
	
	curAngle = 135;
	curAngle2 = (curAngle + 180)%360;
	double* pixels4 = pixelValues(angles, curAngle, curAngle2, numAngles);
	double** featImage4 = extractFeatureImage(gfilter,pixels4, angleIndices, numAngles, transformed,  hsize,  gridSize, false );

	// feature image for endpoints
	double** featImage5 = extractFeatureImage(gfilter,pixels1, angleIndices, numAngles, transformed,  hsize,  gridSize, true );

	// merge & linearize the feature images in a 720-element array
	for (int i = 0; i < gridSize; ++i) {
		for (int j = 0; j < gridSize; ++j) {
			idmFeature[i*gridSize+j] = featImage1[i][j];
			idmFeature[fimgSize+i*gridSize+j] = featImage2[i][j];
			idmFeature[2*fimgSize+i*gridSize+j] = featImage3[i][j];
			idmFeature[3*fimgSize+i*gridSize+j] = featImage4[i][j];
			idmFeature[4*fimgSize+i*gridSize+j] = featImage5[i][j];
		}
		delete [] featImage1[i];
		delete [] featImage2[i];
		delete [] featImage3[i];
		delete [] featImage4[i];
		delete [] featImage5[i];
	}
	
	for (int i = 0; i <= hsize; ++i) {
		delete [] gfilter[i];
	}
	
	delete resampled;
	delete normalized;
	delete transformed;
	
	delete [] gfilter;
	delete [] pixels1;
	delete [] pixels2;
	delete [] pixels3;
	delete [] pixels4;
	delete [] angleIndices;
	delete [] featImage1;
	delete [] featImage2;
	delete [] featImage3;
	delete [] featImage4;
	delete [] featImage5;
	delete [] angles;
	
    return idmFeature;
}

// extraction of a single feature image
double** FeatureExtractor::extractFeatureImage(double** gfilter, double* pixels, int* angleIndices, int numAngles, Sketch* transformed, int hsize, int gridSize, bool endpt )
{
	// init. the feat. image on the memory
	// and the variables
	double** featim = init2Darray(2*gridSize);;
	int* sIndices = transformed->getStrokeIndices();
	double** sCoords = transformed->getCoords();
	int numOfPoints = transformed->getNumPoints();  
	int numOfStrokes = transformed->getNumStrokes();
	if(!endpt)
	{
		// if we will extract the feature image based on a reference angle
		int strokeStart, strokeEnd, angleStart, angleEnd;
		for(int i = 0; i < numOfStrokes; ++i)
		{
			strokeStart = sIndices[i];
			angleStart = angleIndices[i];
			
			if(i < numOfStrokes - 1 )
			{
				strokeEnd = sIndices[i+1];
				angleEnd = angleIndices[i+1]; 
			}
			else
			{
				angleEnd = numAngles;
				strokeEnd = numOfPoints;
			}
			
			// construct the feature image
			pointsToImage(pixels, angleIndices, sCoords, gridSize, strokeStart, strokeEnd, angleStart, angleEnd, featim);
		}
	}
	else
	{
		// if we will extract the feature image based on the endpoints of all strokes
		int strokeStart, strokeEnd;
		for(int i = 0; i < numOfStrokes; ++i)
		{
			strokeStart = sIndices[i];
			if(i < numOfStrokes - 1 )
			{
				strokeEnd = sIndices[i+1] - 1; 
			}
			else
			{
				strokeEnd = numOfPoints - 1;
			}
			// mark the corresponding pixel of every endpoint of every stroke
			featim[ (int)sCoords[strokeStart][1] ][ (int)sCoords[strokeStart][0] ] = 1;
			featim[ (int)sCoords[strokeEnd][1] ][ (int)sCoords[strokeEnd][0] ] = 1;
		}
	}
	
	// smooth & downsample the feature image
	double** smoothed = smoothim(featim, gfilter, hsize, gridSize);
	double** downsampled = downsample(smoothed, gridSize);
	
	for (int i = 0; i < 2*gridSize; ++i) {
		delete [] smoothed[i];
	}
	delete [] smoothed;
	
	return downsampled; 
}

// a function constructing a gaussian filter given the filter size and sigma
double** FeatureExtractor::gaussianFilter(int hsize, double sigma)
{
	double** gfilter = init2Darray(hsize+1);
	double sum = 0;
	double r, s = 2.0 * sigma * sigma;
	int d = hsize/2;
	
	// construct the filter in every entry (i,j)
	// using bi-Gaussian distribution
	for (int x = -d; x <= d; x++) // Loop to generate 5x5 kernel
    {
        for(int y = -d; y <= d; y++)
        {
            r = sqrt(x*x + y*y);
            gfilter[x + d][y + d] = (exp(-(r*r)/s))/(PI * s);
            sum += gfilter[x + d][y + d];
        }
    }
	
     for(int i = 0; i <= hsize; ++i) // Loop to normalize the kernel
        for(int j = 0; j <= hsize; ++j)
            gfilter[i][j] /= sum;
            
    return gfilter;
}

// image smoothing using a given smoothing kernel (in our case, the filter is Gaussian)
double** FeatureExtractor::smoothim(double **image, double** fgauss, int hsize, int gridSize)
{
	double** result = init2Darray(2*gridSize);
	int ax, ay;
	double sum = 0, maximum = 0;
	for(int i = 0; i < 2*gridSize; ++i)
	{
		for(int j = 0; j < 2*gridSize; ++j)
		{
			sum = 0;
			for(int fi = 0; fi < hsize + 1; ++fi)
			{
				for(int fj = 0; fj < hsize + 1; ++fj)
				{
					ax = i + (fi - hsize/2);
                    ay = j + (fj - hsize/2);
                    if(ax >= 0 && ay >= 0 && ax < 2*gridSize && ay < 2*gridSize)
                    {
                        sum = sum + image[ax][ay] * fgauss[fi][fj];
                    }  
				}
			}
			result[i][j] = sum;
			maximum = max(maximum, sum);
		}
	}
	if(maximum != 0)
	{
		// normalization
		for(int i = 0; i < 2*gridSize; ++i)
		{
			for(int j = 0; j < 2*gridSize; ++j)
			{
				result[i][j] /= maximum;
			}
		}
	}

	
	return result;
}

// downsampling of an image through maxpooling
double** FeatureExtractor::downsample(double** image, double gridSize)
{
	double** downsampled = init2Darray(gridSize);
	double curMax;
	for(int i = 0; i < 2*gridSize; i += 2)
	{
		for(int j = 0; j < 2*gridSize; j += 2)
		{
			curMax = max(image[i][j], image[i+1][j]);
			curMax = max(curMax, image[i][j+1]);
			curMax = max(curMax, image[i+1][j+1]);
			downsampled[i/2][j/2] = curMax;
		}
	}

	return downsampled;
}

// mapping from sketch points to feature images
void FeatureExtractor::pointsToImage(double* pixels, int* angleIndices, double** sCoords, int gridSize, int strokeStart, int strokeEnd, int angleStart, int angleEnd, double** &image)
{
	int numOfAngles = angleEnd - angleStart;
	
	if(numOfAngles == 0)
	{
		return;
	}
	int angleIndex, pointIndex;
	double x1, x2, y1, y2;
	for(int i = 0; i < numOfAngles; ++i)
	{
		angleIndex = angleStart + i;
		pointIndex = strokeStart + i;
		x1 = sCoords[pointIndex][0];
		x2 = sCoords[pointIndex+1][0];
		y1 = sCoords[pointIndex][1];
		y2 = sCoords[pointIndex+1][1];
		
		if(pixels[angleIndex] > 0)
		{
			drawBresenham(  x1,  y1, x2, y2, pixels, angleIndex, image );
		}
	}
	
}

// 2d array initialization
double** FeatureExtractor::init2Darray(int size)
{
	double** array = new double*[size];
	for(int i =0; i < size; ++i)
	{
		array[i] = new double[size];
	}

	for(int i = 0; i < size; ++i)
	{
		for(int j = 0; j < size; ++j)
		{
			array[i][j] = 0;
		}
	}
	return array;
}

// the same as range function in python
vector<int> FeatureExtractor::arange(int a, int b, int step)
{
  vector<int> ans;
  if(step == 1)
  {
    for(int i = a; i <= b; ++i)
    {
      ans.push_back(i);
    }
  }
  else
  {
    for(int i = a; i >= b; --i)
    {
      ans.push_back(i);
    }
  }
  return ans;
}

// taking cumulative in an array based on a given operation
vector<int> FeatureExtractor::cum(vector<int> q,  int x, char op)
{
  vector<int> ans;
  if(op == '+')
  {
    for(int i = 0; i < q.size(); ++i)
      {
        ans.push_back(q[i] + x );
      }
  }
  else
  {
    for(int i = 0; i < q.size(); ++i)
      {
        ans.push_back( x - q[i] );
      }
  }
  return ans;
}

// bresenham drawing algorithm
void FeatureExtractor::drawBresenham(  double x1,  double y1, double x2, double y2, double* pixels, int angleIndex, double** &image )
{
  x1 = round(x1);
  x2 = round(x2);
  y1 = round(y1);
  y2 = round(y2);
  double dx = abs(x2 - x1);
  double dy = abs(y2 - y1);
  bool steep = abs(dy) > abs(dx);
  vector<int> q;

  if(steep)
  {
    int tmp = dx;
    dx = dy;
    dy = tmp;
  }
  if(dy == 0)
  {
    for(int  i = 0; i <= dx; ++i)
    {
      q.push_back(0);
    }
  }
  else
  {
    q.push_back(0);
    vector<int> arr;
    int cur;
    int cumSum = 0;
    for(int i = floor(dx/2); i >= floor(dx/2) - dy*dx; i -= dy)
    {
      cur = (i + (int) dx)%((int) dx);
      while(cur < 0)
      {
        cur = (cur + (int) dx)%((int) dx);
      }
      arr.push_back(cur);
    }
    for(int i = 1; i < arr.size(); ++i)
      {
        cur = arr[i] - arr[i-1];
        cur = cur >=0?1:0;
        cumSum += cur;
        q.push_back(cumSum);
      }
  }
  vector<int> y;
  vector<int> x;
  if(steep)
  {
    if(y1 <= y2)
    {
      y = arange(y1, y2, 1);
    }
    else
    {
      y = arange(y1, y2, -1);
    }

    if(x1 <= x2)
    {
      x = cum(q, x1, '+');
    }
    else
    {
      x = cum(q, x1, '-');
    }
  }
  else
  {
    if(x1 <= x2)
    {
      x = arange(x1, x2, 1);
    }
    else
    {
      x = arange(x1, x2, -1);
    }

    if(y1 <= y2)
    {
      y = cum(q, y1, '+');
    }
    else
    {
      y = cum(q, y1, '-');
    }
  }
  
  for(int i = 0; i < x.size(); ++i)
  {
        if(image[ y[i] ][ x[i] ] < pixels[angleIndex])
        {
        	image[ y[i] ][ x[i] ] = pixels[angleIndex];
        }
  }
}

// output the angles to the x-axis between every two consecutive point
void FeatureExtractor::coords2angles(int *&angleIndices, double *&angles, int &numOfAngles, double &maxX, double &maxY, double &minX, double &minY) {
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
	minX = sCoords[0][0];
	maxX = sCoords[0][0];
	minY = sCoords[0][1];
	maxY = sCoords[0][1];
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
		//update maxs and mins
		minX = min(minX, sCoords[sIndices[str]][0]);
		minY = min(minY, sCoords[sIndices[str]][1]);
		maxX = max(maxX, sCoords[sIndices[str]][0]);
		maxY = max(maxY, sCoords[sIndices[str]][1]);
		for (int pt = sIndices[str]+1; pt < lastIndex; ++pt) {
			//Get differences both in x and y directions
			diffy = sCoords[pt][1] - sCoords[pt-1][1];
			diffx = sCoords[pt][0] - sCoords[pt-1][0];
			minX = min(minX, sCoords[pt][0]);
			minY = min(minY, sCoords[pt][1]);
			maxX = max(maxX, sCoords[pt][0]);
			maxY = max(maxY, sCoords[pt][1]);
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
	double curDiff, curDiff2;
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


double* FeatureExtractor::pixelValues(double* angles, double curAngle, double curAngle2, int numAngles)
{	/*
        Assign pixel values are calculated as a difference between stroke angle and the reference angle
        and vary linearly between 1.0(if the two are equal) and 0.0(if they differ by more than 45 degrees)
      */
    double* minDist = getMinAngleDistance(angles, curAngle, curAngle2, numAngles );
    
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

