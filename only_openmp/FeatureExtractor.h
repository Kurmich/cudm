#include "SketchIO.h"
#include <math.h>
#include <cmath>
#include <vector>
#include <limits>

#define PI 3.14159265

class FeatureExtractor {
	private:
		Sketch* sketch;
	public:
		FeatureExtractor(Sketch* sketch);
		void coords2angles(int *&strokeIndices, double *&angles, int &numAngles);
		double* getMinAngleDistance(double *angles, double curAngle, double curAngle2, int numAngles);
		double truncate(double curDiff);
		double* pixelValues(double *angles, double curAngle, double curAngle2, int numAngles);
		void findExtremum(double &maxX, double &maxY, double &minX, double &minY);
		double* extract();
		void setSketch(Sketch* newSketch);
		vector<int> arange(int a, int b, int step);
		vector<int> cum(vector<int> q,  int x, char op);
		void drawBresenham(  double x1,  double y1, double x2, double y2, double* pixels, int angleIndex, double** &image );
		double** init2Darray(int size);
		double** extractFeatureImage(double** gfilter, double* pixels, int* angleIndices, int numAngles, Sketch* transformed,  int hsize, int gridSize, bool endpt);
		void pointsToImage(double* pixels, int* angleIndices, double** sCoords, int gridSize, int strokeStart, int strokeEnd, int angleStart, int angleEnd, double** &image);
		double** downsample(double** image, double gridSize);
		double** gaussianFilter(int hsize, double sigma);
		double** smoothim(double **image, double** fgauss, int hsize, int gridSize);
};

FeatureExtractor::FeatureExtractor(Sketch* sketch) : sketch(sketch) {}

void FeatureExtractor::setSketch(Sketch* newSketch)
{
	sketch = newSketch;
}

void FeatureExtractor::findExtremum(double &maxX, double &maxY, double &minX, double &minY) {
	double** sCoords = sketch->getCoords();
	int numPoints = sketch->getNumPoints();
	double my_minX,my_minY,my_maxX,my_maxY;
	
	minX = numeric_limits<double>::max();
	minY = numeric_limits<double>::max();
	maxX = numeric_limits<double>::min();
	maxY = numeric_limits<double>::min();
	
	#pragma omp parallel private(my_minX,my_minY,my_maxX,my_maxY) num_threads(32)
	{
		my_minX = numeric_limits<double>::max();
		my_minY = numeric_limits<double>::max();
		my_maxX = numeric_limits<double>::min();
		my_maxY = numeric_limits<double>::min();
		
		#pragma omp for
		for (int i = 0; i < numPoints; ++i) {
			if (my_minX > sCoords[i][0]) {
				my_minX = sCoords[i][0];
			}
			
			if (my_minY > sCoords[i][1]) {
				my_minY = sCoords[i][1];
			}
			
			if (my_maxX < sCoords[i][0]) {
				my_maxX = sCoords[i][0];
			}
			
			if (my_maxY < sCoords[i][1]) {
				my_maxY = sCoords[i][1];
			}
		}
		
		#pragma omp critical
		{
			if (my_minX < minX) {
				minX = my_minX;
			}
			
			if (my_minY < minY) {
				minY = my_minY;
			}
			
			if (my_maxX > maxX) {
				maxX = my_maxX;
			}
			
			if (my_maxY > maxY) {
				maxY = my_maxY;
			}
		}
	}
}

double* FeatureExtractor::extract()
{
	
	double resampleInterval = 50.0;
    double sigma = 10.0;
    int hsize = 4.0;
    int gridSize = 12;
	int fimgSize = gridSize*gridSize;
    int idmFeatureSize = fimgSize*5;
    double* idmFeature = new double[idmFeatureSize];
    Sketch* resampled = sketch->resample(resampleInterval);
    Sketch* normalized = resampled->normalized();
    int *angleIndices;
	double *angles;
	int numAngles;
	double minX, minY, maxX, maxY;
	int numOfStrokes = normalized->getNumStrokes();
	setSketch(normalized);
	findExtremum(maxX,maxY,minX,minY);
	coords2angles(angleIndices,angles,numAngles);
	Sketch* transformed = normalized->transform(minX, minY, maxX, maxY);
	double** gfilter = gaussianFilter(hsize, sigma);

	int curAngle, curAngle2;
	double *pixels1, *pixels2, *pixels3, *pixels4;
	double **featImage1, **featImage2, **featImage3, **featImage4, **featImage5;
	
	#pragma omp parallel num_threads(8) private(curAngle,curAngle2)
	{
		#pragma omp sections
		{
			#pragma omp section
			{
				curAngle = 0;
				curAngle2 = (curAngle + 180)%360;
				pixels1 = pixelValues(angles, curAngle, curAngle2 , numAngles);
				featImage1 = extractFeatureImage(gfilter,pixels1, angleIndices, numAngles, transformed,  hsize,  gridSize, false );
			}
			
			#pragma omp section
			{
				curAngle = 45;
				curAngle2 = (curAngle + 180)%360;
				pixels2 = pixelValues(angles, curAngle, curAngle2, numAngles);
				featImage2 = extractFeatureImage(gfilter,pixels2, angleIndices, numAngles, transformed,  hsize,  gridSize, false );
			}
			
			#pragma omp section
			{
				curAngle = 90;
				curAngle2 = (curAngle + 180)%360;
				pixels3 = pixelValues(angles, curAngle, curAngle2, numAngles);
				featImage3 = extractFeatureImage(gfilter,pixels3, angleIndices, numAngles, transformed,  hsize,  gridSize, false );
			}
			
			#pragma omp section
			{
				curAngle = 135;
				curAngle2 = (curAngle + 180)%360;
				pixels4 = pixelValues(angles, curAngle, curAngle2, numAngles);
				featImage4 = extractFeatureImage(gfilter,pixels4, angleIndices, numAngles, transformed,  hsize,  gridSize, false );
			}
			
			#pragma omp section
			{
				featImage5 = extractFeatureImage(gfilter,pixels1, angleIndices, numAngles, transformed,  hsize,  gridSize, true );
			}
		}
	}

	#pragma omp parallel for num_threads(32)
	for (int i = 0; i < gridSize; ++i) {
		for (int j = 0; j < gridSize; ++j) {
			idmFeature[i*gridSize+j] = featImage1[i][j];
			idmFeature[fimgSize+i*gridSize+j] = featImage2[i][j];
			idmFeature[2*fimgSize+i*gridSize+j] = featImage3[i][j];
			idmFeature[3*fimgSize+i*gridSize+j] = featImage4[i][j];
			idmFeature[4*fimgSize+i*gridSize+j] = featImage5[i][j];
		}
	}
	
    return idmFeature;
}

double** FeatureExtractor::extractFeatureImage(double** gfilter, double* pixels, int* angleIndices, int numAngles, Sketch* transformed, int hsize, int gridSize, bool endpt )
{
	double** featim = init2Darray(2*gridSize);;
	int* sIndices = transformed->getStrokeIndices();
	double** sCoords = transformed->getCoords();
	int numOfPoints = transformed->getNumPoints();  
	int numOfStrokes = transformed->getNumStrokes();
	if(!endpt)
	{
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
			
			pointsToImage(pixels, angleIndices, sCoords, gridSize, strokeStart, strokeEnd, angleStart, angleEnd, featim);
		}
	}
	else
	{
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
			featim[ (int)sCoords[strokeStart][1] ][ (int)sCoords[strokeStart][0] ] = 1;
			featim[ (int)sCoords[strokeEnd][1] ][ (int)sCoords[strokeEnd][0] ] = 1;
		}
	}
	
	double** smoothed = smoothim(featim, gfilter, hsize, gridSize);

	double** downsampled = downsample(smoothed, gridSize);
	
	return featim; 
}

double** FeatureExtractor::gaussianFilter(int hsize, double sigma)
{
	double** gfilter = init2Darray(hsize+1);
	double sum = 0;
	double r, s = 2.0 * sigma * sigma;
	int d = hsize/2;
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

double** FeatureExtractor::smoothim(double **image, double** fgauss, int hsize, int gridSize)
{
	double** result = init2Darray(2*gridSize);
	int ax, ay;
	double sum = 0, maximum = 0;
	
	#pragma omp parallel for num_threads(32)
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
		#pragma omp parallel for num_threads(32)
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
	double curDiff, curDiff2;
	
	#pragma omp parallel for num_threads(32)
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
	
	#pragma omp parallel for num_threads(32)
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

