#include "SketchIO.h"
#include <math.h>
#include <cmath>
#include <vector>

#define PI 3.14159265
#define BLOCK_SIZE 512

class FeatureExtractor {
	private:
		Sketch* sketch;
	public:
		FeatureExtractor(Sketch* sketch);
		void coords2angles(int *&strokeIndices, double *&angles, int &numAngles, double &maxX, double &maxY, double &minX, double &minY);
		double* getMinAngleDistance(double *angles, double curAngle, double curAngle2, int numAngles);
		double truncate(double curDiff);
		double* pixelValues(double *angles, double* curAngles, double* curAngles2, int numAngles);
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
    //cout<<"resampled"<<endl;
    //resampled->printContents();
    //cout << "rs ns=" << resampled->getNumStrokes() << ", np=" << resampled->getNumPoints() << endl;
    Sketch* normalized = resampled->normalized();
    //cout<<"normalized"<<endl;
    //normalized->printContents();
    int *angleIndices;
	double *angles;
	int numAngles;
	double minX, minY, maxX, maxY;
	int numOfStrokes = normalized->getNumStrokes();
	setSketch(normalized);
	coords2angles(angleIndices,angles,numAngles, maxX, maxY, minX, minY);
	Sketch* transformed = normalized->transform(minX, minY, maxX, maxY);
	//cout<<"transformed"<<endl;
	//transformed->printContents();
	double** gfilter = gaussianFilter(hsize, sigma);
	
	double* curAngles = new double[4];
	double* curAngles2 = new double[4];
	
	for (int i = 0; i < 4; ++i) {
		curAngles[i] = (double) i*45;
		curAngles2[i] = (double) ((i*45 + 180) % 360);
	}
	
	double *pixels = pixelValues(angles,curAngles,curAngles2,numAngles);
	
	for (int i = 0; i < 4*numAngles; ++i) {
		if ( i % numAngles == 0) {
			cout << endl << endl;
		}
		
		cout << pixels[i] << " ";
	}
	
	//int curAngle = 0;
	//int curAngle2 = (curAngle + 180)%360;
	//double* pixels1 = pixelValues(angles, curAngle, curAngle2 , numAngles);
	double** featImage1 = extractFeatureImage(gfilter,pixels, angleIndices, numAngles, transformed,  hsize,  gridSize, false );
	
	//curAngle = 45;
	//curAngle2 = (curAngle + 180)%360;
	//double* pixels2 = pixelValues(angles, curAngle, curAngle2, numAngles);
	double** featImage2 = extractFeatureImage(gfilter,&pixels[numAngles], angleIndices, numAngles, transformed,  hsize,  gridSize, false );
	
	//curAngle = 90;
	//curAngle2 = (curAngle + 180)%360;
	//double* pixels3 = pixelValues(angles, curAngle, curAngle2, numAngles);
	double** featImage3 = extractFeatureImage(gfilter,&pixels[2*numAngles], angleIndices, numAngles, transformed,  hsize,  gridSize, false );
	
	//curAngle = 135;
	//curAngle2 = (curAngle + 180)%360;
	//double* pixels4 = pixelValues(angles, curAngle, curAngle2, numAngles);
	double** featImage4 = extractFeatureImage(gfilter,&pixels[3*numAngles], angleIndices, numAngles, transformed,  hsize,  gridSize, false );

	double** featImage5 = extractFeatureImage(gfilter,pixels, angleIndices, numAngles, transformed,  hsize,  gridSize, true );
	/*
	for(int i = 0; i < 2*gridSize; ++i)
	{
		for(int j = 0; j < 2*gridSize; ++j)
		{
			cout<<featImage1[i][j]<< " ";
		}
		cout<<endl;
	}

	cout << "pixels = " << endl;
	for (int i = 0; i < numAngles; ++i) {
		cout << pixels1[i] << endl;
	}*/
	
	delete [] pixels;
	delete [] curAngles;
	delete [] curAngles2;

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
			//cout<<"pointsToImage "<<i<<endl;
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


    /*cout<<"featureImage: "<<endl;
	for(int i = 0; i < 2*gridSize; ++i)
	{
		for(int j = 0; j < 2*gridSize; ++j)
		{
			cout<<featim[i][j]<< " ";
		}
		cout<<endl;
	}*/
	
	double** smoothed = smoothim(featim, gfilter, hsize, gridSize);
	/*cout<<"smoothed: "<<endl;
	for(int i = 0; i < 2*gridSize; ++i)
	{
		for(int j = 0; j < 2*gridSize; ++j)
		{
			cout<<smoothed[i][j]<< " ";
		}
		cout<<endl;
	}*/

	double** downsampled = downsample(smoothed, gridSize);
	/*cout<<"downsampled: "<<endl;
	for(int i = 0; i < gridSize; ++i)
	{
		for(int j = 0; j < gridSize; ++j)
		{
			cout<<downsampled[i][j]<< " ";
		}
		cout<<endl;
	}

	cout << "pixels = " << endl;
	for (int i = 0; i < numAngles; ++i) {
		cout << pixels[i] << endl;
	}*/
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
    //cout<<"gaussianFilter size "<<hsize+1<<endl;
    return gfilter;
}

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
	//cout<<"number of angles"<<numOfAngles<<endl;
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
		//cout<< " x1 " << x1 << " x2 " <<x2 << " y1 "<<y1 << " y2 " << y2 << endl;
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
  int dx = abs(x2 - x1);
  int dy = abs(y2 - y1);
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
      cur = (i + dx)%dx;
      while(cur < 0)
      {
        cur = (cur + dx)%dx;
      }
      arr.push_back(cur);
     // cout<<cur<<endl;
    }
    for(int i = 1; i < arr.size(); ++i)
      {
        cur = arr[i] - arr[i-1];
        cur = cur >=0?1:0;
        cumSum += cur;
        q.push_back(cumSum);
       // cout<<cumSum<<endl;
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
  //cout<<"Hello"<<endl;
 // cout<<"bresenham angle index "<<angleIndex<<endl;
  for(int i = 0; i < x.size(); ++i)
  {
		//cout<<"bresenham "<<x[i]<<" "<<y[i]<<endl;
        if(image[ y[i] ][ x[i] ] < pixels[angleIndex])
        {
        	image[ y[i] ][ x[i] ] = pixels[angleIndex];
        }
  }
  //cout<<dx<< " " <<dy<<endl;
}


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

__global__ void pixelValues_kernel(double *angles, double *curAngles, double *curAngles2, int *numAngles, double *pixValues, double *angleThreshold) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	if (i < 4*(*numAngles)) {
		double curAngle,curAngle2;
		
		if (i < *numAngles) {
			curAngle = curAngles[0];
			curAngle2 = curAngles2[0];
		}
		else if (i < 2*(*numAngles)) {
			curAngle = curAngles[1];
			curAngle2 = curAngles2[1];
		}
		else if (i < 3*(*numAngles)) {
			curAngle = curAngles[2];
			curAngle2 = curAngles2[2];
		}
		else {
			curAngle = curAngles[3];
			curAngle2 = curAngles2[3];
		}
		
		double angle = angles[i % (*numAngles)];
		double curPixel, diff;
		
		//get angle distances relative to current angles
		double curDiff = angle - curAngle;
		double curDiff2 = angle - curAngle2;
		//truncate angle distance between 0 and 180
		if ( curDiff >= 180) {
			curDiff -= 360;
		}
		else if ( curDiff <= -180) {
			curDiff += 360;
		}
		
		curDiff = abs(curDiff);
		
		//curDiff = truncate(curDiff);
		if ( curDiff2 >= 180) {
			curDiff2 -= 360;
		}
		else if ( curDiff2 <= -180) {
			curDiff2 += 360;
		}
		
		curDiff2 = abs(curDiff2);
		
		//curDiff2 = truncate(curDiff2);
		//Assign minimum angle distance
		if (curDiff < curDiff2) {
			diff = curDiff;
		}
		else {
			diff = curDiff2;
		}
		
		//diff[i] = min(curDiff, curDiff2);
		
		if(diff <= (*angleThreshold))
		{
			curPixel = 1 - (diff / (*angleThreshold));
		}
		else
		{
			curPixel = 0;
		}
		pixValues[i] = curPixel;
	}
}

double* FeatureExtractor::pixelValues(double* angles, double* curAngles, double* curAngles2, int numAngles)
{	/*
        Assign pixel values are calculated as a difference between stroke angle and the reference angle
        and vary linearly between 1.0(if the two are equal) and 0.0(if they differ by more than 45 degrees)
      */
    //double* minDist = getMinAngleDistance(angles, curAngle, curAngle2, numAngles );
    /*cout << "diffs = " << endl;
	for (int i = 0; i < numAngles; ++i) {
		cout << minDist[i] << endl;
	}*/
	double* pixValues = new double[4*numAngles];
	
	//double curPixel;
	double angleThreshold = 45;
	//bool valid;
	/*for(int i = 0; i < numAngles; ++i)
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
	}*/
	
	double *angles_device;
	cudaMalloc( &angles_device, sizeof(double)*numAngles);
	cudaMemcpy( angles_device, angles, sizeof(double)*numAngles, cudaMemcpyHostToDevice);
	double *curAngles_device;
	cudaMalloc( &curAngles_device, sizeof(double)*4);
	cudaMemcpy( curAngles_device, curAngles, sizeof(double)*4, cudaMemcpyHostToDevice);
	double *curAngles2_device;
	cudaMalloc( &curAngles2_device, sizeof(double)*4);
	cudaMemcpy( curAngles2_device, curAngles2, sizeof(double)*4, cudaMemcpyHostToDevice);
	int *numAngles_device;
	cudaMalloc( &numAngles_device, sizeof(int));
	cudaMemcpy( numAngles_device, &numAngles, sizeof(int), cudaMemcpyHostToDevice);
	double *pixValues_device;
	cudaMalloc( &pixValues_device, sizeof(double)*numAngles*4);
	double *angleThreshold_device;
	cudaMalloc( &angleThreshold_device, sizeof(double));
	cudaMemcpy( angleThreshold_device, &angleThreshold, sizeof(double), cudaMemcpyHostToDevice);
	
	int itemsPerBlock = (numAngles*4) / BLOCK_SIZE + ((numAngles*4) % BLOCK_SIZE == 0 ? 0 : 1);
	pixelValues_kernel<<<itemsPerBlock,BLOCK_SIZE>>>( angles_device, curAngles_device, curAngles2_device, numAngles_device, pixValues_device, angleThreshold_device);
	
	cudaError_t execErr = cudaGetLastError();
	if (execErr != cudaSuccess) cout << "ERROR: " << cudaGetErrorString(execErr) << endl;
	
	cudaMemcpy( pixValues, pixValues_device, sizeof(double)*numAngles*4, cudaMemcpyDeviceToHost);
	
	cudaFree(angles_device);
	cudaFree(curAngles_device);
	cudaFree(curAngles2_device);
	cudaFree(numAngles_device);
	cudaFree(pixValues_device);
	cudaFree(angleThreshold_device);
	
	/*for (int k = 0; k < numAngles; ++k) {
		cout << pixValues[k] << " ";
	}
	cout << endl;*/

	return pixValues;
}

