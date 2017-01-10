#include "SketchIO.h"
#include <math.h>
#include <cmath>
#include <vector>
#include <limits>

#define PI 3.14159265
#define BLOCK_SIZE 1024

class FeatureExtractor {
	private:
		Sketch* sketch;
	public:
		// constructor
		FeatureExtractor(Sketch* sketch);
		// a method that finds the angle of a line relative to the x axis
		// between every two consecutive sketch point in every stroke and maps them to the pixels in a feature image
		// using CUDA
		void coords2pixels(int *&angleIndices, double *curAngles, double *curAngles2, int &numOfAngles, double *&pixValues);
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
		// this method finds the max. & min. of x&y coordinates 
		// in the given sketch
		void findExtremum(double &maxX, double &maxY, double &minX, double &minY);
		// downsampler for 24x24 images
		double** downsample(double** image, double gridSize);
		double** gaussianFilter(int hsize, double sigma);
		double** smoothim(double **image, double** fgauss, int hsize, int gridSize);
};

// constructor
FeatureExtractor::FeatureExtractor(Sketch* sketch) : sketch(sketch) {}

// this method sets the sketch to be extracted to the sketch given in the parameter
void FeatureExtractor::setSketch(Sketch* newSketch)
{
	sketch = newSketch;
}

// this method finds the max. & min. of x&y coordinates 
// in the given sketch
void FeatureExtractor::findExtremum(double &maxX, double &maxY, double &minX, double &minY) {
	// get the coordinates of sketch points
	double** sCoords = sketch->getCoords();
	// get the number of points
	int numPoints = sketch->getNumPoints();
	double my_minX,my_minY,my_maxX,my_maxY;
	
	// init. min & max x & y
	minX = numeric_limits<double>::max();
	minY = numeric_limits<double>::max();
	maxX = numeric_limits<double>::min();
	maxY = numeric_limits<double>::min();
	
	// parallelized min & max reduction algorithm 
	#pragma omp parallel private(my_minX,my_minY,my_maxX,my_maxY) num_threads(32)
	{
		// each thread finds their local mins and maxs
		my_minX = numeric_limits<double>::max();
		my_minY = numeric_limits<double>::max();
		my_maxX = numeric_limits<double>::min();
		my_maxY = numeric_limits<double>::min();
		
		#pragma omp for
		for (int i = 0; i < numPoints; ++i) {
			if (my_minX > sCoords[0][i]) {
				my_minX = sCoords[0][i];
			}
			
			if (my_minY > sCoords[1][i]) {
				my_minY = sCoords[1][i];
			}
			
			if (my_maxX < sCoords[0][i]) {
				my_maxX = sCoords[0][i];
			}
			
			if (my_maxY < sCoords[1][i]) {
				my_maxY = sCoords[1][i];
			}
		}
		
		// here we reduce the min & maxs
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
	int numAngles;
	double minX, minY, maxX, maxY;
	int numOfStrokes = normalized->getNumStrokes();
	// we'll use the normalized sketch in feature extraction
	setSketch(normalized);
	// find the extreme coordinates
	findExtremum(maxX,maxY,minX,minY);
	
	// set the angles to be referenced during feature image construction
	double* curAngles = new double[4];
	double* curAngles2 = new double[4];
	double* pixels;
	
	for (int i = 0; i < 4; ++i) {
		curAngles[i] = (double) i*45;
		curAngles2[i] = (double) ((i*45 + 180) % 360);
	}
	
	// construct the pixels
	coords2pixels(angleIndices,curAngles,curAngles2,numAngles,pixels);
	Sketch* transformed = normalized->transform(minX, minY, maxX, maxY);
	
	double** gfilter = gaussianFilter(hsize, sigma);
	
	double **featImage1, **featImage2, **featImage3, **featImage4, **featImage5;
	
	// here, we extract 5 different feature images
	// using OpenMP task parallelism
	// a single task is defined as extraction of a single 12x12 feature image
	#pragma omp parallel num_threads(8)
	{
		#pragma omp sections
		{
			#pragma omp section
			featImage1 = extractFeatureImage(gfilter,pixels, angleIndices, numAngles, transformed,  hsize,  gridSize, false );
			
			#pragma omp section
			featImage2 = extractFeatureImage(gfilter,&pixels[numAngles], angleIndices, numAngles, transformed,  hsize,  gridSize, false );
			
			#pragma omp section
			featImage3 = extractFeatureImage(gfilter,&pixels[2*numAngles], angleIndices, numAngles, transformed,  hsize,  gridSize, false );
			
			#pragma omp section
			featImage4 = extractFeatureImage(gfilter,&pixels[3*numAngles], angleIndices, numAngles, transformed,  hsize,  gridSize, false );
			
			#pragma omp section
			featImage5 = extractFeatureImage(gfilter,pixels, angleIndices, numAngles, transformed,  hsize,  gridSize, true );
		}
	}
	
	delete [] pixels;
	delete [] curAngles;
	delete [] curAngles2;

	// merge & linearize the feature images in a 720-element array
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
	
	for (int i = 0; i <= hsize; ++i) {
		delete [] gfilter[i];
	}
	
	delete resampled;
	delete normalized;
	delete transformed;
	
	delete [] gfilter;
	delete [] featImage1;
	delete [] featImage2;
	delete [] featImage3;
	delete [] featImage4;
	delete [] featImage5;
	
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
			featim[ (int)sCoords[1][strokeStart] ][ (int)sCoords[0][strokeStart] ] = 1;
			featim[ (int)sCoords[1][strokeEnd] ][ (int)sCoords[0][strokeEnd] ] = 1;
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
	
	// we parallelized the convolution operation using OpenMP
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
			
			if (sum > maximum) {
				maximum = sum;
			}
		}
	}
	
	if(maximum != 0)
	{
		// normalization
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

// downsampling of an image through maxpooling
double** FeatureExtractor::downsample(double** image, double gridSize)
{
	double** downsampled = init2Darray(gridSize);
	double curMax;
	
	for(int i = 0; i < 2*gridSize; i += 2)
	{
		for(int j = 0; j < 2*gridSize; j += 2)
		{
			if (image[i][j] < image[i+1][j]) {
				curMax = image[i+1][j];
			}
			else {
				curMax = image[i][j];
			}
			
			if (curMax < image[i][j+1]) {
				curMax = image[i][j+1];
			}
			
			if (curMax < image[i+1][j+1]) {
				curMax = image[i+1][j+1];
			}
			
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
		x1 = sCoords[0][pointIndex];
		x2 = sCoords[0][pointIndex+1];
		y1 = sCoords[1][pointIndex];
		y2 = sCoords[1][pointIndex+1];
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

// CUDA kernel to calculate pixel values for each sketch point
__global__ void coords2pixels_kernel(double *sCoords_x, double *sCoords_y, double *curAngle, double *curAngle2, double *angleThreshold, int *numAngles, int *strokes, int *strokeIndices, int *numOfPoints, double *pixValues) {
	int pt = blockIdx.x*blockDim.x + threadIdx.x;
	
	if (pt < *numOfPoints && strokeIndices[strokes[pt]] < pt) {
		//Get differences both in x and y directions
		double diffy = sCoords_y[pt] - sCoords_y[pt-1];
		double diffx = sCoords_x[pt] - sCoords_x[pt-1];
		double diff,curPixel;
		//Compute angle
		double angle = atan2(diffy,diffx);
		angle = fmod( (angle + 2*PI), (2*PI));
		angle *= 180.0/PI;
		double curAng = *curAngle;
		double curAng2 = *curAngle2;
		
		double curDiff = angle - curAng;
		double curDiff2 = angle - curAng2;
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
		
		pixValues[(((int) curAng)/45)*(*numAngles) + pt-(1+strokes[pt])] = curPixel;
	}
}

// mapping from sketch points to pixel values
void FeatureExtractor::coords2pixels(int *&angleIndices, double *curAngles, double *curAngles2, int &numOfAngles, double *&pixValues) {
	//Get stroke coordinates and indices
	int *sIndices = sketch->getStrokeIndices();
	double **sCoords = sketch->getCoords();
	int numOfStrokes = sketch->getNumStrokes();
	int numOfPoints = sketch->getNumPoints();
	numOfAngles = numOfPoints - numOfStrokes;
	angleIndices = new int[numOfAngles];
	double angleThreshold = 45;
	
	pixValues = new double[numOfAngles*4];
	
	int lastIndex;
	int *strokes = new int[numOfPoints];
	
	for (int str = 0; str < numOfStrokes; ++str) {
		angleIndices[str] = sIndices[str] - str;
		
		if (str == numOfStrokes - 1) {
			lastIndex = numOfPoints;
		}
		else {
			lastIndex = sIndices[str+1];
		}
		
		#pragma omp parallel for num_threads(32)
		for (int pt = sIndices[str]; pt < lastIndex; ++pt) {
			strokes[pt] = str;
		}
	}
	
	// transfer the necessary contents to device memory
	double *sCoords_x_device;
	cudaMalloc( &sCoords_x_device, sizeof(double)*numOfPoints);
	cudaMemcpy( sCoords_x_device, sCoords[0], sizeof(double)*numOfPoints, cudaMemcpyHostToDevice);
	double *sCoords_y_device;
	cudaMalloc( &sCoords_y_device, sizeof(double)*numOfPoints);
	cudaMemcpy( sCoords_y_device, sCoords[1], sizeof(double)*numOfPoints, cudaMemcpyHostToDevice);
	double *curAngles_device;
	cudaMalloc( &curAngles_device, sizeof(double)*4);
	cudaMemcpy( curAngles_device, curAngles, sizeof(double)*4, cudaMemcpyHostToDevice);
	double *curAngles2_device;
	cudaMalloc( &curAngles2_device, sizeof(double)*4);
	cudaMemcpy( curAngles2_device, curAngles2, sizeof(double)*4, cudaMemcpyHostToDevice);
	double *angleThreshold_device;
	cudaMalloc( &angleThreshold_device, sizeof(double));
	cudaMemcpy( angleThreshold_device, &angleThreshold, sizeof(double), cudaMemcpyHostToDevice);
	double *pixValues_device;
	cudaMalloc( &pixValues_device, sizeof(double)*numOfAngles*4);
	int *numOfPoints_device;
	cudaMalloc( &numOfPoints_device, sizeof(int));
	cudaMemcpy( numOfPoints_device, &numOfPoints, sizeof(int), cudaMemcpyHostToDevice);
	int *strokes_device;
	cudaMalloc( &strokes_device, sizeof(int)*numOfPoints);
	cudaMemcpy( strokes_device, strokes, sizeof(int)*numOfPoints, cudaMemcpyHostToDevice);
	int *strokeIndices_device;
	cudaMalloc( &strokeIndices_device, sizeof(int)*numOfStrokes);
	cudaMemcpy( strokeIndices_device, sIndices, sizeof(int)*numOfStrokes, cudaMemcpyHostToDevice);
	int *numAngles_device;
	cudaMalloc( &numAngles_device, sizeof(int));
	cudaMemcpy( numAngles_device, &numOfAngles, sizeof(int), cudaMemcpyHostToDevice);
	
	int itemsPerBlock = numOfPoints / BLOCK_SIZE + (numOfPoints % BLOCK_SIZE == 0 ? 0 : 1);
	
	// task parallelism in here
	// each task is defined as finding pixel values using a specific reference angle
	#pragma omp parallel num_threads(4)
	{
		#pragma omp sections
		{
			#pragma omp section
			coords2pixels_kernel<<<itemsPerBlock,BLOCK_SIZE>>>(sCoords_x_device,sCoords_y_device,&curAngles_device[1],&curAngles2_device[1],angleThreshold_device,numAngles_device,strokes_device,strokeIndices_device,numOfPoints_device,pixValues_device);
			
			#pragma omp section
			coords2pixels_kernel<<<itemsPerBlock,BLOCK_SIZE>>>(sCoords_x_device,sCoords_y_device,&curAngles_device[0],&curAngles2_device[0],angleThreshold_device,numAngles_device,strokes_device,strokeIndices_device,numOfPoints_device,pixValues_device);	
			
			#pragma omp section
			coords2pixels_kernel<<<itemsPerBlock,BLOCK_SIZE>>>(sCoords_x_device,sCoords_y_device,&curAngles_device[2],&curAngles2_device[2],angleThreshold_device,numAngles_device,strokes_device,strokeIndices_device,numOfPoints_device,pixValues_device);
			
			#pragma omp section
			coords2pixels_kernel<<<itemsPerBlock,BLOCK_SIZE>>>(sCoords_x_device,sCoords_y_device,&curAngles_device[3],&curAngles2_device[3],angleThreshold_device,numAngles_device,strokes_device,strokeIndices_device,numOfPoints_device,pixValues_device);
		}
	}
	
	// print CUDA execution error if any
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) cout << "CUDA Execution Error: " << cudaGetErrorString(err) << endl;
	
	// get the pixel contents from the device memory
	cudaMemcpy( pixValues, pixValues_device, sizeof(double)*numOfAngles*4, cudaMemcpyDeviceToHost);
	
	// freeing the arrays resided in the device memory
	cudaFree( sCoords_x_device);
	cudaFree( sCoords_y_device);
	cudaFree( curAngles_device);
	cudaFree( curAngles2_device);
	cudaFree( angleThreshold_device);
	cudaFree( numAngles_device);
	cudaFree( numOfPoints_device);
	cudaFree( strokes_device);
	cudaFree( strokeIndices_device);
	
	delete [] strokes;
}
