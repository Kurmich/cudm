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
		void coords2pixels(int *&angleIndices, double *curAngles, double *curAngles2, int &numOfAngles, double *&pixValues);
		double* extract();
		void setSketch(Sketch* newSketch);
		vector<int> arange(int a, int b, int step);
		vector<int> cum(vector<int> q,  int x, char op);
		void drawBresenham(  double x1,  double y1, double x2, double y2, double* pixels, int angleIndex, double** &image );
		double** init2Darray(int size);
		double** extractFeatureImage(double** gfilter, double* pixels, int* angleIndices, int numAngles, Sketch* transformed,  int hsize, int gridSize, bool endpt);
		void pointsToImage(double* pixels, int* angleIndices, double** sCoords, int gridSize, int strokeStart, int strokeEnd, int angleStart, int angleEnd, double** &image);
		void findExtremum(double &maxX, double &maxY, double &minX, double &minY);
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
	
	minX = sCoords[0][0];
	minY = sCoords[1][0];
	maxX = sCoords[0][0];
	maxY = sCoords[1][0];
	
	for (int i = 1; i < numPoints; ++i) {
		if (minX > sCoords[0][i]) {
			minX = sCoords[0][i];
		}
		
		if (minY > sCoords[1][i]) {
			minY = sCoords[1][i];
		}
		
		if (maxX < sCoords[0][i]) {
			maxX = sCoords[0][i];
		}
		
		if (maxY < sCoords[1][i]) {
			maxY = sCoords[1][i];
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
    //cout<<"resampled"<<endl;
    //resampled->printContents();
    //cout << "rs ns=" << resampled->getNumStrokes() << ", np=" << resampled->getNumPoints() << endl;
    Sketch* normalized = resampled->normalized();
    //cout<<"normalized"<<endl;
    //normalized->printContents();
    int *angleIndices;
	int numAngles;
	double minX, minY, maxX, maxY;
	int numOfStrokes = normalized->getNumStrokes();
	setSketch(normalized);
	findExtremum(maxX,maxY,minX,minY);
	
	double* curAngles = new double[4];
	double* curAngles2 = new double[4];
	double* pixels;
	
	for (int i = 0; i < 4; ++i) {
		curAngles[i] = (double) i*45;
		curAngles2[i] = (double) ((i*45 + 180) % 360);
	}
	
	coords2pixels(angleIndices,curAngles,curAngles2,numAngles,pixels);
	Sketch* transformed = normalized->transform(minX, minY, maxX, maxY);
	//cout<<"transformed"<<endl;
	//transformed->printContents();
	double** gfilter = gaussianFilter(hsize, sigma);
	
	/*for (int i = 0; i < 4*numAngles; ++i) {
		if ( i % numAngles == 0) {
			cout << endl << endl;
		}
		
		cout << pixels[i] << " ";
	}*/
	
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
			featim[ (int)sCoords[1][strokeStart] ][ (int)sCoords[0][strokeStart] ] = 1;
			featim[ (int)sCoords[1][strokeEnd] ][ (int)sCoords[0][strokeEnd] ] = 1;
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
		x1 = sCoords[0][pointIndex];
		x2 = sCoords[0][pointIndex+1];
		y1 = sCoords[1][pointIndex];
		y2 = sCoords[1][pointIndex+1];
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
	//double angle,diffy,diffx;
	//int curAngleIndex;
	int *strokes = new int[numOfPoints];
	
	/*for (int str = 0; str < sketch->getNumStrokes(); ++str) {
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
			diffy = sCoords[1][pt] - sCoords[1][pt-1];
			diffx = sCoords[0][pt] - sCoords[0][pt-1];
			//Compute angle
			angle = atan2(diffy,diffx);
			angle = fmod( (angle + 2*PI), (2*PI));
			angle *= 180.0/PI;

			//Assign current angle
			angles[curAngleIndex++] = angle;
		}
	}*/
	
	for (int str = 0; str < numOfStrokes; ++str) {
		angleIndices[str] = sIndices[str] - str;
		
		if (str == numOfStrokes - 1) {
			lastIndex = numOfPoints;
		}
		else {
			lastIndex = sIndices[str+1];
		}
		
		for (int pt = sIndices[str]; pt < lastIndex; ++pt) {
			strokes[pt] = str;
		}
	}
	
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
	coords2pixels_kernel<<<itemsPerBlock,BLOCK_SIZE>>>(sCoords_x_device,sCoords_y_device,&curAngles_device[1],&curAngles2_device[1],angleThreshold_device,numAngles_device,strokes_device,strokeIndices_device,numOfPoints_device,pixValues_device);
	coords2pixels_kernel<<<itemsPerBlock,BLOCK_SIZE>>>(sCoords_x_device,sCoords_y_device,&curAngles_device[0],&curAngles2_device[0],angleThreshold_device,numAngles_device,strokes_device,strokeIndices_device,numOfPoints_device,pixValues_device);	
	coords2pixels_kernel<<<itemsPerBlock,BLOCK_SIZE>>>(sCoords_x_device,sCoords_y_device,&curAngles_device[2],&curAngles2_device[2],angleThreshold_device,numAngles_device,strokes_device,strokeIndices_device,numOfPoints_device,pixValues_device);
	coords2pixels_kernel<<<itemsPerBlock,BLOCK_SIZE>>>(sCoords_x_device,sCoords_y_device,&curAngles_device[3],&curAngles2_device[3],angleThreshold_device,numAngles_device,strokes_device,strokeIndices_device,numOfPoints_device,pixValues_device);
	
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) cout << "CUDA Execution Error: " << cudaGetErrorString(err) << endl;
	
	cudaMemcpy( pixValues, pixValues_device, sizeof(double)*numOfAngles*4, cudaMemcpyDeviceToHost);
	
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
