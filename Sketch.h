#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>

#define BLOCK_SIZE 16

using namespace std;

// sketch abstraction
// IMPORTANT: To use this abstraction, please give the exact # of strokes and points
// in the constructor and add all the strokes & points before extracting the features.
class Sketch {
	private:
		int numPoints,numStrokes;				// # points and strokes
		int ptAdded,strAdded;					// # of points and strokes added so far
		double **coords;						// array of point coordinates
		int *strokeIndices;						// array of starting indices of each stroke
		void getCentroid(double &x, double &y);	// a method to get the centroid
		void getStd(double &x, double &y);		// a method to get the std. in both x & y axes
		double findMaxDistance();				// a method to find the max. distance from the centroid
	public:
		Sketch(int,int);						// constructor
		~Sketch();								// destructor
		void addPoint(double,double);			// point adder method
		void openStroke();						// to open a stroke
		void printContents();					// to print the coordinates of points contained in each stroke
		int getNumPoints();						// returns # of points
		int getNumStrokes();					// returns # of strokes
		int *getStrokeIndices();				// returns the array containing starting points of every stroke
		double **getCoords();					// returns the array of point coordinates
		
		Sketch* resample(double rate);			// sketch resampler
		Sketch* normalized();					// normalization of a sketch
		Sketch* transform(double minX, double minY, double maxX, double maxY);
};


// constructor
Sketch::Sketch(int numPoints,int numStrokes) : numPoints(numPoints), numStrokes(numStrokes), ptAdded(0), strAdded(0) {
	coords = new double*[numPoints];			// allocation of coordinate array
	for ( int i = 0; i < numPoints; ++i) {
		coords[i] = new double[2];
	}
	
	strokeIndices = new int[numStrokes];		// allocation of stroke indices array
}

// destructor
Sketch::~Sketch() {
	for ( int i = 0; i < numPoints; ++i) {		// deallocation of coordinate array
		delete [] coords[i];
	}
	
	delete [] coords;
	
	delete [] strokeIndices;					// deallocation of stroke indices array
}

// returns the centroid of the sketch
void Sketch::getCentroid(double &x, double &y) {
	if ( numPoints > 0) {						// if there are some points
		double xsum = 0;						// then take the avg. of coordinates
		double ysum = 0;						// and return them in parameters
		
		for (int i = 0; i < ptAdded; ++i) {
			xsum += coords[i][0];
			ysum += coords[i][1];
		}
		
		xsum /= numPoints;
		ysum /= numPoints;
		
		x = xsum;
		y = ysum;
	}
	else {
		x = y = 0;
	}
}

// returns the std. in both axes
void Sketch::getStd(double &x, double &y) {
	// if there are some points
	if ( numPoints > 0) {
		double cx,cy;
		getCentroid(cx,cy);										// take the std. in both axes
																// by just applying the defn. of
		double xsqsum = 0,ysqsum = 0;							// standard deviation
		for ( int i = 0; i < ptAdded; ++i) {
			xsqsum += (coords[i][0] - cx)*(coords[i][0] - cx);
			ysqsum += (coords[i][1] - cy)*(coords[i][1] - cy);
		}
		
		xsqsum /= numPoints;
		ysqsum /= numPoints;
		
		x = sqrt(xsqsum);
		y = sqrt(ysqsum);
	}
	else {
		x = y = 0;
	}
}

double Sketch::findMaxDistance() {
	// find the maximum distance from the centroid that can be possible,
	double x,y;
	getCentroid(x,y);
	
	double maxdist = -1;
	double curdist;
	for (int i = 0; i < ptAdded; ++i) {
		curdist = (coords[i][1] - y)*(coords[i][1] - y) + (coords[i][0] - x)*(coords[i][0] - x);
		
		if (curdist > maxdist) {
			maxdist = curdist;
		}
	}
	
	// return this distance (before returning that I took the sqrt of it here)
	return sqrt(maxdist);
}

// normalization of a sketch
Sketch* Sketch::normalized() {
	double cx,cy;
	double stdx,stdy;
	
	getCentroid(cx,cy);
	getStd(stdx,stdy);
	
	// store the normalized sketch in a new one
	Sketch *newSketch = new Sketch(numPoints,numStrokes);
	
	int upperBound;
	
	for ( int i = 0; i < strAdded; ++i) {
		newSketch->openStroke();
		
		if ( i == numStrokes - 1) {
			upperBound = ptAdded;
		}
		else {
			upperBound = strokeIndices[i+1];
		}
		
		// for each pt, translate the point to the origin
		// and normalize coordinates in both axes by their 
		// corresponding std.
		for ( int j = strokeIndices[i]; j < upperBound; ++j) {
			newSketch->addPoint((coords[j][0]-cx)/stdx,(coords[j][1]-cy)/stdy);
		}
	}
	
	// return the resulting sketch
	return newSketch;
}

// sketch resampler
Sketch* Sketch::resample(double rate) {
	int newNumPoints = 0;
	int upperBound,dist;
	
	// take the sampling interval
	double samplingInterval = findMaxDistance()*1.01 / rate;
	
	// we need to create a new sketch, however, we don't yet know
	// how many points will be created and added to this sketch
	// here I count the # of points needed
	for (int i = 0; i < strAdded; ++i) {
		++newNumPoints;
		
		if (i < numStrokes - 1) {
			upperBound = strokeIndices[i+1];
		}
		else {
			upperBound = ptAdded;
		}
		
		for (int j = strokeIndices[i]+1; j < upperBound; ++j) {
			// count the # of points between every point
			dist = sqrt((coords[j][0]-coords[j-1][0])*(coords[j][0]-coords[j-1][0]) + (coords[j][1]-coords[j-1][1])*(coords[j][1]-coords[j-1][1]));
			//cout << newNumPoints << endl;
			newNumPoints += (int) (ceil(dist / samplingInterval));
		}
	}
	
	// create the resampled sketch
	Sketch* resampled = new Sketch(newNumPoints,numStrokes);
	
	double prevx,prevy,sampdistance,cx,cy,angle,newx,newy;
	for (int i = 0; i < strAdded; ++i) {
		// before I start resampling, I need to open a stroke
		resampled->openStroke();
		
		prevx = coords[strokeIndices[i]][0];
		prevy = coords[strokeIndices[i]][1];
		
		resampled->addPoint(prevx,prevy);
		
		if (i < numStrokes - 1) {
			upperBound = strokeIndices[i+1];
		}
		else {
			upperBound = ptAdded;
		}
		
		// keep adding sample points far from a distance of sampling
		// interval, unless the lastly added point is closer than this
		// interval
		for (int j = strokeIndices[i]+1; j < upperBound; ++j) {
			sampdistance = sqrt((coords[j][0]-prevx)*(coords[j][0]-prevx) + (coords[j][1]-prevy)*(coords[j][1]-prevy));
			
			while (sampdistance > samplingInterval) {
				cx = prevx;
				cy = prevy;
				angle = atan2(coords[j][1]-cy, coords[j][0]-cx);
				newx = cos(angle)*samplingInterval + cx;
				newy = sin(angle)*samplingInterval + cy;
				prevx = newx;
				prevy = newy;
				
				// add the new sampled point
				resampled->addPoint(newx,newy);
				
				sampdistance = sqrt((coords[j][0]-prevx)*(coords[j][0]-prevx) + (coords[j][1]-prevy)*(coords[j][1]-prevy));
			}
		}
	}
	
	return resampled;
}

__global__ void transform_kernel_x(double *&coords_transformed, double *ptAdded, double *newRange, double *oldRangeX, double *newMin, double *minX) {
	int el = blockIdx.x*blockDim.x + threadIdx.x;
	
	if ( el < *ptAdded) {
		coords_transformed[2*el] = floor((coords_transformed[2*el] - *minX) * (*newRange) / (*oldRangeX) + *newMin);
	}
}

__global__ void transform_kernel_y(double *&coords_transformed, double *ptAdded, double *newRange, double *oldRangeY, double *newMin, double *minY) {
	int el = blockIdx.x*blockDim.x + threadIdx.x;
	
	if ( el < *ptAdded) {
		coords_transformed[2*el+1] = floor((coords_transformed[2*el+1] - *minY) * (*newRange) / (*oldRangeY) + *newMin);
	}
}

Sketch* Sketch::transform(double minX, double minY, double maxX, double maxY, int *&sIndices, double **&sCoords)
{
	double newMax = 23, newMin = 0;
	double newRange = newMax - newMin;
	double oldRangeX = maxX - minX;
	double oldRangeY = maxY - minY;
	Sketch *transformed= new Sketch(numPoints,numStrokes);
	transformed->ptAdded = this->ptAdded;
	transformed->strAdded = this->strAdded;
	for(int i = 0; i < strAdded; ++i)
	{
		transformed->strokeIndices[i] = strokeIndices[i];
	}

	/*for(int i = 0; i < ptAdded; ++i)
	{
		if(oldRangeX != 0)
		{
			transformed->coords[i][0] = ((coords[i][0] - minX)*newRange/oldRangeX) + newMin;
			transformed->coords[i][0] = floor(transformed->coords[i][0]);
		}
		if(oldRangeY != 0)
		{
			transformed->coords[i][1] = ((coords[i][1] - minY)*newRange/oldRangeY) + newMin;
			transformed->coords[i][1] = floor(transformed->coords[i][1]);
		}
	}*/
	
	double *coords_device;
	cudaMalloc( &coords_device, sizeof(double*)*ptAdded*2);
	
	for (int i = 0; i < ptAdded; ++i) {
		cudaMemcpy( &coords_device[2*i], transformed->coords[i], sizeof(double)*2, cudaMemcpyHostToDevice);
	}
	
	double* ptAdded_device;
	cudaMalloc( &ptAdded_device, sizeof(double));
	cudaMemcpy( ptAdded_device, &ptAdded, sizeof(double), cudaMemcpyHostToDevice);
	double* newRange_device;
	cudaMalloc( &newRange_device, sizeof(double));
	cudaMemcpy( newRange_device, &newRange, sizeof(double), cudaMemcpyHostToDevice);
	double* newMin_device;
	cudaMalloc( &newMin_device, sizeof(double));
	cudaMemcpy( newMin_device, &newMin, sizeof(double), cudaMemcpyHostToDevice);
	
	if (oldRangeX != 0) {
		double* oldRangeX_device;
		cudaMalloc( &oldRangeX_device, sizeof(double));
		cudaMemcpy( oldRangeX_device, &oldRangeX, sizeof(double), cudaMemcpyHostToDevice);
		double* minX_device;
		cudaMalloc( &minX_device, sizeof(double));
		cudaMemcpy( minX_device, &minX, sizeof(double), cudaMemcpyHostToDevice);
		
		transform_kernel_x<<<ptAdded / BLOCK_SIZE + (ptAdded % BLOCK_SIZE == 0 ? 0 : 1),BLOCK_SIZE>>>(coords_device, ptAdded_device, newRange_device, oldRangeX_device, newMin_device, minX_device);
		
		cudaFree(oldRangeX_device);
		cudaFree(minX_device);
	}
	
	if (oldRangeY != 0) {
		double* oldRangeY_device;
		cudaMalloc( &oldRangeY_device, sizeof(double));
		cudaMemcpy( oldRangeY_device, &oldRangeY, sizeof(double), cudaMemcpyHostToDevice);
		double* minY_device;
		cudaMalloc( &minY_device, sizeof(double));
		cudaMemcpy( minY_device, &minY, sizeof(double), cudaMemcpyHostToDevice);
	
		transform_kernel_x<<<ptAdded / BLOCK_SIZE + (ptAdded % BLOCK_SIZE == 0 ? 0 : 1),BLOCK_SIZE>>>(coords_device, ptAdded_device, newRange_device, oldRangeY_device, newMin_device, minY_device);
		
		cudaFree(oldRangeY_device);
		cudaFree(minY_device);
	}
	
	for (int i = 0; i < ptAdded; ++i) {
		cudaMemcpy( &(transformed->coords[i]), &coords_device[2*i], sizeof(double)*2, cudaMemcpyDeviceToHost);
	}
	
	cudaFree(ptAdded_device);
	cudaFree(newRange_device);
	cudaFree(newMin_device);
	
	return transformed;
}

// point adder
void Sketch::addPoint(double x, double y) {
	// if we haven't completely filled up our sketch
	if (ptAdded < numPoints) {
		// add the given point
		coords[ptAdded][0] = x;
		coords[ptAdded++][1] = y;
	}
	else {
		// otherwise give an error message
		cout << "ERROR: Sketch is full!" << endl;
	}
}

// to open a new stroke
void Sketch::openStroke() {
	// if there are un-opened strokes
	if (strAdded < numStrokes) {
		// add the starting index of the new stroke
		strokeIndices[strAdded++] = ptAdded;
	}
	else {
		// otherwise give an error message
		cout << "ERROR: Sketch is full!" << endl;
	}
}

// print the contents of a sketch
void Sketch::printContents() {
	int upperBound;
	
	for ( int i = 0; i < strAdded; ++i) {
		cout << i << "=>";
		
		if ( i == numStrokes - 1) {
			upperBound = ptAdded;
		}
		else {
			upperBound = strokeIndices[i+1];
		}
		
		for ( int j = strokeIndices[i]; j < upperBound; ++j) {
			cout << "(" << coords[j][0] << "," << coords[j][1] << ")";
			
			if ( j < upperBound - 1) {
				cout << "-";
			} 
		}
		
		cout << endl;
	}
} 

// getter methods
int Sketch::getNumPoints() {
	return ptAdded;
}

int Sketch::getNumStrokes() {
	return strAdded;
}

double** Sketch::getCoords() {
	return coords;
}

int* Sketch::getStrokeIndices() {
	return strokeIndices;
}
