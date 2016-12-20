#include <string>
#include <iostream>
#include <vector>
using namespace std;

class Point
{
	string pid;
	double time;
	double x;
	double y;

public:
	Point(string newPid, double newTime, double newX, double newY);
	void printContents();
	double getX() const;
	double getY() const;
};


Point::Point(string newPid, double newTime, double newX, double newY)
{
	pid = newPid;
	time = newTime;
	x = newX;
	y = newY;
}

void Point::printContents()
{
	cout <<"point id = " << pid << " time = " << time << " x = " << getX() << " y = " << y << endl;
}

double Point::getX() const
{
	return x;
}

double Point::getY() const
{
	return y;
}


class Stroke
{
	string sid;
	vector<Point> points; 
public:
	Stroke(string newSid);
	void addPoint(Point point);
	vector<pair<double, double> > listCoordinates(); 
	
};


Stroke::Stroke(string newSid)
{
	sid = newSid;
}

void Stroke::addPoint(Point point)
{
	points.push_back(point);
}

vector<pair<double, double> > Stroke::listCoordinates()
{
	vector<pair<double, double> > result;
	typedef vector<Point> vector_type;
	for(vector_type::const_iterator pos = points.begin(); pos != points.end(); ++pos)
	{
		//cout<<pos->getX()<< " " << pos->getY()<<endl;
		result.push_back(make_pair(pos->getX(), pos->getY()));
	}
	return result;
}



class Sketch
{
	string sketch_id;
	vector<Stroke> strokes;
public:
	Sketch(string newSketch_id);
	Sketch(string newSketch_id, const vector<Stroke> &newStrokes);
	void addStrokes(Stroke stroke);
};


Sketch::Sketch(string newSketch_id)
{
	sketch_id = newSketch_id;
}

Sketch::Sketch(string newSketch_id, const vector<Stroke> &newStrokes)
{
	sketch_id = newSketch_id;
	strokes = newStrokes;
}

void Sketch::addStrokes(Stroke stroke)
{
	strokes.push_back(stroke);
}