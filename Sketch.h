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
	void printContents() const;
	string getPointId() const;
	double getTime() const;
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

void Point::printContents() const
{
	cout <<"point id = " << getPointId() << " time = " << getTime() << " x = " << getX() << " y = " << getY() << endl;
}

string Point::getPointId() const
{
	return pid;
}

double Point::getTime() const
{
	return time;
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
	string getStrokeId() const;
	vector<Point> getPoints() const;
	vector<pair<double, double> > listCoordinates() const; 
	
};


Stroke::Stroke(string newSid)
{
	sid = newSid;
}

void Stroke::addPoint(Point point)
{
	points.push_back(point);
}

string Stroke::getStrokeId() const
{
	return sid;
}

vector<Point> Stroke::getPoints() const
{
	return points;
}

vector<pair<double, double> > Stroke::listCoordinates() const
{
	//returns coordinates of all points as (x,y) pairs
	vector<pair<double, double> > result;
	typedef vector<Point> vector_type;
	for(vector_type::const_iterator pos = points.begin(); pos != points.end(); ++pos)
	{
		//cout<<(*pos).getX()<< " " << (*pos).getY()<<endl;
		result.push_back(make_pair((*pos).getX(), (*pos).getY()));
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
	vector<pair<double, double> > listCoordinates();
	void printContents() const;
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

vector<pair<double, double> > Sketch::listCoordinates()
{
	//returns coordinates of all points of all strokes as (x,y) pairs
	vector<pair<double, double> > result;
	typedef vector<Stroke> vector_type;
	for(vector_type::const_iterator pos = strokes.begin(); pos != strokes.end(); ++pos)
	{
		//Insert all points of current stroke to resulting array
		result.insert(result.end(), (*pos).listCoordinates().begin(), (*pos).listCoordinates().end() );
	}
	return result;
}


void Sketch::printContents() const
{
	typedef vector<Stroke> stroke_type;
	typedef vector<Point> point_type;
	for(stroke_type::const_iterator pos = strokes.begin(); pos != strokes.end(); ++pos)
	{
		cout<<"Stroke id = " << (*pos).getStrokeId()<<endl;
		for(point_type::const_iterator pos2 = (*pos).getPoints().begin(); pos2 != (*pos).getPoints().end(); ++pos2 )
		{
			cout<< " time = "<< (*pos2).getTime() << " x = " << (*pos2).getX() << " y = " << (*pos2).getY() << endl;
		}
	}
}