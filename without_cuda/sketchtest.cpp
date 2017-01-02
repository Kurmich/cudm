
#include "Sketch.h"
#include <string>
using namespace std;

int main()
{
	Point p("14", 1.0, 2.0, 3.0);
	Point p2("134", 12.0, 24.0, 3.0);
	Point p3("1341", 122.0, 224.0, 3.0);
	Point p4("1342", 122.0, 244.0, 23.0);
	p.printContents();
	Stroke st("stroke1");
	st.addPoint(p);
	st.addPoint(p2);
	Stroke st2("stroke2");
	st2.addPoint(p3);
	st2.addPoint(p4);
	vector<pair<double, double> > c = st.listCoordinates();
	cout << "success "<<endl;
	Sketch sk("sketch");
	sk.addStrokes(st);
	sk.addStrokes(st2);
	sk.printContents();
	return 0;
}