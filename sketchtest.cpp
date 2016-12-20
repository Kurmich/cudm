
#include "Sketch.h"
#include <string>
using namespace std;

int main()
{
	Point p("14", 1.0, 2.0, 3.0);
	Point p2("134", 12.0, 24.0, 3.0);
	p.printContents();
	Stroke st("stroke");
	st.addPoint(p);
	st.addPoint(p2);
	vector<pair<double, double> > c = st.listCoordinates();
	cout << "success "<<endl;
	Sketch sk("sketch");
	sk.addStrokes(st);
	return 0;
}