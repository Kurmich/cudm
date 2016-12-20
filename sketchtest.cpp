
#include "Sketch.h"
#include <string>
using namespace std;

int main()
{
	Point p("14", 1.0, 2.0, 3.0);
	p.printContents();
	Stroke st("stroke");
	st.addPoint(p);
	vector<pair<double, double> > c = st.listCoordinates();
	cout << "success "<<endl;
	Sketch sk("sketch");
	return 0;
}