#include "FeatureExtractor.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/time.h>

int is_regular_file(const char *path)
{
    struct stat path_stat;
    stat(path, &path_stat);
    return S_ISREG(path_stat.st_mode);
}

int is_dir(const char *path)
{
    struct stat path_stat;
    stat(path, &path_stat);
    return S_ISDIR(path_stat.st_mode);
}

void writeFeatures(double *features) {
	cout << "[";
	
	for (int i = 0; i < 720; ++i) {
		cout << features[i];
		
		if (i < 719) {
			cout << ",";
		}
	}
	
	cout << "]" << endl;
}

// returns the current time
static const double kMicro = 1.0e-6;
double get_time() {
	struct timeval TV;
	struct timezone TZ;
	const int RC = gettimeofday(&TV, &TZ);
	if(RC == -1) {
		printf("ERROR: Bad call to gettimeofday\n");
		return(-1);
	}
	return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );
}

int main(int argc, char *argv[]) {
	if ( argc < 2) {
		cout << "Usage: " << argv[0] << " <sketch_file> or <path_to_directory>" << endl;
		return -1;
	}
	
	if (is_dir(argv[1])) {
		string ext = ".sketch";
		string dirr(argv[1]);
		
		if (dirr[dirr.size()-1] != '/') {
			dirr += "/";
		}
		
		DIR *dir;
		struct dirent *ent;
		double *features;
		SketchIO sio("sio");
		Sketch* sketch;
		
		if ((dir = opendir (argv[1])) != NULL) {
			double startTime = get_time();
			
			while ((ent = readdir (dir)) != NULL) {
				string fname(ent->d_name);
				
				if (fname.size() > 7) {
					if (fname.substr(fname.size()-7,7).compare(ext) == 0) {
						sio.setFileName(dirr+fname);
						cout << fname << endl;
						sketch = sio.read();
						FeatureExtractor fextractor(sketch);
						features = fextractor.extract();
						delete sketch;
						//writeFeatures(features);
						delete [] features;
					}
				}
			}
			
			closedir(dir);
			
			double endTime = get_time();
			
			cout << "Elapsed time: " << (endTime - startTime) << " secs." << endl;
		} 
		else {
			cout << argv[0] << ": Can't open directory" << endl;
			return -1;
		}
	}
	else if (is_regular_file(argv[1])) {
		string ext = ".sketch";
		
		double *features;
		FeatureExtractor fextractor(NULL);
		SketchIO sio("sio");
		Sketch* sketch;
		
		double startTime = get_time();
		
		string fname(argv[1]);
				
		if (fname.size() > 7) {
			if (fname.substr(fname.size()-7,7).compare(ext) == 0) {
				sio.setFileName(fname);
				cout << fname << endl;
				sketch = sio.read();
				fextractor.setSketch(sketch);
				features = fextractor.extract();
				delete sketch;
				writeFeatures(features);
				delete [] features;
			}
		}
		
		double endTime = get_time();
		
		cout << "Elapsed time: " << (endTime - startTime) << " secs." << endl;
	}
	else {
		cout << argv[0] << ": The given path is neither a directory nor a path" << endl;
		return -1;
	}
	
	return 0;
}
