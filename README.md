# cuDM
cuDM = CUDA + [IDM](http://rationale.csail.mit.edu/publications/Ouyang2009IJCAI.pdf)

To practice different parallel programming models (OpenMP and CUDA), we (Kurmanbek and Ozan) implemented the cuDM, it is a CUDA-accelerated IDM feature extractor.

## How to Use
* Choose a version to use. The *openmp* folder contains pure OpenMP implementation and *cuda_openmp* folder includes a hybrid (OpenMP + CUDA) implementation. If you don't want to use any parallel version, you can opt for the serial implementation available under the *serial* directory.
* Once you `cd`ed into one of the above directories, build the code using `make` command. If you prefer to use the CUDA version, make sure that your graphics card supports CUDA, and nvcc is installed with an OpenMP support. For pure OpenMP version use g++ compiler which supports OpenMP.
* Now type `./cudm` to execute. Enjoy recognizing sketches! Please note that there are some example sketches in the *benchmark_sketches* directory to check if your executable is working properly.

## Which version is the fastest?
One would expect CUDA to be the fastest. However, a vectorized representation of sketches is considered as a "small data" (in relation to the image representation), therefore transferring it to the GPU memory becomes a major bottleneck. The fastest versions are pure OpenMP and serial implementations (sometimes there might be some thread startup delays on the OpenMP one). All in all, it's good to note that all of these versions are faster than the [Python implementation](http://github.com/ozymaxx/sketchfe).

Feel free to open an issue if you find a bug.
