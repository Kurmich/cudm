# cuDM
cuDM = CUDA + [IDM](http://rationale.csail.mit.edu/publications/Ouyang2009IJCAI.pdf)

To practice different parallel programming models (OpenMP and CUDA), we wrote cuDM, a CUDA-accelerated IDM feature extractor. 

## How to Use
* Choose which version to use. *openmp* folder contains the pure OpenMP implementation and *cuda_openmp* version includes the hybrid (OpenMP + CUDA) implementation. If you don't want to use a parallel implementation, you can use the serial implementation available under *serial* directory.
* Once you 'cd'ed into one of the above directories, build the code using 'make' command. Make sure your graphics card supports CUDA, nvcc is installed with OpenMP support if you prefer using CUDA version. If you would like to use pure OpenMP version, your g++ compiler must support OpenMP too.
* Type './cudm' to see the usage. Enjoy recognizing sketches! :P Please note that there are some example sketches in *benchmark_sketches* directory to see your executable working.

## Which version is faster?
One would expect CUDA is the fastest. However, the reality is different. As far as we experimented, sketch is not a big data, therefore sending it to the GPU memory takes longer than processing it on the GPU. The fastest ones are pure OpenMP and the serial implementations (sometimes there might be some thread startup delays on the OpenMP one though). But it's good to note that all these are faster than the [Python implementation](http://github.com/ozymaxx/sketchfe). 

Feel free to open an issue if you find a bug.
