pswr-charm
==========

Pipeline Schwarz Waveform relaxation using Charm++

### Compiling from source

Code can be build using the commands ```make swr``` or ```make pswr```
from the ```./code``` directory. These commands respectively build a
schwarz waveform relaxation executable and a pipeline schwarz waveform
executable.  A full suite of tests is available in ```./code/tests```,
and can be run using the command ```make test```.  Simple SWR and PSWR
test cases can be run using the commands ```make test-swr``` and
```make test-pswr```. The environment variables ```PETSC_DIR```,
```PETSC_ARCH```, and ```CHARMC``` must be set before running
make. You must also ensure that ```LD_LIBRARY_PATH``` contains the
petsc and charm libraries. You may need to build petsc without mpi.

### Running SWR/PSWR

Both executable can be run using the standard charmrun syntax, and
accept command line arguments controlling the number of domains,
domain size and number of iterations. Pass the ```--help``` option to
view a list of accepted parameters. Currently, the subdomain time
stepping logic (relative number of time steps and time adaptive
behavior) is hard coded in the ```SWRDomain``` and ```PSWRDomain```
classes, and needs to be manually adjusted for different load
balancing tests.

### Running on taub

The directory ```./taub``` contains a sample qsub script. It expects
there to exist a ```setup.sh``` script in this directory that set the
environment variables described above.