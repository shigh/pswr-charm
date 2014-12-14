pswr-charm
==========

Pipeline Schwarz Waveform relaxation using Charm++

### Compiling from source

Code can be build using the command ```make swr``` or ```make pswr```
from the ```./code``` directory. A full suite of tests is available in
```./code/tests```, and can be run using the command ```make test```.
Simple SWR and PSWR test cases can be run using the commands ```make
test-swr``` and ```make test-swr```. The environment variables
```PETSC_DIR```, ```PETSC_ARCH```, and ```CHARMC``` must be set before
running make. You must also ensure that ```LD_LIBRARY_PATH``` contains
the petsc and charm libraries.

### Running on taub

The directory ```./taub``` contains a sample qsub script. It expects
there to exist a setup.sh script in this directory that set the
environment varibles described above.