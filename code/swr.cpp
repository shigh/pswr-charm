#include <stdlib.h>
#include <vector>
#include <string>
#include "pup_stl.h"
#include "swr.hpp"
#include "SWRCharm.decl.h"

/*readonly*/ CProxy_Main mainProxy;

class Main: public CBase_Main
{
private:

public:

    Main(CkArgMsg* m)
	{
		// doneCells = 0;

		// mainProxy = thisProxy;
		// particlesPerCell = atoi(m->argv[1]);
		// cellDimension = atoi(m->argv[2]);
		// delete m;

		// n_itr = 0;
		// iteration = 0;
		// cells = CProxy_Cell::ckNew(cellDimension, cellDimension);

		// thisProxy.run_simulation();

    }

};

class SWRDomain: public CBase_SWRDomain
{

public:

	SWRDomain(CkMigrateMessage* M) {}

    SWRDomain()
	{

    }

};

#include "SWRCharm.def.h"
