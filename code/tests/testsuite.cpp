#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE pswrcharm
#include <boost/test/unit_test.hpp>
#include <petscmat.h>

class GlobalTestFixture
{

public:

	GlobalTestFixture()
	{
		PetscInitializeNoArguments();
	}

	~GlobalTestFixture()
	{
		PetscFinalize();
	}

};

BOOST_GLOBAL_FIXTURE( GlobalTestFixture );
