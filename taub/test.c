/* getenv example: getting path */
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* getenv */

int main ()
{
	char* pPath;
	pPath = getenv ("PATH");
	if (pPath!=NULL)
		printf ("\n%s\n",pPath);
	return 0;
}
