#include "hypergeom.h"
#include "py_hypergeom.h"

 double pvalue_enrichment(int pop,int popt,int sam,int samt)
{
	 int suc=popt;
	 int x=samt;
	 int ceiling=MIN(popt,sam);
	 	int ifault;
	 return chyper (false, sam, ceiling, pop, suc, &ifault ) - chyper (false, sam, x-1, pop, suc, &ifault );
}
 