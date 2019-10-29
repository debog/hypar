#include <stdlib.h>
#include <physicalmodels/burgers.h>

int BurgersCleanup(void *s)
{
  Burgers *physics = (Burgers*) s;

  return(0);
}
