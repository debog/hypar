/*
  Burgers Equation

Reference:
+ Hundsdorfer & Verwer, "Numerical Solution of Time Dependent Advection-
  Diffusion-Reaction Euqations", Springer-Verlag Berlin Heidelberg, 2010

  du     d [0.5 u^2]
  --  + ------------  =  0
  dt        dx_i

*/

#define _BURGERS_  "burgers"

typedef struct burgers_parameters {

} Burgers;

int    BurgersInitialize (void*,void*);
int    BurgersCleanup    (void*);
