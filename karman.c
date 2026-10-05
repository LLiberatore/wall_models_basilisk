#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "view.h"
 
scalar f[];
scalar * tracers = {f};
 
double Reynolds = 160.;
int maxlevel = 10;  
 
face vector muv[];
 
double D = 0.125, U0 = 1.;
double xc = 0., yc = 0., zc=0;
 
int local_maxlevel (double x, double y, double z) {
  return (x < xc) ? 9 : 10;
}
 
#include "adapt_wavelet_leave_interface_limited.h"
 
int main(){
  L0 = 8. [1];
  origin (-0.5, -L0/2., -L0/2.);
 
  init_grid (1 << 8);
 
  mu = muv;
 
  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);
 
  run();
}
 
event properties (i++){
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}
 
u.n[left] = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0.);
 
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
 
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);
 
event init (t = 0){
  refine (fabs(sqrt(sq(x - xc) + sq(y - yc) + sq(z - zc)) - D/2.) < 0.15*D &&
          level < local_maxlevel(x, y, z));
  solid (cs, fs, sqrt(sq(x - xc) + sq(y - yc) + sq(z - zc)) - D/2.);
  foreach()
    u.x[] = cs[] ? U0 : 0.;
}
 
event logfile (i++)
  fprintf(stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
 
event snapshot (t = 0.5){
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (name);
}
 
scalar lev[];

event mesh_png (i += 4, t <= 15.)
{
  char name[256];

  foreach()
    lev[] = level;

  view (fov = 5,
        tx = 0,
        ty = 0.,
        width = 1800,
        height = 1200);

  clear();

  // Colora la slice z = 0 con il livello locale della cella
  squares ("lev",
           n = {0,0,1}, alpha = 0.,
           min = 6, max = 12,
           linear = false);

  // Disegna i bordi delle celle sulla stessa slice
  cells (n = {0,0,1}, alpha = 0.);

  // Disegna l'intersezione della geometria embedded con la slice
  draw_vof ("cs");

  sprintf (name, "frames/mesh-%06d.png", i);
  save (name);
}
event adapt (i++){
  adapt_wavelet_leave_interface_limited (
    slist = {u, f},                            // criteri wavelet standard
    vol_frac = {cs},                           // interfaccia da preservare
    max = (double[]){3e-2, 3e-2, 3e-2, 3e-2},  // tolleranze u.x,u.y,f
    minlevel = 5,
    padding = 0                                // margine attorno all'interfaccia
  );
}
