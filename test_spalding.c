/**
# Spalding wall model on inclined embedded walls (2D)

A velocity field parallel to a plane wall inclined by an angle
`theta` and satisfying Spalding's law with a prescribed friction
velocity `utau` is imposed. *spalding_update()* must then recover
$\tau_w = \rho u_\tau^2$ from the viscous flux scale `wmr*mu*Gt`,
independently of the inclination of the wall. */

#include "embed_Luca.h"
#include "navier-stokes/centered.h"

double nu = 1e-5, utau = 0.05;
face vector muv[];

/**
Inverse of Spalding's law: the dimensionless velocity $u^+$ for a given
$y^+$ (Newton iterations on $y^+ - F(u^+)$). */

static double uplus (double yp)
{
  double kk = spalding_kappa, E = exp (- kk*spalding_B), u = min (yp, 25.);
  for (int it = 0; it < 100; it++) {
    double k = kk*u, e = exp (k) - 1. - k - k*k/2.;
    double R = u + E*(e - k*k*k/6.) - yp, dR = 1. + E*kk*e;
    u -= R/dR;
    if (fabs (R/dR) < 1e-12)
      break;
  }
  return u;
}

int main()
{
  spalding_nu = nu;
  mu = muv;
  init_grid (64);
  run();
}

event init (t = 0)
{
  double angles[3] = {0., 30., 45.};
  int fail = 0;
  for (int a = 0; a < 3; a++) {
    double th = angles[a]*M_PI/180.;
    coord nf = {sin(th), cos(th)}, t = {cos(th), - sin(th)}; // nf: towards the fluid
    solid (cs, fs, (x - 0.5)*nf.x + (y - 0.3)*nf.y);
    foreach_face()
      muv.x[] = fm.x[]*nu;
    foreach()
      if (cs[] > 0.) {
	double d = max ((x - 0.5)*nf.x + (y - 0.3)*nf.y, 0.);
	double ut = utau*uplus (d*utau/nu);
	u.x[] = ut*t.x, u.y[] = ut*t.y;
      }
      else
	u.x[] = u.y[] = 0.;
    u.x.spalding = u.y.spalding = true;
    spalding_update (u, muv, rho);

    /**
    We check the orientation of the normal and compare $\tau_w$
    with the analytical value. */

    double err = 0., tw = utau*utau;
    int n = 0;
    foreach (reduction(max:err) reduction(+:n))
      if (cs[] > 0. && cs[] < 1.) {
	coord nn, p;
	embed_geometry (point, &p, &nn);
	if (cs[] > 0.5 && p.x*nn.x + p.y*nn.y <= 0.) {
	  fprintf (stderr, "p.n <= 0 at %g %g\n", x, y);
	  fail = 1;
	}
	double gx = 0., gy = 0., c = 0.;
	gx = dirichlet_gradient (point, u.x, cs, nn, p, 0., &c) + c*u.x[];
	gy = dirichlet_gradient (point, u.y, cs, nn, p, 0., &c) + c*u.y[];
	double gn = gx*nn.x + gy*nn.y;
	double Gt = sqrt (sq(gx - gn*nn.x) + sq(gy - gn*nn.y));
	if (wmr[] != 1. && wmr[] < spalding_rmax) {
	  double e = fabs (wmr[]*nu*Gt - tw)/tw;
	  if (e > err) err = e;
	  n++;
	}
      }
    fprintf (stderr, "theta = %g: %d cells, max relative error on tau_w = %g\n",
	     angles[a], n, err);
    if (n == 0 || err > 0.1)
      fail = 1;
  }
  fprintf (stderr, fail ? "FAILED\n" : "OK\n");
  exit (fail);
}
