// 
//
//      H = -1/2 (\nabla_1^2+\nabla_2^2 - 2/|r_1| - 2/|r_2| + 1/|r_1-r_2|
//              Ψ(R) = [ e^{-β(r1A + r2B)} + e^{-β(r1B + r2A)} ] * e^{-      u(r12)}
//              where riA = |ri - d1|, riB = |ri - d2|
//              Jastrow: u(r12) = A*(1 - exp(-r12/F))/r12
//            u' = [ -u + A*exp(-r12/F)/F ] / r12
//               u'' = [ -2 u' - A*exp(-r12/F)/F^2 ] / r12
//
// ∇i φ = -β ( (ri-d1)/|ri-d1| e^{-β|ri-d1|} + (ri-d2)/|ri-d2| e^{-β|ri-d2|} )
// ∇^2 φ = ( (β^2 - 2β/|ri-d1|) e^{-β|ri-d1|} + (β^2 - 2β/|ri-d2|) e^{-β|ri-d2|} )
//        
//     Eloc = -1/2 sum_i [ ∇^2_i lnΨ + (∇_i lnΨ)·(∇_i lnΨ) ] + V(R)        
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double u_r12(double r12, double A, double F) {
    return (A / r12) * (1.0 - exp(-r12 / F));
}

void u_derivs(double r12, double A, double F, double *up, double *upp) {
    if (r12 == 0.0) {
        *up = 0.0;
        *upp = 0.0;
        return;
    }
    double u = u_r12(r12, A, F);
    double expo = exp(-r12 / F);
    double term = A * expo / F;
    *up  = (-u + term) / r12;
    *upp = (-2.0 * (*up) - A * expo / (F*F)) / r12;
}

double compute_lnpsi(double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     double d, double beta, double A, double F) {

    double d1z =  d / 2.0;
    double d2z = -d / 2.0;

    double r1A = sqrt(x1*x1 + y1*y1 + (z1 - d1z)*(z1 - d1z));
    double r1B = sqrt(x1*x1 + y1*y1 + (z1 - d2z)*(z1 - d2z));
    double r2A = sqrt(x2*x2 + y2*y2 + (z2 - d1z)*(z2 - d1z));
    double r2B = sqrt(x2*x2 + y2*y2 + (z2 - d2z)*(z2 - d2z));
    double r12  = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

    double phi1 = exp(-beta * r1A) + exp(-beta * r1B);
    double phi2 = exp(-beta * r2A) + exp(-beta * r2B);

    double mu = u_r12(r12, A, F);

    if (phi1 <= 0.0 || phi2 <= 0.0) return -1e300;
    return log(phi1) + log(phi2) - mu;
}

void phi_grad_lap(double xi, double yi, double zi,
                  double d1z, double d2z, double beta,
                  double *phi, double grad[3], double *lap) {

    double rx1 = xi, ry1 = yi, rz1 = zi - d1z;
    double rx2 = xi, ry2 = yi, rz2 = zi - d2z;

    double r1 = sqrt(rx1*rx1 + ry1*ry1 + rz1*rz1);
    double r2 = sqrt(rx2*rx2 + ry2*ry2 + rz2*rz2);

    double e1 = exp(-beta * r1);
    double e2 = exp(-beta * r2);

    *phi = e1 + e2;

    if (r1 > 0.0) {
        grad[0] = -beta * (rx1 / r1) * e1;
        grad[1] = -beta * (ry1 / r1) * e1;
        grad[2] = -beta * (rz1 / r1) * e1;
    } else {
        grad[0] = grad[1] = grad[2] = 0.0;
    }
    if (r2 > 0.0) {
        grad[0] += -beta * (rx2 / r2) * e2;
        grad[1] += -beta * (ry2 / r2) * e2;
        grad[2] += -beta * (rz2 / r2) * e2;
    }

    double lap1 = 0.0, lap2 = 0.0;
    if (r1 > 0.0) lap1 = (beta*beta - 2.0*beta / r1) * e1;
    if (r2 > 0.0) lap2 = (beta*beta - 2.0*beta / r2) * e2;
    *lap = lap1 + lap2;
}

double local_energy(double x1,double y1,double z1,
                    double x2,double y2,double z2,
                    double d, double beta, double A, double F) {

    double d1z =  d / 2.0;
    double d2z = -d / 2.0;

    /* distances */
    double r1A = sqrt(x1*x1 + y1*y1 + (z1 - d1z)*(z1 - d1z));
    double r1B = sqrt(x1*x1 + y1*y1 + (z1 - d2z)*(z1 - d2z));
    double r2A = sqrt(x2*x2 + y2*y2 + (z2 - d1z)*(z2 - d1z));
    double r2B = sqrt(x2*x2 + y2*y2 + (z2 - d2z)*(z2 - d2z));
    double rx = x1 - x2, ry = y1 - y2, rz = z1 - z2;
    double r12 = sqrt(rx*rx + ry*ry + rz*rz);

    /* potential energy */
    double epot = -1.0/r1A - 1.0/r1B - 1.0/r2A - 1.0/r2B + 1.0/r12 + 1.0/d;

    /* phi and derivatives for each electron */
    double phi1, phi2, lap1, lap2;
    double gradphi1[3], gradphi2[3];
    phi_grad_lap(x1,y1,z1,d1z,d2z,beta,&phi1,gradphi1,&lap1);
    phi_grad_lap(x2,y2,z2,d1z,d2z,beta,&phi2,gradphi2,&lap2);


    double gradlnphi1[3], gradlnphi2[3];
    for (int k=0;k<3;k++){
        gradlnphi1[k] = gradphi1[k] / (phi1>0.0 ? phi1 : 1e-300);
        gradlnphi2[k] = gradphi2[k] / (phi2>0.0 ? phi2 : 1e-300);
    }
    double laplnphi1 = lap1 / (phi1>0.0 ? phi1 : 1e-300);
    double laplnphi2 = lap2 / (phi2>0.0 ? phi2 : 1e-300);
    double g2_1 = gradlnphi1[0]*gradlnphi1[0] + gradlnphi1[1]*gradlnphi1[1] + gradlnphi1[2]*gradlnphi1[2];
    double g2_2 = gradlnphi2[0]*gradlnphi2[0] + gradlnphi2[1]*gradlnphi2[1] + gradlnphi2[2]*gradlnphi2[2];
    laplnphi1 -= g2_1;
    laplnphi2 -= g2_2;

    /* Jastrow derivatives */
    double up, upp;
    u_derivs(r12, A, F, &up, &upp);

    double grad_u1[3], grad_u2[3];
    if (r12 > 0.0) {
        grad_u1[0] =  rx * up / r12; grad_u1[1] =  ry * up / r12; grad_u1[2] =  rz * up / r12;
        grad_u2[0] = -grad_u1[0];    grad_u2[1] = -grad_u1[1];    grad_u2[2] = -grad_u1[2];
    } else {
        grad_u1[0]=grad_u1[1]=grad_u1[2]=0.0;
        grad_u2[0]=grad_u2[1]=grad_u2[2]=0.0;
    }
    double lap_u1 = (r12>0.0) ? (upp + 2.0*up/r12) : 0.0;
    double lap_u2 = lap_u1;

    double gradlnpsi1[3], gradlnpsi2[3];
    for (int k=0;k<3;k++) {
        gradlnpsi1[k] = gradlnphi1[k] - grad_u1[k];
        gradlnpsi2[k] = gradlnphi2[k] - grad_u2[k];
    }
    double laplnpsi1 = laplnphi1 - lap_u1;
    double laplnpsi2 = laplnphi2 - lap_u2;

    
    double sq1 = gradlnpsi1[0]*gradlnpsi1[0] + gradlnpsi1[1]*gradlnpsi1[1] + gradlnpsi1[2]*gradlnpsi1[2];
    double sq2 = gradlnpsi2[0]*gradlnpsi2[0] + gradlnpsi2[1]*gradlnpsi2[1] + gradlnpsi2[2]*gradlnpsi2[2];
    double kinetic = -0.5 * ( (laplnpsi1 + sq1) + (laplnpsi2 + sq2) );

    double eloc = kinetic + epot;
    return eloc;
}
  
int main () {
  int i,n;
  double p,delta,beta,lnpsi,lnpsip;
  double x1,y1,z1,xp1,yp1,zp1;
  double x2,y2,z2,xp2,yp2,zp2;
 
  srand48((unsigned)time(NULL));        // initialization of rng

  n=1e5;                              // number of monte carlo steps
  beta=1.0 ;                            // variational parameter
  delta=2.5;                            // size of the move
  
  x1=0.3;                               // initial position x1
  y1=0.1;                                // initial position y1
  z1=0.01;                               // initial position z1
  x2=0.1;                               // initial position x2
  y2=0.2;                                // initial position y2
  z2=0.3;                               // initial position z2
  
  double d = 10.0;   // internuclear distance (a.u.)
  double A = .0005;  // Jastrow parameters
  double F = sqrt(A);   //
  
  lnpsi = compute_lnpsi(x1, y1, z1, x2, y2, z2, d, beta, A, F); // log psi at R=(x1,y1,z1,x2,y2,z2)


  for (i=0; i<n ; i++) {
    xp1=x1+delta*(drand48()-0.5);       //
    yp1=y1+delta*(drand48()-0.5);       //
    zp1=z1+delta*(drand48()-0.5);       // propose a move to (xp1,yp1,zp1)
    xp2=x2+delta*(drand48()-0.5);       //                   (xp2,yp2,zp2)
    yp2=y2+delta*(drand48()-0.5);       //
    zp2=z2+delta*(drand48()-0.5);       //
    
    lnpsip = compute_lnpsi(xp1, yp1, zp1, xp2, yp2, zp2, d, beta, A, F);  // log psi at Rp
    
    p=fmin(1.0,exp(2*(lnpsip-lnpsi)));  // probability to accept the move
    
    if(p>drand48()){
      x1=xp1;                             //  *****************************
      y1=yp1;                             //  * update electron positions,*
      z1=zp1;                             //  * distances e-n,            *
                                          //  * and log psi               *
      x2=xp2;                             //  *                           *
      y2=yp2;                             //  *                           *
      z2=zp2;                             //  *                           *
      lnpsi=lnpsip;                       //  *****************************
      
    }
    
    double etot = local_energy(x1,y1,z1,x2,y2,z2,d,beta,A,F);
    printf("%f %f\n", etot, p);
  }
}
