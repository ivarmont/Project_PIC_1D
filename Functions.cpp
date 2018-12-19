#include "Functions.h"
#include <cstdlib>
#include <cmath>

//random number between 0 and 1
double rnd()
{
	return rand()/(double)RAND_MAX;
}

double InvErrorFunction(double yy){

     double ierf=0.0;
     double  s, t, u, w, x, z,y;


     double qa = 9.16461398268964e-01;
     double qb = 2.31729200323405e-01;
     double qc = 4.88826640273108e-01;
     double qd = 1.24610454613712e-01; 
     double q0 = 4.99999303439796e-01; 
     double q1 = 1.16065025341614e-01; 
     double q2 = 1.50689047360223e-01; 
     double q3 = 2.69999308670029e-01; 
     double q4 = -7.28846765585675e-02;

     double pa = 3.97886080735226000e+00; 
     double pb = 1.20782237635245222e-01; 
     double p0 = 2.44044510593190935e-01; 
     double p1 = 4.34397492331430115e-01; 
     double p2 = 6.86265948274097816e-01; 
     double p3 = 9.56464974744799006e-01; 
     double p4 = 1.16374581931560831e+00; 
     double p5 = 1.21448730779995237e+00; 
     double p6 = 1.05375024970847138e+00; 
     double p7 = 7.13657635868730364e-01; 
     double p8 = 3.16847638520135944e-01; 
     double p9 = 1.47297938331485121e-02;
     double p10 = -1.05872177941595488e-01; 
     double p11 = -7.43424357241784861e-02;

     double p12 = 2.20995927012179067e-02; 
     double p13 = 3.46494207789099922e-02; 
     double p14 = 1.42961988697898018e-02; 
     double p15 = -1.18598117047771104e-02; 
     double p16 = -1.12749169332504870e-02; 
     double p17 = 3.39721910367775861e-03; 
     double p18 = 6.85649426074558612e-03; 
     double p19 = -7.71708358954120939e-04; 
     double p20 = -3.51287146129100025e-03; 
     double p21 = 1.05739299623423047e-04; 
     double p22 = 1.12648096188977922e-03;



      y=1.0-yy;
	if(yy <= -0.9999999 ) y=1.0+0.9999998;
	if(yy >=  0.9999999 ) y=1.0-0.9999998;

      z = y;
	if (y > 1.0){  
		z = 2.0 - y ;
	}

      w = qa - log(z);
      u = sqrt(w);
      s = (qc + log(u))/w;
      t = 1.0/(u + qb);
      x = u*(1.0 - s*(0.5 + s*qd)) - ((((q4*t + q3)*t + q2)*t + q1)*t + q0)*t;
      t = pa/(pa + x);
      u = t - 0.5;
      s = (((((((((p22*u + p21)*u + p20)*u +  p19)*u + p18)*u + p17)*u + p16)*u + p15)*u + p14)*u + p13)*u + p12;
      s = ((((((((((((s*u + p11)*u + p10)*u + p9)*u + p8)*u + p7)*u + p6)*u + p5)*u + p4)*u + p3)*u + p2)*u + p1)*u + p0)*t - z*exp(x*x - pb);
      x = x + s*(1.0 + x*s);
      if (y > 1.0) x = -x;

      ierf = x;
	return ierf;
}
  
/* samples random velocity from Maxwellian distribution using Birdsall's method*/
double SampleVel(double v_th)
{
    const int M = 12;    
    double sum = 0;
	for (int i=0;i<M;i++) sum+=rnd();

    return sqrt(0.5)*v_th*(sum-M/2.0)/sqrt(M/12.0);
}

double SampleVel2(double v_th)
{
    double auxx=0.0;
    double vel=0.0;

    auxx=rnd();		

    auxx=auxx*2.0-1.0;	


    vel=v_th*InvErrorFunction(auxx);	

    return vel;
   
}

double dot_product(double *a,double *b, int N)
{
	double sum=0.0;
	for(int i =0;i<=N;i++) { sum+=a[i]*b[i]; }
	return sum;
}

