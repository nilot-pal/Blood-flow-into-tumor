/* UDF for Carreau-Yasuda viscosity model of RBC */
#include "udf.h"
DEFINE_PROPERTY(rbc_viscosity,cell,thread)
{
   /* real alpha_plasma = 0.99;  */
   /* real alpha_leukocyte = 1.; */
   real visc_plasma = 0.001;
   /* real visc_leukocyte = 1.; */
   real lambda = 0.11;
   
   real k_not;
   real visc_rbc;
   real m;
   real n;
   real rate = C_STRAIN_RATE_MAG(cell,thread);
   real alpha_rbc = C_VOF(cell,thread);
 
   if (rate > 6.)
       n = 0.8092*pow(alpha_rbc,3) - 0.8246*pow(alpha_rbc,2) - 0.3503*alpha_rbc + 1 ;
       m = 122.28*pow(alpha_rbc,3) - 51.213*pow(alpha_rbc,2) + 16.305*alpha_rbc + 1 ;
   if (rate <= 6.)
       k_not = log(log(1 + pow(lambda*rate,2)))/log(1 + pow(lambda*rate,2)) ;
       n = k_not*(-0.8913*pow(alpha_rbc,3) + 2.0679*pow(alpha_rbc,2) - 1.7814*alpha_rbc) + 1 ;
       m = 70.782*pow(alpha_rbc,3) - 22.454*pow(alpha_rbc,2) + 9.7193*alpha_rbc + 1 ;
       
   visc_rbc = (visc_plasma*m*(pow((1 + pow(lambda*rate,2)),(n-1)/2)) - (1-alpha_rbc)*visc_plasma)/alpha_rbc ;

   return visc_rbc;
}

/* UDF for shape factor of RBCs */
#include "udf.h"
#define lambda 0.11
DEFINE_EXCHANGE_PROPERTY(shape_factor,cell,mix_thread,s_col,f_col)
{
   Thread *thread_l , *thread_s;
   real phi_sf;
   real rate;

/* find the threads for the liquid (primary) */
/* and solid (secondary phases) */

   thread_l = THREAD_SUB_THREAD(mix_thread, s_col) ;  /* liquid phase */
   thread_s = THREAD_SUB_THREAD(mix_thread, f_col) ;  /* solid phase */

   rate = C_STRAIN_RATE_MAG(cell,thread_l);

   if (rate <= 300)
       phi_sf = 1.5*pow((1 + pow(lambda*rate,2)),0.058697);
   else
       phi_sf = 1;

   return phi_sf;
}

/* UDF for pulsatile velocity at inlet of artery */
#include "udf.h"
DEFINE_PROFILE(inlet_velocity,thread,index)
{
   real tp;
   real w = 8.572;
   real a0 = 0.04152;
   real a1 = -0.005889;
   real b1 = 0.01303;
   real a2 = -0.009072;
   real b2 = 0.005077;
   real a3 = -0.003488;
   real b3 = -0.004537;
   real a4 = -0.001737;
   real b4 = -0.0002039;
   real a5 = -4.077e-05;
   real b5 = -0.004262;
   real a6 = 0.001224;
   real b6 = 0.0006473;
   real a7 = 0.0001488;
   real b7 = -0.002382;
   real a8 = 0.001851;
   real b8 = 0.001263;
   face_t f ;
   begin_f_loop(f,thread)
   {
    tp = RP_Get_Real("flow-time") ;
    F_PROFILE(f,thread,index) = a0 + a1*cos(tp*w) + b1*sin(tp*w) +
               a2*cos(2*tp*w) + b2*sin(2*tp*w) + a3*cos(3*tp*w) + b3*sin(3*tp*w) + 
               a4*cos(4*tp*w) + b4*sin(4*tp*w) + a5*cos(5*tp*w) + b5*sin(5*tp*w) + 
               a6*cos(6*tp*w) + b6*sin(6*tp*w) + a7*cos(7*tp*w) + b7*sin(7*tp*w) +
               a8*cos(8*tp*w) + b8*sin(8*tp*w) ;
   }
  end_f_loop(f,thread)
}
