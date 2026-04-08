// 2021/2/1 last editted
// Standard model code

// Impactor composition�喞hondritic 
// Impact velocity�哥eedingZone (Raylie distribution)
// Planet : Earth 
// MO model : T = 1500 K
// LA model : T = 288K, CO2, H2O condense at the critical partial pressure
// Solubility & Partitioning coefficients : oxidized set
// H solubility : Moore model

#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define DBG_LA(tag) \
  printf("DBG_LA %s %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", \
    tag, A, dA, m_a, N1, N2, N3, N1cum_mass, N2cum_mass, N3cum_mass, E_A_ave, E_p_ave, m_mean, H, d_0)

// Global varieties

// Input parameters
double R_t = 6370e5;  // planetary radius[cm] 
double R_t_LA = 6370e5; // Earth
double const M = 5.97e27; // final planetary mass[g]   Earth : 5.97e24 kg, Venus : 4.87e24 kg, Mars : 6.39e23 kg
double const a_start = 0.1; // initial accreted(planetary) mass fraction[planetary mass] late accretion : 0.99, from Mars size to Earth : 0.1
double const a_end = 0.895; // final accreted mass fraction [final planetary mass] 
double const a_GIstart = 0.895; // initial accreted(planetary) mass fraction[planetary mass] 
double const a_GIend = 0.995; // total accreted fraction after GI stage[final planetary mass] 
double const a_LAstart = 0.995; // initial accreted(planetary) mass fraction[planetary mass] 
double const a_LAend = 1; // final accreted mass fraction[final planetary mass] 
double const a_final = 1.01; // final accreted mass fraction[final planetary mass]
double const T_MO = 1500; // surface temperature[K]
double const T_LA = 288; // surface temperature[K]
double const P1_crit = 1e7; // CO2 partial pressure upper limit in CGS��1 [b] = 0.1 [Pa] = 1e-6 bar 
double const X_1_MO = 1.54e-2;  // CO2 fraction in impactor 
double const X_2_MO = 7.45e-3;   // H2O fraction in impactor 
double const X_3_MO = 1.80e-4;  // N2 fraction in impactor 
double const X_1_LA = 1.54e-2;  // CO2 fraction in impactor 
double const X_2_LA = 7.45e-3;   // H2O fraction in impactor 
double const X_3_LA = 1.80e-4;  // N2 fraction in impactor 
double const X_m = 0.325; // metal fraction in impactor
double const X_mantle = 0.675; // Earth mantle mass 67.5% * M (crust + mantle)
double const X_mo = 0.3; // magma ocean / planetary mass fraction 竏� MO thickness setting = 30-40% of mantle
double const P_initial = 1e5; // Initial total pressure[CGS:b]

// Ginant Impactor's composition 
double M_imp = 5.97e26; // impactor's mass
double R_imp = 3390e5; // impactor's radius
double R_impa = 1; // for calculation
double V_imp = 11.2e5; // Impact velocity of Giant Impacts
double V_mesc = 11.2e5; // Mutual escape velocity of Target and Impactor
double R_c_imp = 1864e5; // impactor's mass
double M_c_imp = 1; // Impactor's core mass
double const N_GI = 1;   // Number of Giant Impacts
double const Mc_a_imp = 7.181130e+22; // stable Mars sized differentiated planet (5% Earth mass started partitioning)
double const Mh_a_imp = 2.081862e+21;
double const Mn_a_imp = 6.124103e+21;
double const Mc_s_imp = 2.286164e+21;
double const Mh_s_imp = 5.481392e+22;
double const Mn_s_imp = 1.044456e+20;
double const Mc_mc_imp = 1.741299e+24;
double const Mh_mc_imp = 4.026896e+23;
double const Mn_mc_imp = 1.171397e+22;

// dry + CC12% : impactor composition
double const Xc_cc12 = 4.2e-3;
double const Xh_cc12 = 8.28e-4;
double const Xn_cc12 = 3.09e-4;

// Physical Constants
double pi; // ﾏ calculated in main function
double const k_B = 1.3806581e-16; // Boltzmann constants CGS [erg/K]
double const N_A = 6.0221367e23;  // Avogadro number CGS [/mol]
double const G = 6.6725985e-8;  // Gravitational constants CGS

// Target planets, Impactors
double g; // gravitational acceleration[cm/s^2] 竊� calculated in main function
double d_planet = 5.51; // mean density of planetary material 
double const d_t = 2.63; // density of planetary material on the surface 2.63 g/cc(Earth's granites cf.Shuvalov2009)
double const d_p = 3.32; // Impactor's density 3.32 g/cc (carbonaceous chondrites, basaltic asteroids cf.Shuvalov2009)
double const t_acc = 9e14; // Planet's accretion timescale = 30 Myr
double const t_rain = 6e8; // rain-out timescale ~ 20 years
double v_esc = 11.2e5; // planet's escape velocity[cm/s] calculated in main function
double P2_crit = 1e6; // upper limit on H2O partial pressure[b] 竊� calculated in main function (saturated vapor pressure in case of that temperature)
double A_tot = 6e25; // Impactors' total mass[g] ex.1% of planetary mass縲竊� calculated in main function
double dA = 6e22;   // one step of impactor's mass[g] ﾎ釆｣_imp calculated in main function
double const dA_accuracy = 1e-4;    // impactor's step bin number
double const dA_accuracy_LA = 1e-5;    // impactor's step bin number
double const dA_steprate = 1e-3;  // impactor's mass step accuracy 

double m1 = 44 ; // CO2 molecule mass g/molecule縲竊偵calculated in main function
double m2 = 18 ; // H2O molecule mass g/molecule
double m3 = 28 ; // N2 molecule mass g/molecule

// Other variaties
double X_1 = 1.46667e-2;  // CO2 fraction in impactor 
double X_2 = 0.18e-2;   // H2O fraction in impactor 
double X_3 = 0.04e-2;  // N2 fraction in impactor 
double T = 1500; // surface temperature
double A = 0;     // cumulative impactor mass[g]
double H = 8.5e5;    // planet atmosphere scaleheight[cm] present Earth H(t=0) = 8.5 km 竊� cm, by deNiem+2012
double d_0 = 0.0012;   // atmoshperic density 
double m_a = 1e22;   // atmospheric mass 竊� calculated in main function
double m_1 = 1e20; // CO2 mass in atmosphere 
double m_2 = 1e20; // H2O mass in atmosphere
double m_3 = 1e20; // N2 mass in atmosphere
double m_1_t0 = 1e20; // Initial CO2 mass in the atmosphere 
double m_2_t0 = 1e20; // Initial H2O mass in the atmosphere
double m_3_t0 = 1e20; // Initial N2 mass in the atmosphere
double N1;  // CO2 molecular number 竊� calculated in main function
double N2;  // H2O molecular number
double N3;  // N2 molecular number
double N1_temp;  // temporal CO2 molecule number
double N2_temp;  // temporal H2O molecule number
double N3_temp;  // temporal N2 molecule number
double dN1 = 0;
double dN2 = 0;
double dN3 = 0;
double N1er_mass = 0; // eroded mass
double N2er_mass = 0;
double N3er_mass = 0;
double N1cum = 0; // cumulative CO2 mass
double N2cum = 0;
double N3cum = 0;
double N1cum_mass = 0; // cumulative CO2 mass in all reservoirs : supply - loss
double N2cum_mass = 0; // cumulative H2O mass in all reservoirs : supply - loss
double N3cum_mass = 0; // cumulative N2 mass in all reservoirs : supply - loss
double Ccum_mass = 0; // cumulative C mass in all reservoirs : supply - loss
double Hcum_mass = 0; // cumulative H mass in all reservoirs : supply - loss
double Ncum_mass = 0; // cumulative N mass in all reservoirs : supply - loss
double dN1er = 0;
double dN2er = 0;
double dN3er = 0;
double dN1er_mass = 0;
double dN2er_mass = 0;
double dN3er_mass = 0;
double dCer_mass = 0;
double dHer_mass = 0;
double dNer_mass = 0;
double dMc_c = 0;
double dMh_c = 0;
double dMn_c = 0;
double m_mean = 6.980655e-23; // mean molecular mass <g/molecule> calculated in main function
double m_meanf = 6.980655e-20; // mean molecular mass <g/molecule> calculated in main function
double P1 = 1e5;  // CO2 partial pressure 竊� calculated in main function
double P2 = 1e5;  // H2O partial pressure 竊� calculated in main function
double P3 = 1e5;  // N2 partial pressure 竊� calculated in main function
double P = 1e5; // total pressure

// Erosion efficiency series
double E_A_ave = 0; // atmospheric erosion efficiency averaged over impact velocity and size distribution 
double E_p_ave = 0; // planetary vapor escape efficiency averaged over impact velocity and size distribution 
  double E_A_Vave = 1; // average atmospheric erosion efficiency over impact velocity distribution
  double E_p_Vave = 1; // average planetary vapor escape efficiency over impact velocity distribution
double Ci_V; // normalizing factor for Velocity distribution
double E_A_GI = 0; // atmospheric loss fraction by Ginant impacts


double EZ = 0.5; 
double EZ_f = 1;
  double EZ1 = 0.5;
  double EZ1_f = 1;
  double EZ2 = 0.5;
  double EZ2_f = 1;
  double EZ3 = 0.5;
  double EZ3_f = 1;


// constants
double mn_mean = 30; // average mass number
double mn_a; // mass number of atmosphereic components
  double const mn_C = 12; // mass number of C
  double const mn_CO = 28; // mass number of CO
  double const mn_CO2 = 44; // mass number of CO2
  double const mn_H = 1; // mass number of H
  double const mn_H2O = 18; // mass number of H2O
  double const mn_N = 14; // mass number of N
  double const mn_N2 = 28; // mass number of N2

// partitioning parameters
double Dc, Dh, Dn, Sc, Sh, Sn;

// variables
double R_ms = 1e-6; // metal/silicate fraction GI 
double R_ms_i, R_mmo;
double const CN_asm = 25; // Initial C/N ratio
double CN_s, CN_as; // C/N ration of silicate, silicate + atmosphere
double M_a, M_mo, M_s, M_s_tar, M_s_GI, M_m, dM_m, M_m_temp, M_c, dM_c, Mc_t, Mc_a, Mc_s, Mc_m, Mc_m_temp, Mc_c, Mc_mc, Mh_t, Mh_a, Mh_s, Mh_m, Mh_m_temp, Mh_c, Mh_mc, Mn_t, Mn_a, Mn_s, Mn_m, Mn_m_temp, Mn_c, Mn_mc; // atmospheric, silicate, metal, total mass, magmaocean 
double Mc_t_as, Mc_t_sm, Mc_s_GI, Mh_t_as, Mh_t_sm, Mh_s_GI, Mn_t_as, Mn_t_sm, Mn_s_GI;
double a, a1, a2, b, c, d, R_ta, R_tb; // variables for calculation
int count; // for iteration calculation

// time
clock_t start, end;
time_t t = 0;

double f_P2_crit(double T); // partial pressure limit (temperature T[K]) : saturated vapor pressure�哺urphy+Koop2005_Eq.(7)(10)

double f_m_mean(double N1, double N2, double N3); // each molecules number 竊� mean molecular mass

double f_Pi(double Ni, double m_mean); // molecular number 竊� partial pressure
double f_Ni(double Pi, double m_mean); // partial pressure + temporal mean molecular number 竊� molecular number

double f_N1_crit(double N2, double N3); // CO2 critical partial pressure 竊� CO2 molecular number
double f_N2_crit(double N1, double N3);

void f_N123_check(double *N1_temp, double *N2_temp, double *N3_temp);  // function to check partial pressure >=< critical pressure and adjust molecular numbers

// ﾎｷ, ﾎｶ average functions
void f_etazeta_ave(double *E_A_ave, double *E_p_ave, double d_p, double d_t, double v_esc, double H, double d_0, double (f_mass) (double y), double (f_Velocity) (double y)); // calculate average efficiency over V distributions
void f_etazeta_Dave(double *E_A_Dave, double *E_p_Dave, double V, double d_p, double d_t, double v_esc, double H, double d_0, double (f_mass) (double y)); // calculate average efficiency over D distributions

double f_Vnormal(double v_esc); // calculating normalizing constant for impact distribution

double f_mass3(double y); // impactor mass distribution
double f_mass2(double y); // impactor mass distribution
double f_Velocity_original(double y); // impact velocity raylie distribution : feeding zone
double f_Velocity(double y); // impact velocity distribution

double f1(double y); // impact velocity distribution-1
double f2(double y); // impact velocity distribution-2
double f3(double y); // impact velocity distribution-3
double f4(double y); // impact velocity distribution-4
double f5(double y); // impact velocity distribution-5

void f_etazeta(double *E_A, double *E_p, double D, double V, double d_p, double d_t, double v_esc, double H, double d_0); // calculating atmospheric erosion efficiency and impactor vapor escape efficiency
void f_eta_GI(double *E_A_GI, double V_imp, double M_imp, double v_esc, double M_t); // Schlichting 2015 Atmospheric loss fraction by Ginant Impacts

void ASMpartitioning_H(double R_ms, double Mi_t, double *Mi_a, double *Mi_s, double *Mi_m, double Si, double Di, double mn_a, double mn_i); // equilibration partitioning for Henry's law C (index: 1)
void ASMpartitioning_notH(double R_ms, double Mi_t, double *Mi_a, double *Mi_s, double *Mi_m, double Si, double Di, double mn_a, double mn_i); // equilibration partitioning for H, N (index: 1/2)
void ASpartitioning_H(double Mi_t_as, double *Mi_a, double *Mi_s, double Si, double mn_a, double mn_i); // equilibration between atmosphere and magmaocean for Henry's law
void ASpartitioning_notH(double Mi_t_as, double *Mi_a, double *Mi_s, double Si, double mn_a, double mn_i); // equilibration between atmosphere and magmaocean for Hydrogen (H2O -> H)
void SMpartitioning(double R_ms, double Mi_t_sm, double *Mi_s_GI, double *Mi_m, double Di);

void ASMpartitioning_H_Moore(double Mh_t, double *Mh_a, double *Mh_s, double *Mh_m, double R_ms, double Dh, double R_t, double m_mean);
void ASpartitioning_H_Moore(double Mh_t_as, double *Mh_a, double *Mh_s, double R_t, double m_mean);
double f_Mh_s_Moore_ASM(double Mh_t, double Mh_s, double R_ms, double Dh, double R_t, double m_mean);
double f_Mh_s_Moore_AS(double Mh_t_as, double Mh_s, double R_t, double m_mean);


int main(){

T = T_MO; // set surface temperature for MO stage
X_1 = X_1_MO; // impactor's composition for MO stage
X_2 = X_2_MO;
X_3 = X_3_MO;

printf("Temperature during MO stage: %e K\n", T);
printf("Final planetary mass M: %e kg\n", M * 1e-3);
printf("Impactor mass step accuracy dA_accuracy: %e \n", dA_accuracy);

  pi = 4 * atan(1); // calculate ﾏ 
  P2_crit = f_P2_crit(T); // H2O saturate pressure as a function of Temperature : constant
  m1 = mn_CO2 / N_A; // CO2 molecule mass g/molecule
  m2 = mn_H2O / N_A; // H2O molecule mass g/molecule
  m3 = mn_N2 / N_A; // N2 molecule mass g/molecule

  A = a_start * M; // initial accreted mass
printf("Initial planetary mass: %e g \n", A);
  A_tot = a_end * M;  // Total impactor mass 
  dA = A_tot * dA_accuracy;   // bin width of impactor mass ﾎ釆｣_imp 竊� update in each step :1000 bins constant (10/9 edited)
  dM_m = dA * X_mo * R_ms / (1 + R_ms); // metal liquid increase mass in each step
  dM_c = dA * X_m - dM_m; // metal segregated mass into core in each step
printf("metal liquid mass step: %e g, metal core segregation step: %e g\n", dM_m, dM_c);
    R_ta = - 0.20945 + (1.0 / 3.0) * log10(A / 5.97e27 / 6.41) - 0.0804 * pow((A / 5.97e27 / 6.41), 0.394); // Seager201?(?)
printf("A: %e, R_ta: %e g\n", A, R_ta);
    R_t = pow(10, R_ta) * 3.29 * 6370e5;
    printf("Initial planetary radius R_t: %e km\n", R_t * 1e-5);
  M_s = X_mantle * A; // Initial magma ocean = 30-40% of mantle, Earth mantle mass 67.5% * M (crust + mantle), ~4e24 kg
  M_m = X_m * A;
  M_mo = M_s + M_m; // magma ocean mass
printf("Initial silicate mass: %e g \n", M_s);
  v_esc = pow((2 * G * A) / R_t, 0.5);     // Escape vilocity Earth : v_esc = 11.2 km/s 
printf("Initial escaping velocity: %e cm/s\n", v_esc);
  g = G * A / pow(R_t, 2); // calculate g (gravitational acceleration)
printf("Initial gravitational acceleration: %e cm/s^2\n", g);
  mn_mean = m_mean * N_A;
printf("Initial mean molecular mass of atmosphere: %e g/molecule, mean mass number of atmosphere: %e\n", m_mean, mn_mean);

// Initial C/H/N abundances
  N1cum_mass = M * a_start * X_1; 
  N2cum_mass = M * a_start * X_2; 
  N3cum_mass = M * a_start * X_3; 
  Ccum_mass = M * a_start * X_1 * mn_C / mn_CO2;
  Hcum_mass = M * a_start * X_2 * (mn_H * 2) / mn_H2O;
  Ncum_mass = M * a_start * X_3;

// Equilibrium calculation - partitioning constants

      Dc = 500;
      Dh = 6.5;
      Dn = 20;
      Sc = 1.6e-13; // 1.6 ppm/MPa
      Sh = 1.897e-8;
      Sn = 1e-13; // 1 ppm/MPa
      R_ms_i = M_m / M_s; // initial metal/silicate ratio completely melting
        printf("Metal / silicate ratio : %e \n", R_ms_i);    // metal/silicate ratio

      Mc_t = Ccum_mass;
      Mh_t = Hcum_mass;
      Mn_t = Ncum_mass;

        printf("Equilibrium C: %e, H: %e, N: %e\n", Mc_t, Mh_t, Mn_t);

      ASMpartitioning_H(R_ms_i, Mc_t, &Mc_a, &Mc_s, &Mc_m, Sc, Dc, mn_mean, mn_C);
      ASMpartitioning_H_Moore(Mh_t, &Mh_a, &Mh_s, &Mh_m, R_ms_i, Dh, R_t, m_mean);
      ASMpartitioning_H(R_ms_i, Mn_t, &Mn_a, &Mn_s, &Mn_m, Sn, Dn, mn_mean, mn_N);

        printf("Carbon: %e, %e, %e, Hydrogen: %e, %e, %e, Nitrogen: %e, %e, %e.\n", Mc_a, Mc_s, Mc_m, Mh_a, Mh_s, Mh_m, Mn_a, Mn_s, Mn_m);

      M_m_temp = M_m;
      M_mo = A * X_mo;
      M_m = M_mo * (R_ms / (1 + R_ms));
      Mc_c = Mc_m * (M_m_temp - M_m) / M_m_temp; 
      Mh_c = Mh_m * (M_m_temp - M_m) / M_m_temp;
      Mn_c = Mn_m * (M_m_temp - M_m) / M_m_temp;

      Mc_mc = Mc_m; // C mass in total metal (metal liquid + segregated core)
      Mc_m = Mc_m * M_m / M_m_temp; // C mass in metal liquid in MO
      Mh_mc = Mh_m; // H mass in total metal (metal liquid + segregated core)
      Mh_m = Mh_m * M_m / M_m_temp; // H mass in metal liquid in MO
      Mn_mc = Mn_m; // N mass in total metal (metal liquid + segregated core)
      Mn_m = Mn_m * M_m / M_m_temp; // N mass in metal liquid in MO

      m_1 = Mc_a * mn_CO2 / mn_C;
      N1 = m_1 / m1;
      m_2 = Mh_a * mn_H2O / (mn_H * 2);
      N2 = m_2 / m2;
      m_3 = Mn_a;
      N3 = m_3 / m3;
      m_a = m_1 + m_2 + m_3;

// planet atmosphere setting
      m_mean = f_m_mean(N1, N2, N3); // initial mean molecular mass  
        printf("Initial N1 = %e, N2 = %e, N3 = %e, mean molecular mass = %e\n", N1, N2, N3, m_mean);
      mn_mean = m_mean * N_A;
        printf("Initial mass number of atmosphere after equilibration = %e\n", mn_mean);

      H = ( k_B * T )/ ( g * m_mean );    // Earth's scaleheight H = 8.5 km 竊� cm by deNiem2012 竊� no condensation case : constant
    printf("Scaleheight H: %e km\n", H * 1e-5);
      d_0 = m_a / (4 * pi * (double) pow(R_t, 2) * H); // atmospheric total mass 竊� atmospheric density

      P1 = f_Pi(N1, m_mean); // calculate the partial pressure
      P2 = f_Pi(N2, m_mean);
      P3 = f_Pi(N3, m_mean);

FILE *fp;
  fp=fopen("MOevolution.dat", "wt"); // make and write data file
      printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass, Mc_a, Mh_a, Mn_a, Mc_s, Mh_s, Mn_s, Mc_mc, Mh_mc, Mn_mc, dA, dCer_mass, dHer_mass, dNer_mass, dMc_c, dMh_c, dMn_c, H, Mc_m, Mh_m, Mn_m, Mc_c, Mh_c, Mn_c);  
      fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass, Mc_a, Mh_a, Mn_a, Mc_s, Mh_s, Mn_s, Mc_mc, Mh_mc, Mn_mc, dA, dCer_mass, dHer_mass, dNer_mass, dMc_c, dMh_c, dMn_c, H, Mc_m, Mh_m, Mn_m, Mc_c, Mh_c, Mn_c); 
  fclose(fp);


for( A = a_start * M ; A < A_tot ; A = A + dA) // continue the calculation during A < A_tot (total impactor mass)
  {
      fp=fopen("MOevolution.dat", "a"); // open and write data file

      R_ta = - 0.20945 + (1.0 / 3.0) * log10(A / 5.97e27 / 6.41) - 0.0804 * pow((A / 5.97e27 / 6.41), 0.394);
      R_t = pow(10, R_ta) * 3.29 * 6370e5;
          printf("Planetary radius R_t: %e km\n", R_t * 1e-5);
      M_mo = X_mo * A; // Magma ocean = 30-40% of mantle, Earth mantle mass 67.5% * M (crust + mantle), 4e24 kg
          printf("Magma ocean mass: %e g \n", M_mo);
      M_s = M_mo * (1 / (1 + R_ms));
      M_m = M_mo * (R_ms / (1 + R_ms));
          printf("Silicate mass: %e g, Metal mass: %e g\n", M_s, M_m);
      v_esc = pow((2 * G * A) / R_t, 0.5);     // escape vilocity Earth : v_esc = 11.2 km/s 
          printf("Escaping velocity: %e cm/s\n", v_esc);
      g = G * A / pow(R_t, 2); // calculate g (gravitational acceleration)
          printf("Gravitational acceleration: %e cm/s^2\n", g);
      m_mean = f_m_mean(N1, N2, N3);
          printf("Average molecule mass : %e g / molecule\n", m_mean);
      mn_mean = m_mean * N_A;
        printf("Mean mass number of atmosphere = %e\n", mn_mean);
      H = ( k_B * T ) / (m_mean * g); // update the scaleheight
          printf("Scale height : %e km\n", H * 1e-5);

      d_0 = m_a / (4 * pi * (double) pow(R_t, 2) * H); // atmospheric total mass 竊� atmospheric density
          printf("Air density : %e g/cc\n", d_0);
          if ( isnan(d_0) ) // atmospheric density d_0 = nan? check
          {
            printf("ERROR : Air density d_0 is NAN !!!\n");
            exit(1); // quit the calculation
          }

      f_etazeta_ave( &E_A_ave, &E_p_ave, d_p, d_t, v_esc, H, d_0, f_mass2, f_Velocity); // calculate the average erosion efficiencies over D & V distribution
          printf("Average Eta : %e \nAverage Zeta : %e\n", E_A_ave, E_p_ave);

      dN1er = ( E_p_ave * X_1 / m1 + E_A_ave * N1 / m_a) * dA;
      N1er_mass = N1er_mass + dN1er * m1; // eroded CO2 mass
      N1cum_mass = A * X_1 - N1er_mass; // including core
      Ccum_mass = N1cum_mass * mn_C / mn_CO2; // atm + sil + met + core
      dCer_mass = dN1er * m1 * mn_C / mn_CO2;
    
      dN2er = ( E_p_ave * X_2 / m2 + E_A_ave * N2 / m_a) * dA;
      N2er_mass = N2er_mass + dN2er * m2;
      N2cum_mass = A * X_2 - N2er_mass;
      Hcum_mass = N2cum_mass * (mn_H * 2) / mn_H2O;
      dHer_mass = dN2er * m2 * (mn_H * 2) / mn_H2O;
      
      dN3er = ( E_p_ave * X_3 / m3 + E_A_ave * N3 / m_a) * dA;
      N3er_mass = N3er_mass + dN3er * m3;
      N3cum_mass = A * X_3 - N3er_mass;
      Ncum_mass = N3cum_mass;
      dNer_mass = dN3er * m3;

// Equilibrium calculation
      Mc_t = Ccum_mass - Mc_c; // equivalent C = supplied C - segregated C
      Mh_t = Hcum_mass - Mh_c;
      Mn_t = Ncum_mass - Mn_c;
        printf("Equilibrium C: %e, H: %e, N: %e\n", Mc_t, Mh_t, Mn_t);
        printf("Metal / silicate ratio : %e \n", R_ms);    // output metal/silicate ratio

      ASMpartitioning_H(R_ms, Mc_t, &Mc_a, &Mc_s, &Mc_m, Sc, Dc, mn_mean, mn_C);
      ASMpartitioning_H_Moore(Mh_t, &Mh_a, &Mh_s, &Mh_m, R_ms, Dh, R_t, m_mean);
      ASMpartitioning_H(R_ms, Mn_t, &Mn_a, &Mn_s, &Mn_m, Sn, Dn, mn_mean, mn_N);

        printf("Carbon: %e, %e, %e, Hydrogen: %e, %e, %e, Nitrogen: %e, %e, %e.\n", Mc_a, Mc_s, Mc_m, Mh_a, Mh_s, Mh_m, Mn_a, Mn_s, Mn_m);

      N1 = Mc_a * mn_CO2 / mn_C / m1;
      N2 = Mh_a * mn_H2O / (mn_H * 2) / m2;
      N3 = Mn_a / m3;
      m_mean = f_m_mean(N1, N2, N3);
        printf("N1 = %e, N2 = %e, N3 = %e, mean molecular mass = %e\n", N1, N2, N3, m_mean);

      m_1 = m1 * N1;
      m_2 = m2 * N2;
      m_3 = m3 * N3;
      
      m_a = m_1 + m_2 + m_3;   // atmospheric total mass in the new step

          printf("Atmospheric total mass : %e g\n", m_a);

      Mc_mc = Mc_m + Mc_c; // carbon mass in all metal = C in equivalent metal liquid + C in former segregated core
      dMc_c = Mc_m * dM_c / M_m;
      Mc_c = Mc_c + dMc_c; // cumulative C mass segregated into core 
      Mc_m = Mc_m * (M_m + dM_m) / M_m; 

      Mh_mc = Mh_m + Mh_c;
      dMh_c = Mh_m * dM_c / M_m;
      Mh_c = Mh_c + dMh_c; // next core segregation
      Mh_m = Mh_m * (M_m + dM_m) / M_m;

      Mn_mc = Mn_m + Mn_c;
      dMn_c = Mn_m * dM_c / M_m;
      Mn_c = Mn_c + dMn_c; // next core segregation
      Mn_m = Mn_m * (M_m + dM_m) / M_m;

      P1 = f_Pi(N1, m_mean); // CO2 partial pressure
      P2 = f_Pi(N2, m_mean); // H2O partial pressure
      P3 = f_Pi(N3, m_mean); // N2 partial pressure

      P = g * m_a / (4 * pi * pow(R_t, 2)); // total pressure

      printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass, Mc_a, Mh_a, Mn_a, Mc_s, Mh_s, Mn_s, Mc_mc, Mh_mc, Mn_mc, dA, dCer_mass, dHer_mass, dNer_mass, dMc_c, dMh_c, dMn_c, H, Mc_m, Mh_m, Mn_m, Mc_c, Mh_c, Mn_c);  
      fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass, Mc_a, Mh_a, Mn_a, Mc_s, Mh_s, Mn_s, Mc_mc, Mh_mc, Mn_mc, dA, dCer_mass, dHer_mass, dNer_mass, dMc_c, dMh_c, dMn_c, H, Mc_m, Mh_m, Mn_m, Mc_c, Mh_c, Mn_c);    

  fclose(fp);

     end = clock();
printf("%lu seconds have passed from start\n", (end - start) / CLOCKS_PER_SEC);
  t = time(NULL);
  printf("%s", ctime(&t));
  }


printf("Planetesimal accretion phase end! with planetary mass %e g\n", A);

// Giant Impact phase
      Mc_c = Mc_mc; // all metal fall into core
      Mc_m = 0;
      Mh_c = Mh_mc;
      Mh_m = 0;
      Mn_c = Mn_mc;
      Mn_m = 0;

  A = a_GIstart * M; // initial accreted mass
printf("Giant Impact phase start! with a planetary mass %e g \n", A);
  A_tot = a_GIend * M;  // Total impactor mass at the end of GI phase
  M_imp = (a_GIend - a_GIstart) * M / N_GI; // 10% Earth mass ; Mars-sized impactor
  dA = M_imp;   // each Ginant impactor mass

  R_impa = - 0.20945 + (1.0 / 3.0) * log10(M_imp / 5.97e27 / 6.41) - 0.0804 * pow((M_imp / 5.97e27 / 6.41), 0.394);
  R_imp = pow(10, R_impa) * 3.29 * 6370e5;
      printf("Planetary radius R_t: %e km\n", R_imp * 1e-5);
  R_c_imp = 0.55 * R_imp; // impactor's core radius
  M_c_imp = X_m * dA; // Impactor's core mass
      printf("Impactor core radius R_c_imp: %e km\n", R_c_imp * 1e-5);

  V_mesc = pow((2 * G * (A + M_imp)) / (R_t + R_imp), 0.5); // Mutual escape velocity of Target and Impactor
  V_imp = 1.1 * V_mesc;     // mutual escape velocity
      printf("Impact velocity: %e km/s\n", V_imp * 1e-5);

  M_s_tar = X_mantle * A;
      printf("Target MO mass: %e g\n", M_s_tar);

      A = A + dA; // calculate atmospheric erosion with the total mass with impactor
      fp=fopen("MOevolution.dat", "a"); // open and write data file

      R_ta = - 0.20945 + (1.0 / 3.0) * log10(A / 5.97e27 / 6.41) - 0.0804 * pow((A / 5.97e27 / 6.41), 0.394);
      R_t = pow(10, R_ta) * 3.29 * 6370e5;
          printf("Planetary radius R_t: %e km\n", R_t * 1e-5);
      M_m = M_c_imp;

      R_mmo = pow(1 + 0.25 * 0.45 * R_t / R_c_imp, -3); // metal fraction in the equilibrating plume Rubie+2015
      R_ms = R_mmo / (1 - R_mmo); // metal/silicate ratio in the plume
      printf("GI meta/MO ratio: %e, GI metal/silicate ratio: %e \n", R_mmo, R_ms);

      M_s_GI = M_m / R_ms; // silicate mass in the plume
      M_s = X_mantle * A; // silicate mass in the completely melted magma ocean
      M_mo = M_s + M_m; // completely melted magma ocean total mass
          printf("Magma ocean mass (completely melted): %e g \n", M_mo);
          printf("Silicate mass: %e g, Metal mass: %e g\n", M_s, M_m);
      v_esc = pow((2 * G * A) / R_t, 0.5);     // escape vilocity of the combined body Earth : v_esc = 11.2 km/s for calculation of atmospheric erosion
          printf("Escaping velocity: %e cm/s\n", v_esc);
      g = G * A / pow(R_t, 2); // calculate g (gravitational acceleration)
          printf("Gravitational acceleration: %e cm/s^2\n", g);

// atmospheric combination: target + impactor
      Mc_a = Mc_a + Mc_a_imp; 
      Mh_a = Mh_a + Mh_a_imp; 
      Mn_a = Mn_a + Mn_a_imp; 

      N1 = Mc_a * mn_CO2 / mn_C / m1;
      N2 = Mh_a * mn_H2O / (mn_H * 2) / m2;
      N3 = Mn_a / m3;
      Ccum_mass = Ccum_mass + (Mc_a_imp + Mc_s_imp + Mc_mc_imp);
      N1cum_mass = Ccum_mass * mn_CO2 / mn_C;
      Hcum_mass = Hcum_mass + (Mh_a_imp + Mh_s_imp + Mh_mc_imp);
      N2cum_mass = Hcum_mass * mn_H2O / mn_H;
      Ncum_mass = Ncum_mass + (Mn_a_imp + Mn_s_imp + Mn_mc_imp);
      N3cum_mass = Ncum_mass;

      m_1 = m1 * N1;
      m_2 = m2 * N2;
      m_3 = m3 * N3;
      
      m_a = m_1 + m_2 + m_3;   // atmospheric total mass in the new step
          printf("Atmospheric total mass : %e g\n", m_a);

      m_mean = f_m_mean(N1, N2, N3);
          printf("Average molecule mass : %e g / molecule\n", m_mean);
      mn_mean = m_mean * N_A;
        printf("Mean mass number of atmosphere = %e\n", mn_mean);
      H = ( k_B * T ) / (m_mean * g); // update the scaleheight
          printf("Scale height : %e km\n", H * 1e-5);

      d_0 = m_a / (4 * pi * (double) pow(R_t, 2) * H); // atmospheric total mass 竊� atmospheric density
          printf("Air density : %e g/cc\n", d_0);
          if ( isnan(d_0) ) // atmospheric density d_0 = nan? check
          {
            printf("ERROR : Air density d_0 is NAN !!!\n");
            exit(1); // quit the calculation
          }

      f_eta_GI( &E_A_GI, V_imp, M_imp, v_esc, A); // Schlichting model
          printf("Atmospheric loss fraction in Giant Impacts : %e \n", E_A_GI);

      dN1er = Mc_a * E_A_GI * mn_CO2 / mn_C; 
      N1er_mass = N1er_mass + dN1er * m1; // eroded mass
      N1cum_mass = N1cum_mass - dN1er * m1;
      Ccum_mass = N1cum_mass * mn_C / mn_CO2; // atm + sil + met + core
      dCer_mass = dN1er * m1 * mn_C / mn_CO2;
    
      dN2er = Mh_a * E_A_GI * mn_H2O / (mn_H * 2); 
      N2er_mass = N2er_mass + dN2er * m2; // eroded mass
      N2cum_mass = N2cum_mass - dN2er * m2;
      Hcum_mass = N2cum_mass * (mn_H * 2) / mn_H2O;
      dHer_mass = dN2er * m2 * (mn_H * 2) / mn_H2O;
      
      dN3er = Mn_a * E_A_GI; 
      N3er_mass = N3er_mass + dN3er * m3; // eroded mass
      N3cum_mass = N3cum_mass - dN3er * m3;
      Ncum_mass = N3cum_mass;
      dNer_mass = dN3er * m3;

// Equilibrium calculation in the silicate-metal plume
      Mc_t_sm = Mc_s * (M_s_GI / M_s_tar) + Mc_mc_imp; // total C mass in the silicate-metal plume
      Mh_t_sm = Mh_s * (M_s_GI / M_s_tar) + Mh_mc_imp; // total C mass in the silicate-metal plume
      Mn_t_sm = Mn_s * (M_s_GI / M_s_tar) + Mn_mc_imp; // total C mass in the silicate-metal plume
        printf("Equilibrium C in the plume: %e, H: %e, N: %e\n", Mc_t_sm, Mh_t_sm, Mn_t_sm);
        printf("Metal / silicate ratio : %e \n", R_ms);    // metal/silicate ratio

      SMpartitioning(R_ms, Mc_t_sm, &Mc_s_GI, &Mc_m, Dc);
      SMpartitioning(R_ms, Mh_t_sm, &Mh_s_GI, &Mh_m, Dh);
      SMpartitioning(R_ms, Mn_t_sm, &Mn_s_GI, &Mn_m, Dn);

      Mc_mc = Mc_m + Mc_c; // carbon mass in all metal = C in equivalent metal liquid + C in former segregated core
      dMc_c = Mc_m;
      Mc_c = Mc_c + dMc_c; // Impactor's core would be combined with the target's core
      Mc_m = 0; 

      Mh_mc = Mh_m + Mh_c;
      dMh_c = Mh_m;
      Mh_c = Mh_c + dMh_c; 
      Mh_m = 0;

      Mn_mc = Mn_m + Mn_c;
      dMn_c = Mn_m;
      Mn_c = Mn_c + dMn_c; 
      Mn_m = 0;

// Equilibrium calculation between the remnant magma ocean and atmosphere

      Mc_t_as = Mc_a * (1 - E_A_GI) + Mc_s * (1 - M_s_GI / M_s_tar) + Mc_s_GI + Mc_s_imp; // eroded atmosphere + remnant MO + plume equilibrated silicate + impactor's silicate
      Mh_t_as = Mh_a * (1 - E_A_GI) + Mh_s * (1 - M_s_GI / M_s_tar) + Mh_s_GI + Mh_s_imp;
      Mn_t_as = Mn_a * (1 - E_A_GI) + Mn_s * (1 - M_s_GI / M_s_tar) + Mn_s_GI + Mn_s_imp;

      for(count = 0; count < 5; count++) // iteration for m_mean
      {
      m_meanf = m_mean;

      ASpartitioning_H(Mc_t_as, &Mc_a, &Mc_s, Sc, mn_mean, mn_C); 
      ASpartitioning_H_Moore(Mh_t_as, &Mh_a, &Mh_s, R_t, m_mean);
      ASpartitioning_H(Mn_t_as, &Mn_a, &Mn_s, Sn, mn_mean, mn_N); 

      N1 = Mc_a * mn_CO2 / mn_C / m1;
      N2 = Mh_a * mn_H2O / (mn_H * 2) / m2;
      N3 = Mn_a / m3;
      m_mean = f_m_mean(N1, N2, N3);
      mn_mean = m_mean * N_A;
	  }

      printf("Former mean molecular mass: %e, New mean molecular mass: %e. after loop\n", m_meanf, m_mean);
      printf("N1 = %e, N2 = %e, N3 = %e, mean molecular mass = %e\n", N1, N2, N3, m_mean);

      m_1 = m1 * N1;
      m_2 = m2 * N2;
      m_3 = m3 * N3;
      
      m_a = m_1 + m_2 + m_3;   // atmospheric total mass in the new step

          printf("Atmospheric total mass : %e g\n", m_a);

      P1 = f_Pi(N1, m_mean); // CO2 partial pressure
      P2 = f_Pi(N2, m_mean); // H2O partial pressure
      P3 = f_Pi(N3, m_mean); // N2 partial pressure

      P = g * m_a / (4 * pi * pow(R_t, 2)); // total pressure

      printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass, Mc_a, Mh_a, Mn_a, Mc_s, Mh_s, Mn_s, Mc_mc, Mh_mc, Mn_mc, dA, dCer_mass, dHer_mass, dNer_mass, dMc_c, dMh_c, dMn_c, H, Mc_m, Mh_m, Mn_m, Mc_c, Mh_c, Mn_c);   
      fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass, Mc_a, Mh_a, Mn_a, Mc_s, Mh_s, Mn_s, Mc_mc, Mh_mc, Mn_mc, dA, dCer_mass, dHer_mass, dNer_mass, dMc_c, dMh_c, dMn_c, H, Mc_m, Mh_m, Mn_m, Mc_c, Mh_c, Mn_c);    

  fclose(fp);

    printf("Giant Impact phase end! with a planetary mass %e g \n", A);

     end = clock();
printf("%lu seconds have passed from start\n", (end - start) / CLOCKS_PER_SEC);
  t = time(NULL);
  printf("%s", ctime(&t));


    fp=fopen("MOfinal.dat", "wt"); // make and write data file
      fprintf(fp, "1 %e %e %e %e\n2 %e %e %e %e\n3 %e %e %e %e\n", Xc_cc12, Mc_a / M, Mc_s / M, Mc_mc / M, Xn_cc12, Mn_a / M, Mn_s / M, Mn_mc / M, Xh_cc12, Mh_a / M, Mh_s / M, Mh_mc / M);   // final each abundance
    fclose(fp);

    fp=fopen("ratio.dat", "wt"); // make and write data file
      fprintf(fp, "1 2 3 %e %e %e\n", (Mc_a + Mc_s) / M, (Mn_a + Mn_s) / M, (Mh_a + Mh_s) / M);    // final C/H C/N
    fclose(fp);



// MO end & LA start -----------------------------------------------------
printf("--- MO solidified and LA start ---\n");
     R_t = R_t_LA;  // set planetary radius
     T = T_LA; // set temperature for LA
     X_1 = X_1_LA; // impactor's composition for MO stage
     X_2 = X_2_LA;
     X_3 = X_3_LA;
     N1er_mass = 0; // reset the total eroded mass of CO2 -> for LA model
     N2er_mass = 0;
     N3er_mass = 0;

printf("Temperature during LA: %e K\n", T);
printf("Target diameter R_t: %e km\n", R_t * 1e-5);
printf("CO2 critical partial pressure P1_crit: %e bar\n", P1_crit * 1e-6);
printf("Impactor mass step accuracy dA_steprate: %e \n", dA_steprate);
printf("%e %e %e %e %e %e\n %e %e %e %e %e %e\n %e %e %e %e %e %e\n %e %e %e %e %e %e\n %e %e %e %e %e %e\n %e %e %e %e %e %e\n", m_1, m_2, m_3, m_1_t0, m_2_t0, m_3_t0, N1, N2, N3, N1_temp, N2_temp, N3_temp, dN1, dN2, dN3, N1er_mass, N2er_mass, N3er_mass, N1cum, N2cum, N3cum, N1cum_mass, N2cum_mass, N3cum_mass, Ccum_mass, Hcum_mass, Ncum_mass, dN1er, dN2er, dN3er, dN1er_mass, dN2er_mass, dN3er_mass, dCer_mass, dHer_mass, dNer_mass);

  P2_crit = f_P2_crit(T); // H2O saturate pressure as a function of Temperature : constant

  v_esc = pow((2 * G * M) / R_t, 0.5);  // Escape vilocity of 100%-mass planet
  g = G * M / pow(R_t, 2); // calculate g (gravitational acceleration) 

  A = a_LAstart * M; // LA initial accreted mass = end of MO phase
printf(" LAInitial planetary mass: %e g \n", A); 
  A_tot = a_LAend * M;  // Total impactor mass : end of the calculation
  dA = A_tot * (a_LAend - a_LAstart) * dA_accuracy_LA;   // bin width of impactor mass ﾎ釆｣_imp 

// Initial atmospheric mass INPUT! from MO calculation (total abundance in BSE)
  m_1_t0 = (Mc_a + Mc_s) * mn_CO2 / mn_C; // CO2 mass at the end of MO // full scenario 7.76e23
  m_2_t0 = (Mh_a + Mh_s) * mn_H2O / (mn_H * 2); // H2O mass at the end of MO // full scenario 2.17e24
  m_3_t0 = (Mn_a + Mn_s); // N2 mass at the end of MO // full scenario 6.53e22

// calculate the m_i, Ni, Nicum, Icum_mass, m_a at the MO solidification
  m_1 = m_1_t0; 
  N1 = m_1 / m1;
  N1cum_mass = m_1;
  Ccum_mass = N1cum_mass * 12. / 44.; // reset Ccum_mass
  m_2 = m_2_t0; 
  N2 = m_2 / m2;
  N2cum_mass = m_2;
  Hcum_mass = N2cum_mass * 2. / 18.; // reset Hcum_mass
  m_3 = m_3_t0; 
  N3 = m_3 / m3;
  N3cum_mass = m_3;
  Ncum_mass = N3cum_mass; // reset Ncum_mass
  m_a = m_1 + m_2 + m_3; // total atmosphere mass

  m_mean = f_m_mean(N1, N2, N3);  // initial mean molecule number

      P1 = f_Pi(N1, m_mean);  // calculate partial pressure from number
      P2 = f_Pi(N2, m_mean);
      P3 = f_Pi(N3, m_mean);

  H = k_B * T / ( g * m_mean );    // scale height
  printf("Scale height at the end of MO stage : %e\n", H); 
  
  fp=fopen("LAevolution.dat", "wt");   // make file

      printf("%e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3);     // standard out put (impactor's mass, atmospheric mass, each molecule mass, partial pressure

fclose(fp);

for( A = a_LAstart * M + dA ; A < A_tot ; A = A + dA) // number of impact steps
  {
      fp=fopen("LAevolution.dat", "a");  // open file
      m_mean = f_m_mean(N1, N2, N3);
          printf("Average molecule mass : %e g / molecule\n", m_mean);
      H = (k_B * T ) / (m_mean * g); 
          printf("Scale height : %e km\n", H * 1e-5);
      d_0 = m_a / (4 * pi * (double) pow(R_t, 2) * H);    // m_a -> d_0 atmospheric density (n)
          printf("Air density : %e g/cc\n", d_0);
          if ( isnan(d_0) )   // check d_0
          {
            printf("ERROR  :  Air density d_0 is NAN !!!\n");
            exit(1);  // end up with error
          }
      DBG_LA("begin");

      f_etazeta_ave( &E_A_ave, &E_p_ave, d_p, d_t, v_esc, H, d_0, f_mass3, f_Velocity); //  calculate averaged erosion efficiency
          printf("Average Eta : %e \nAverage Zeta : %e\n", E_A_ave, E_p_ave);

      N1_temp = N1 + ((1 - E_p_ave) * X_1 / m1 - E_A_ave * N1 / m_a) * dA;   // differential equation for N1
      N2_temp = N2 + ((1 - E_p_ave) * X_2 / m2 - E_A_ave * N2 / m_a) * dA;   // differential equation for N2
      N3_temp = N3 + ((1 - E_p_ave) * X_3 / m3 - E_A_ave * N3 / m_a) * dA;   // differential equation for N3

      if (N1_temp < 0)
      {
        N1_temp = 0;
      }
      if (N2_temp < 0)
      {
        N2_temp = 0;
      }
      if (N3_temp < 0)
      {
        N3_temp = 0;
      }
      DBG_LA("erosion");

      dN1 = N1 - N1_temp; // increase molecule number : not used?
      dN1er = ( E_p_ave * X_1 / m1 + E_A_ave * N1 / m_a) * dA; // eroded molecule number in one step
      dN1er_mass = dN1er * m1;
      if (dN1er_mass > m_1)
      {
        dN1er_mass = m_1;
      }
      N1er_mass = N1er_mass + dN1er_mass;
      N1cum_mass = m_1_t0 + (A - M * a_LAstart) * X_1 - N1er_mass; 
      N1cum = N1cum_mass / m1; //cumulative molecule number in surface reservoirs -> for critical pressure limit check
      Ccum_mass = N1cum_mass * 12. / 44.; // total C mass in all surface reservoirs not including core
    
      dN2 = N2 - N2_temp; 
      dN2er = ( E_p_ave * X_2 / m2 + E_A_ave * N2 / m_a) * dA;
      dN2er_mass = dN2er * m2;
      if (dN2er_mass > m_2)
      {
        dN2er_mass = m_2;
      }
      N2er_mass = N2er_mass + dN2er_mass;      
      N2cum_mass = m_2_t0 + (A - M * a_LAstart) * X_2 - N2er_mass;
      N2cum = N2cum_mass / m2;
      Hcum_mass = N2cum_mass * 2. / 18.;
      
      dN3 = N3 - N3_temp; 
      dN3er = ( E_p_ave * X_3 / m3 + E_A_ave * N3 / m_a) * dA;
      dN3er_mass = dN3er * m3;
      if (dN3er_mass > m_3)
      {
        dN3er_mass = m_3;
      }
      N3er_mass = N3er_mass + dN3er_mass;      
      N3cum_mass = m_3_t0 + (A - M * a_LAstart) * X_3 - N3er_mass;
      N3cum = N3cum_mass / m3;
      Ncum_mass = N3cum_mass;
      DBG_LA("cumulative");

      printf("%e %e %e %e %e %e\n %e %e %e %e %e %e\n %e %e %e %e %e %e\n %e %e %e %e %e %e\n %e %e %e %e %e %e\n %e %e %e %e %e %e\n", m_1, m_2, m_3, m_1_t0, m_2_t0, m_3_t0, N1, N2, N3, N1_temp, N2_temp, N3_temp, dN1, dN2, dN3, N1er_mass, N2er_mass, N3er_mass, N1cum, N2cum, N3cum, N1cum_mass, N2cum_mass, N3cum_mass, Ccum_mass, Hcum_mass, Ncum_mass, dN1er, dN2er, dN3er, dN1er_mass, dN2er_mass, dN3er_mass, dCer_mass, dHer_mass, dNer_mass);

      dA = dA_steprate * m_a / fabs((1 - E_p_ave) * X_3 - E_A_ave);   //  calculate dA step 

        printf("Impactor mass step : %e g\n", dA);    // impactor's mass
  
        EZ1_f = EZ1;
        EZ2_f = EZ2;
        EZ3_f = EZ3;
      
      f_N123_check(&N1cum, &N2cum, &N3cum); // check the critical pressure
          printf("Number of molecules : %e %e %e\n", N1cum, N2cum, N3cum);
      N1 = N1cum; //  number of CO2 molecules in atmosphere
      N2 = N2cum;
      N3 = N3cum;

      m_1 = m1 * N1; // CO2 mass in atmosphere
      m_2 = m2 * N2;
      m_3 = m3 * N3;
      
      m_a = m_1 + m_2 + m_3;   // total atmospheric mass in (n+1) step

          printf("Atmospheric total mass : %e g\n", m_a);

      P1 = f_Pi(N1, m_mean);  // Partial pressure
      P2 = f_Pi(N2, m_mean);
      P3 = f_Pi(N3, m_mean);

      P = g * m_a / (4 * pi * pow(R_t, 2)); // total pressure
      DBG_LA("final");

      printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass);     
      fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass);    

  fclose(fp);

}



      fp=fopen("LAfinal.dat", "wt");   // 繝��繧ｿ繝輔ぃ繧､繝ｫ縺ｮ菴懈�

       fprintf(fp, "1 2 3 %e %e %e %e\n", A, Ccum_mass, Ncum_mass, Hcum_mass);    // 陦晉ｪ∝､ｩ菴鍋ｴｯ險郁ｳｪ驥上→螟ｧ豌鈴㍼縺ｮ髢｢菫ゅｒ繝��繧ｿ蛹� 繝�く繧ｹ繝医ヵ繧｡繧､繝ｫ縺ｸ譖ｸ縺崎ｾｼ縺ｿ   

      fclose(fp);

    fp=fopen("LAfinalpattern.dat", "wt"); // make and write data file
      fprintf(fp, "1 %e %e\n2 %e %e\n3 %e %e\n", Xc_cc12, Ccum_mass / M, Xn_cc12, Ncum_mass / M, Xh_cc12, Hcum_mass / M);   // final each abundance
    fclose(fp);


// a_LAend -> a_final : after LA for plotting the evolution of abundances

  A_tot = a_final * M;
    fp=fopen("afterLAevolution.dat", "wt");   // make file

      printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass);   // rewrite A = a_LAend   
      fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass);    

    fclose(fp);

  for( A = a_LAend * M ; A < A_tot ; A = A + dA) // number of impact steps
  {
      fp=fopen("afterLAevolution.dat", "a");  // open file
      m_mean = f_m_mean(N1, N2, N3);
          printf("Average molecule mass : %e g / molecule\n", m_mean);
      H = (k_B * T ) / (m_mean * g); 
          printf("Scale height : %e km\n", H * 1e-5);
      d_0 = m_a / (4 * pi * (double) pow(R_t, 2) * H);    // m_a -> d_0 atmospheric density (n)
          printf("Air density : %e g/cc\n", d_0);
          if ( isnan(d_0) )   // check d_0
          {
            printf("ERROR  :  Air density d_0 is NAN !!!\n");
            exit(1);  // end up with error
          }

      f_etazeta_ave( &E_A_ave, &E_p_ave, d_p, d_t, v_esc, H, d_0, f_mass3, f_Velocity); //  calculate averaged erosion efficiency
          printf("Average Eta : %e \nAverage Zeta : %e\n", E_A_ave, E_p_ave);


      N1_temp = N1 + ((1 - E_p_ave) * X_1 / m1 - E_A_ave * N1 / m_a) * dA;   // differential equation for N1
      N2_temp = N2 + ((1 - E_p_ave) * X_2 / m2 - E_A_ave * N2 / m_a) * dA;   // differential equation for N2
      N3_temp = N3 + ((1 - E_p_ave) * X_3 / m3 - E_A_ave * N3 / m_a) * dA;   // differential equation for N3

      if (N1_temp < 0)
      {
        N1_temp = 0;
      }
      if (N2_temp < 0)
      {
        N2_temp = 0;
      }
      if (N3_temp < 0)
      {
        N3_temp = 0;
      }

      dN1 = N1 - N1_temp;
      dN1er = ( E_p_ave * X_1 / m1 + E_A_ave * N1 / m_a) * dA;
      dN1er_mass = dN1er * m1;
      if (dN1er_mass > m_1)
      {
        dN1er_mass = m_1;
      }
      N1er_mass = N1er_mass + dN1er_mass;
      N1cum_mass = m_1_t0 + (A - M * a_LAstart) * X_1 - N1er_mass;
      N1cum = N1cum_mass / m1;
      Ccum_mass = N1cum_mass * 12. / 44.;
    
      dN2 = N2 - N2_temp;
      dN2er = ( E_p_ave * X_2 / m2 + E_A_ave * N2 / m_a) * dA;
      dN2er_mass = dN2er * m2;
      if (dN2er_mass > m_2)
      {
        dN2er_mass = m_2;
      }
      N2er_mass = N2er_mass + dN2er_mass;      
      N2cum_mass = m_2_t0 + (A - M * a_LAstart) * X_2 - N2er_mass;
      N2cum = N2cum_mass / m2;
      Hcum_mass = N2cum_mass * 2. / 18.;
      
      dN3 = N3 - N3_temp;
      dN3er = ( E_p_ave * X_3 / m3 + E_A_ave * N3 / m_a) * dA;
      dN3er_mass = dN3er * m3;
      if (dN3er_mass > m_3)
      {
        dN3er_mass = m_3;
      }
      N3er_mass = N3er_mass + dN3er_mass;      
      N3cum_mass = m_3_t0 + (A - M * a_LAstart) * X_3 - N3er_mass;
      N3cum = N3cum_mass / m3;
      Ncum_mass = N3cum_mass;

      dA = dA_steprate * m_a / fabs((1 - E_p_ave) * X_3 - E_A_ave);   //  calculate dA step 


        printf("Impactor mass step : %e g\n", dA);    // impactor's mass
  
        EZ1_f = EZ1;
        EZ2_f = EZ2;
        EZ3_f = EZ3;

      
      f_N123_check(&N1cum, &N2cum, &N3cum); // check the critical pressure
          printf("Number of molecules : %e %e %e\n", N1_temp, N2_temp, N3_temp);
      N1 = N1cum;
      N2 = N2cum;
      N3 = N3cum;

      m_1 = m1 * N1;
      m_2 = m2 * N2;
      m_3 = m3 * N3;
      
      m_a = m_1 + m_2 + m_3;   // total atmospheric mass in (n+1) step

          printf("Atmospheric total mass : %e g\n", m_a);

///*
      P1 = f_Pi(N1, m_mean);  // Partial pressure
      P2 = f_Pi(N2, m_mean);
      P3 = f_Pi(N3, m_mean);
//*/

      P = g * m_a / (4 * pi * pow(R_t, 2)); // total pressure

      printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass);     
      fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", A, m_a, m_1, m_2, m_3, P1, P2, P3, E_A_ave, E_p_ave, N1er_mass, N2er_mass, N3er_mass, Ccum_mass, Hcum_mass, Ncum_mass);    

  fclose(fp);



     end = clock();
printf("%lu seconds have passed from start\n", (end - start) / CLOCKS_PER_SEC);
  t = time(NULL);
  printf("%s", ctime(&t));
}

          printf("Mission accomplished\n");
  return 0;

} // end of main function

double f_P2_crit(double T)  // temperature T[K] 竊� saturated vapor pressure�哺urphy+Koop2005_Eq.(7)(10) 竊� [Pa] 竊� [b]
{
  if(T < 273.16)  // Triple point
  {
    return ( exp( 9.550426 - 5723.265 / T + 3.53068 * log(T) - 0.00728332 * T ) );  // Eq.(7) solid - gas equilibrium
  }
  else
  {
    double lnP = 54.842763 - 6763.22 / T - 4.210 * log(T) + 0.000367 * T
                 + tanh( 0.0415 * ( T - 218.8 ) ) * ( 53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T );   // Eq.(10) liquid - gas equilibrium
    return ( exp( lnP ) * 10); // 1 Pa 竊� 10 b
  }
}

double f_m_mean(double N1, double N2, double N3) // each molecules number 竊� mean molecular mass
{
  return ( (m1 * N1 + m2 * N2 + m3 * N3) / (N1 + N2 + N3) );
}

double f_Pi(double Ni, double m_mean) // molecular number 竊� partial pressure
{
  return ( (g * m_mean * Ni) / (4 * pi * pow(R_t, 2)) );
}

double f_Ni(double Pi, double m_mean) // partial pressure + temporal mean molecular number 竊� molecular number
{
  return ( (4 * pi * pow(R_t, 2)) * Pi / (g * m_mean) );
}

double f_N1_crit(double N2, double N3)  // CO2 critical partial pressure 竊� CO2 molecules number
{
  double b = m2 * N2 + m3 * N3 - 4 * pi * pow(R_t, 2) * P1_crit / g;  // coefficience in quadratic equation
  double c = - 4 * pi * pow(R_t, 2) / g * P1_crit * (N2 + N3);  // constant in quadratic equation
  return ((- b + pow((b * b - 4 * m1 * c), 0.5)) / (2 * m1)); // quadratic formula
}

double f_N2_crit(double N1, double N3)  // H2O critical partial pressure 竊� H2O molecules number
{
  double b = m1 * N1 + m3 * N3 - 4 * pi * pow(R_t, 2) * P2_crit / g;
  double c = - 4 * pi * pow(R_t, 2) / g * P2_crit * (N1 + N3);
  return ((- b + pow((b * b - 4 * m2 * c), 0.5)) / (2 * m2)); // quadratic formula
}


void f_N123_check(double *N1_temp, double *N2_temp, double *N3_temp)  // check temporal molecules number >=< critical partial pressure 竊� update molecules number
{
  int i = 0;  // loop number

  P1 = f_Pi(*N1_temp, m_mean);  // temporal calculation
  P2 = f_Pi(*N2_temp, m_mean);
  P3 = f_Pi(*N3_temp, m_mean);
          printf("Partial pressure calculation is done\n");

  if (P1 > P1_crit) // check CO2 critical upper limit of partial pressure
  {
         printf("CO2 condenses\n");
      if (P2 > P2_crit) // check H2O critical upper limit of partial pressure
      {
            printf("H2O condenses\n");

        double Rate = 1;  // changing rate of N
        double N1_temp_f = *N1_temp;  // temporal molecules number
        double N2_temp_f = *N2_temp;

        for ( i = 0 ; Rate > 1e-10 ; ++i) // waiting converging
        {
          m_mean = f_m_mean(N1, N2, *N3_temp);
          *N1_temp = f_Ni(P1_crit, m_mean);
          m_mean = f_m_mean(*N1_temp, N2, *N3_temp);
          *N2_temp = f_Ni(P2_crit, m_mean);
          Rate = fmax((*N1_temp - N1_temp_f) / N1_temp_f , (*N2_temp - N2_temp_f) / N2_temp_f);
          N1_temp_f = *N1_temp;
          N2_temp_f = *N2_temp;
        }
          printf("CO2 & H2O condenses N1 & N2 calculation is done\n");
      }
      else
      {
        *N1_temp = f_N1_crit(*N2_temp, *N3_temp);
        printf("only CO2 condense N1 calculation is done\n");
      }
  }
  else
  {
      if (P2 > P2_crit) // check H2O critical upper limit of partial pressure
      {
        printf("H2O condenses\n");
        *N2_temp = f_N2_crit(*N1_temp, *N3_temp);
        printf("only H2O condense N2 calculation is done\n");
      }
  }
  printf("f_N123check : done\n");
}


void f_etazeta_ave(double *E_A_ave, double *E_p_ave, double d_p, double d_t, double v_esc, double H, double d_0, double (f_mass) (double y), double (f_Velocity) (double y)) // calculate average efficiency over V distributions
{
  double E_A_Dave = 0;  // average ﾎｷ over impactor size distribution 
  double E_p_Dave = 0;  // average ﾎｶ over impactor size distribution (ﾏ㍉pr, E_pr^Shuvalov2009, ﾎｶ)

  Ci_V = f_Vnormal(v_esc); //calculate the normalizing constant for V distribution with v_esc

  double ans_A = 0; // last integral value of ﾎｷ
  double ans_p = 0; // last integral value of ﾎｶ
  int i = 0; // integer
  size_t N = 1; // bin number
  double x = 0; // variable : impact velocity
  double const x_min = v_esc; // start point x
  double const x_max = 4.324 * v_esc; // end point x, Earth��48.42 km/s, Venus��44.79 km/s, Mars��21.75 km/s
  double dx = 0; // bin width ﾎ肺
  double F_A = 0; // ﾎｷ integral value in case of each x
  double F_p = 0; // ﾎｶ integral value in case of each x
  double const F_min = f_Velocity(x_min); // start point integral value
  double const F_max = f_Velocity(x_max); // end point integral value

  double A_A = 0; // Trapezoidal integral law value
  double A_f_A = 0; // one former step Trapezoidal integral value
  double B_A = 1; // Simpson's integral law value
  double B_f_A = 1; // one former step Simpson's integral value
  double T_A = 1;

  double A_p = 0; // Trapezoidal integral law value
  double A_f_p = 0; // one former step Trapezoidal integral value
  double B_p = 1; // Simpson's integral law value
  double B_f_p = 1; // one former step Simpson's integral value
  double T_p = 1;

   f_etazeta_Dave(&E_A_Dave, &E_p_Dave, x_min, d_p, d_t, v_esc, H, d_0, f_mass);
  double const F_min_A = f_Velocity(x_min) * E_A_Dave;
  double const F_min_p = f_Velocity(x_min) * E_p_Dave;
   f_etazeta_Dave(&E_A_Dave, &E_p_Dave, x_max, d_p, d_t, v_esc, H, d_0, f_mass);
  double const F_max_A = f_Velocity(x_max) * E_A_Dave;
  double const F_max_p = f_Velocity(x_max) * E_p_Dave;

  dx = (x_max - x_min) * 2.;    // bin width ﾎ肺 twice of the initial value N = 1 
  A_A = F_min_A + F_max_A;  // start & end point
  A_p = F_min_p + F_max_p;  // start & end point
for( N = 1 ; T_A > 1e-2 || T_p > 1e-2 ; N = N * 2)
{
  B_f_A = B_A;
  A_f_A = A_A; // on former step integral value
    A_A = 1. / 2. * A_A; // reset the integral value
  B_f_p = B_p;
  A_f_p = A_p; // one former step integral value
    A_p = 1. / 2. * A_p; // reset the integral value
    dx = dx / 2; // bin width ﾎ肺
    x = x_min - dx; // initial x(-1)
  for( i = 0 ; i < N / 2 ; i++ )
  {
    x = x + 2 * dx;   // update x
     f_etazeta_Dave(&E_A_Dave, &E_p_Dave, x, d_p, d_t, v_esc, H, d_0, f_mass);
    F_A = f_Velocity(x) * E_A_Dave;  // calculate function value
    A_A = A_A + F_A * dx;  // calculate integral value
    F_p = f_Velocity(x) * E_p_Dave;  // calculate function value
    A_p = A_p + F_p * dx;  // calculate integral value
  }
  B_A = 4. / 3. * A_A - 1. / 3. * A_f_A ;  // Simpson_2N = 4/3 S_2N - 1/3 S_N
  T_A = fabs ( ( B_A - B_f_A ) / B_f_A ); // changing rate
  B_p = 4. / 3. * A_p - 1. / 3. * A_f_p ;  // Simpson_2N = 4/3 S_2N - 1/3 S_N
  T_p = fabs ( ( B_p - B_f_p ) / B_f_p ); // changing rate
 }

 ans_A = ans_A + B_A;
 ans_p = ans_p + B_p;

 *E_A_ave = ans_A;
 *E_p_ave = ans_p;

printf("Eta & Zeta Integral calculation is done\n");

}

void f_etazeta_Dave(double *E_A_Dave, double *E_p_Dave, double V, double d_p, double d_t, double v_esc, double H, double d_0, double (f_mass) (double y)) // calculate average efficiency over D distributions
{
  double E_A = 0; // later input pointer of this into the function of ﾎｷ  
  double E_p = 0; // later input pointer of this into the function of ﾎｶ

  double ans_A = 0; // integral value ﾎｷ
  double ans_p = 0; // integral value ﾎｶ
  int i = 0; // integer
  size_t N = 1; // bin number
  double x = 0; // variable : impactor diameter D
  double dx = 0; // bin width ﾎ肺
  double F_A = 0; // integral calculation
  double F_p = 0; // integral calculation
  double A_A = 0; // Trapezoidal integral law value
  double A_f_A = 0; // one former step Trapezoidal integral law value
  double B_A = 1; // Simpson's integral law value 
  double B_f_A = 1; // one former step Simpson's integral law value 
  double T_A = 1;
  double A_p = 0;   // Trapezoidal integral law value
  double A_f_p = 0; // one former step Trapezoidal integral law value
  double B_p = 1; // Simpson's integral law value 
  double B_f_p = 1; // one former step Simpson's integral law value 
  double T_p = 1;

// function-1 : positive gradient part

  double const x1_min = 3.5; // start point
  double const x1_max = 3.6874278; // end point
   f_etazeta(&E_A, &E_p, pow(10, x1_min), V, d_p, d_t, v_esc, H, d_0);
  double const F1_min_A = f_mass(x1_min) * E_A;
  double const F1_min_p = f_mass(x1_min) * E_p;
   f_etazeta(&E_A, &E_p, pow(10, x1_max), V, d_p, d_t, v_esc, H, d_0);
  double const F1_max_A = f_mass(x1_max) * E_A;
  double const F1_max_p = f_mass(x1_max) * E_p;

  B_A = 1;
  T_A = 1;
  B_p = 1;
  T_p = 1;
  dx = (x1_max - x1_min) * 2.;    // bin width ﾎ肺 twice of the initial number N=1
  A_A = F1_min_A + F1_max_A;  // start & end point
  A_p = F1_min_p + F1_max_p;  // start & end point

for( N = 1 ; T_A > 1e-2 || T_p > 1e-2 ; N = N * 2)
{
  B_f_A = B_A;
  A_f_A = A_A; // one former step integral value
    A_A = 1. / 2. * A_A; // reset the integral value 
  B_f_p = B_p;
  A_f_p = A_p; // one former step integral value
    A_p = 1. / 2. * A_p; // reset the integral value 
    dx = dx / 2;    // bin width ﾎ肺
    x = x1_min - dx; // initial x(-1)
  for( i = 0 ; i < N / 2 ; i++ )
  {
    x = x + 2 * dx;   // update x
     f_etazeta(&E_A, &E_p, pow(10, x), V, d_p, d_t, v_esc, H, d_0); 
    F_A = f_mass(x) * E_A;  // calculate function value
    A_A = A_A + F_A * dx;  // calculate integral value
    F_p = f_mass(x) * E_p;  // calculate function value
    A_p = A_p + F_p * dx;  // calculate integral value
  }
  B_A = 4. / 3. * A_A - 1. / 3. * A_f_A ;  // Simpson_2N = 4/3 S_2N - 1/3 S_N
  T_A = fabs ( ( B_A - B_f_A ) / B_f_A ); // changing rate
  B_p = 4. / 3. * A_p - 1. / 3. * A_f_p ;  // Simpson_2N = 4/3 S_2N - 1/3 S_N
  T_p = fabs ( ( B_p - B_f_p ) / B_f_p ); // changing rate
 }

 ans_A = ans_A + B_A;
 ans_p = ans_p + B_p;

// function-2 negative gradient part

  double const x2_min = 3.6874278; 
  double const x2_max = 8; 
   f_etazeta(&E_A, &E_p, pow(10, x2_min), V, d_p, d_t, v_esc, H, d_0);
  double const F2_min_A = f_mass(x2_min) * E_A;
  double const F2_min_p = f_mass(x2_min) * E_p;
   f_etazeta(&E_A, &E_p, pow(10, x2_max), V, d_p, d_t, v_esc, H, d_0);
  double const F2_max_A = f_mass(x2_max) * E_A;
  double const F2_max_p = f_mass(x2_max) * E_p;

  B_A = 1;
  T_A = 1;
  B_p = 1;
  T_p = 1;
  dx = (x2_max - x2_min) * 2.; 
  A_A = F2_min_A + F2_max_A; 
  A_p = F2_min_p + F2_max_p; 
for( N = 1 ; T_A > 1e-2 || T_p > 1e-2 ; N = N * 2)
{
  B_f_A = B_A;
  A_f_A = A_A;
    A_A = 1. / 2. * A_A; 
  B_f_p = B_p;
  A_f_p = A_p; 
    A_p = 1. / 2. * A_p; 
    dx = dx / 2; 
    x = x1_min - dx; 
  for( i = 0 ; i < N / 2 ; i++ )
  {
    x = x + 2 * dx; 
     f_etazeta(&E_A, &E_p, pow(10, x), V, d_p, d_t, v_esc, H, d_0); 
    F_A = f_mass(x) * E_A; // calculating the value of function 
    A_A = A_A + F_A * dx;  // calculating integral term
    F_p = f_mass(x) * E_p; // calculating the value of function 
    A_p = A_p + F_p * dx;  // calculating integral term
  }
  B_A = 4. / 3. * A_A - 1. / 3. * A_f_A ; // Simpson_2N = 4/3 S_2N - 1/3 S_N
  T_A = fabs ( ( B_A - B_f_A ) / B_f_A ); // changing rate
  B_p = 4. / 3. * A_p - 1. / 3. * A_f_p ; // Simpson_2N = 4/3 S_2N - 1/3 S_N
  T_p = fabs ( ( B_p - B_f_p ) / B_f_p ); // changing rate
 }

 ans_A = ans_A + B_A;
 ans_p = ans_p + B_p;

 *E_A_Dave = ans_A;
 *E_p_Dave = ans_p;

}


  double f_mass3(double y)  // mass distribution function dN/dD 竏� D^-3 (y = logD: 3.5 ~ 8)
  {
    return ( pow(10, y) / 43428071.353536 );  // proportional to 10^(logD)
  }  


  double f_mass2(double y)  // dN/dD 竏� D^-2
  {
    return ( pow(10, (y * 2)) / 2171472233499702.000000 );  
  }  

double f_Vnormal(double v_esc){
  double ans = 0; // integrated value = normalizing constant
  int i = 0;    // integer
  size_t N = 1;    // bin number 
  double x = 0;   // variable : impact velocity V
  double const x_min = v_esc;    // start poing x
  double const x_max = 4.324 * v_esc;    // end point x
  double dx = 0;    // bin width ﾎ肺
  double F = 0;   // integral value in case of each x
  double F_min = f_Velocity_original(x_min); // start point integral value
  double F_max = f_Velocity_original(x_max); //end point integral value
  double A = 0;   // Trapezoidal integral law value
  double A_f = 0; // one former step Trapezoidal integral value
  double B = 0; // Simpson's integral law value
  double B_f = 1; // one former step Simpson's integral value
  double T = 1;

  B = 1;
  T = 1;
  dx = (x_max - x_min) * 2.;  // calculate next bin width (n) * 2
  A = F_min + F_max;  // start & end points
for( N = 1 ; T > 1e-3 ; N = N * 2)
{
  B_f = B;
  A_f = A; // former integral value
    A = 1. / 2. * A; // reset the integral value
    dx = dx / 2;    // bin width ﾎ肺
    x = x_min - dx; // initial x(-1)
  for( i = 0 ; i < N / 2 ; i++ )
  {
    x = x + 2 * dx;   // prepare x
    F = f_Velocity_original(x);  // f(x)
    A = A + F * dx;  // calculate integration
  }
  B = 4. / 3. * A - 1. / 3. * A_f ;  // Simpson_2N = 4/3 S_2N - 1/3 S_N
  T = fabs ( ( B - B_f ) / B_f ); // changing rate of integral value
}

 ans = ans + B;
 return ans;
}

  double f_Velocity_original(double y) // impact velocity raylie distribution : feeding zone
  {
    double Vc = pow(y,2) - pow(v_esc, 2);
    return ( pow(Vc, 0.5) * exp(- Vc / (2 * pow(v_esc, 2 ))) ); // Earth縲11.2 km/s(escaping vilocity), v_max = 48.42 km/s
  }
  double f_Velocity(double y) // impact velocity raylie distribution : feeding zone
  {
    double Vc = pow(y,2) - pow(v_esc, 2);
    return ( pow(Vc, 0.5) * exp(- Vc / (2 * pow(v_esc, 2 ))) / Ci_V); // Earth縲11.2 km/s(escaping vilocity), v_max = 48.42 km/s
  }


void f_etazeta(double *E_A, double *E_p, double D, double V, double d_p, double d_t, double v_esc, double H, double d_0)  // Svetsov2000's ﾎｷ atmospheric erosion efficiency
{

// calculate the parts with input parameters
  double f_D = (double) pow(D / H, 3);    // D,H term in ﾎｷ and ﾎｶ
  double const f_d = d_p * d_t / ( d_0 * (d_t + d_p));   // ﾏ� term in ﾎｷ and ﾎｶ
  double const f_V = (double) pow( V / v_esc , 2) - 1;   // V term in non-dimensional parameter ﾎｾ
  double const f_u = d_t / d_p * V / v_esc;   // constant term in the equation between impactor escape efficiency ﾏ㍉pr and ﾎｾ
  double const E_p2 = 0.07 * f_u;  // 2nd candidate of the impactor vapor escape efficiency, inside of the {min}
  double const F_1 = 256. / 693. ;  // constant integral value x = 1
  double const f_alpha = 5. + pi - 4. / 3. ;// adjustment term fﾎｱ 

// definition
  double x = 0;  // non-dimensional parameters ﾎｾ, ﾏ�
  double lgx = 0;  // log x
  double Vi  = 0; // impact velocity at the ground of the planet
  double Vplume = 0; // expansion speed of vapor plume
  double v_V = 0;
  double F_v_V = 0;
  double f_int = 0;   // integral term in Svetsov ﾎｷ
  double E_p1 = 0;  // 1st candidate of the impactor vapor escape efficiency, inside of the {min}
  
    if (f_V < 0) // impact velocity < escape velocity
    {
      *E_p = 0;
    }
    else
    {
      x = f_D * f_d * f_V;  // non-dimensional parameter by Shuvalov2009_eq(2),deNiem2012_eq(1)
      lgx = (double) log10(x);
      Vi = V * exp(- d_0 / d_p * (( H / D ) + ( 2. * pow(H, 2) * pow(d_0 / d_p, 0.5)) / ( 3. / 4. * pow(D, 2))
                                   + ( pow(H, 3) * d_0) / ( pow(D, 3) / 8. * d_p)) );
      Vplume = pow(26. , 0.5) * Vi;
      if (Vplume < v_esc)
      {
        *E_A = 0;
      }
      else
      {
        v_V = v_esc / Vplume;  // calculating X 
        F_v_V = - 1. / 11. * pow(v_V, 11) + 5. / 9. * pow(v_V, 9) - 10. / 7. * pow(v_V, 7) + 2. * pow(v_V, 5)
                - 5. / 3. * pow(v_V, 3) + v_V;  // F(X) part
        f_int = (F_1 - F_v_V) / F_1;  // integral part
        *E_A = 3. / 4.  * d_0 / d_p * ( 2. * H / D +  16. * pow(H, 2) * pow(d_0 / d_p, 0.5) / ( 3. * pow(D, 2)) 
                                        + 16. * pow(H, 3) * d_0 / ( pow(D, 3) * d_p) ) * f_int * f_alpha;   // ﾎｷ(Svetsov)
      }

      if (*E_A < 0)
      {
        *E_A = 0;
      }

      E_p1 = 0.035 * f_u * (lgx - 1);
      *E_p = (double) fmin( E_p1, E_p2 );   // impactor vapor escape efficiency Shuvalov2009_eq(6)

      if(*E_p > 1)
      {
        *E_p = 1;
      }
      if (*E_p < 0)
      {
        *E_p = 0;
      }

      }
}

void f_eta_GI(double *E_A_GI, double V_imp, double M_imp, double v_esc, double M_t) // Schlichting 2015 Atmospheric loss fraction by Ginant Impacts
{
      d = (V_imp * M_imp) / (v_esc * M_t);
      *E_A_GI = 0.4 * d + 1.4 * pow(d, 2) - 0.8 * pow(d, 3);
}

// Henry's law C(CO2), N2(N2)
void ASMpartitioning_H(double R_ms, double Mi_t, double *Mi_a, double *Mi_s, double *Mi_m, double Si, double Di, double mn_a, double mn_i) 
{
      a = (4 * pi * pow(R_t, 2) * mn_i) / (g * M_s * Si * mn_a);
      b = Di * R_ms;
//      printf("a/s = %e, m/s = %e\n", a, b);

      *Mi_a = Mi_t * (a / (a + 1 + b));
      *Mi_s = Mi_t * (1 / (a + 1 + b));
      *Mi_m = Mi_t * (b / (a + 1 + b));
}

// not following Henry's law H(H2O)
void ASMpartitioning_notH(double R_ms, double Mi_t, double *Mi_a, double *Mi_s, double *Mi_m, double Si, double Di, double mn_a, double mn_i)
{
      a1 = (4 * pi * pow(R_t, 2) * mn_i) / (g * mn_a); 
      a2 = sqrt(a1) / (M_s * Si); // a'
      b = Di * R_ms;
      c = pow((b + 1), 2) + 4 * pow(a2, 2) * Mi_t;
//      printf("a1 =%e, a' = %e, b = %e, c = %e\n", a1, a2, b, c);

      *Mi_s = (sqrt(c) - b - 1) / (2 * pow(a2, 2));
      *Mi_a = (pow(a2, 2) * pow(*Mi_s, 2));
      *Mi_m = b * *Mi_s;
}

// Henry's law C(CO2), N2(N2)
void ASpartitioning_H(double Mi_t_as, double *Mi_a, double *Mi_s, double Si, double mn_a, double mn_i) 
{
      a = (4 * pi * pow(R_t, 2) * mn_i) / (g * M_s * Si * mn_a);
//      printf("a/s = %e, m/s = %e\n", a, b);

      *Mi_a = Mi_t_as * (a / (a + 1));
      *Mi_s = Mi_t_as * (1 / (a + 1));
}

// not following Henry's law H(H2O)
void ASpartitioning_notH(double Mi_t_as, double *Mi_a, double *Mi_s, double Si, double mn_a, double mn_i)
{
      a1 = (4 * pi * pow(R_t, 2) * mn_i) / (g * mn_a); 
      a2 = sqrt(a1) / (M_s * Si); // a'
      c = 1 + 4 * pow(a2, 2) * Mi_t_as;
//      printf("a1 =%e, a' = %e, b = %e, c = %e\n", a1, a2, b, c);

      *Mi_s = (sqrt(c) - 1) / (2 * pow(a2, 2)); // solution of Quadratic equation
      *Mi_a = (pow(a2, 2) * pow(*Mi_s, 2));
}

// closed silicate-metal equilibration
void SMpartitioning(double R_ms, double Mi_t_sm, double *Mi_s_GI, double *Mi_m, double Di) 
{
      b = Di * R_ms;
      *Mi_s_GI = Mi_t_sm * (1 / (1 + b));
      *Mi_m = Mi_t_sm * (b / (1 + b));
}

// Moore 1998's model
void ASMpartitioning_H_Moore(double Mh_t, double *Mh_a, double *Mh_s, double *Mh_m, double R_ms, double Dh, double R_t, double m_mean)
{
  int i = 0;
  double x_min = 1e5;
  double x_max = Mh_t / (1 + (R_ms * Dh)) - 1e5;
  double x_mid = (x_min + x_max) / 2;
  double f_x_min = 0.1;
  double f_x_max = 100;
  double f_x_mid = 50;

  for( i = 1 ; fabs((x_max - x_min) / x_min) > 1e-4 ; i++ ) // stop when range is small enough
  {
    f_x_min = f_Mh_s_Moore_ASM(Mh_t, x_min, R_ms, Dh, R_t, m_mean);
    f_x_max = f_Mh_s_Moore_ASM(Mh_t, x_max, R_ms, Dh, R_t, m_mean);
    f_x_mid = f_Mh_s_Moore_ASM(Mh_t, x_mid, R_ms, Dh, R_t, m_mean);

    if ((f_x_min * f_x_mid) < 0)
    {
      x_max = x_mid;
    }
    else{
      x_min = x_mid;
    }

    x_mid = (x_min + x_max) / 2;

  }

      *Mh_s = x_mid;
      *Mh_m = Dh * R_ms * *Mh_s;
      *Mh_a = Mh_t - *Mh_s * (1 + (R_ms * Dh));
}


void ASpartitioning_H_Moore(double Mh_t_as, double *Mh_a, double *Mh_s, double R_t, double m_mean)
{
  int i = 0;
  double x_min = 1e5;
  double x_max = Mh_t_as - 1e5;
  double x_mid = (x_min + x_max) / 2;
  double f_x_min = 0.1;
  double f_x_max = 100;
  double f_x_mid = 50;

  for( i = 1 ; fabs((x_max - x_min) / x_min) > 1e-4 ; i++ ) // stop when range is small enough
  {
    f_x_min = f_Mh_s_Moore_AS(Mh_t_as, x_min, R_t, m_mean);
    f_x_max = f_Mh_s_Moore_AS(Mh_t_as, x_max, R_t, m_mean);
    f_x_mid = f_Mh_s_Moore_AS(Mh_t_as, x_mid, R_t, m_mean);

    if ((f_x_min * f_x_mid) < 0)
    {
      x_max = x_mid;
    }
    else{
      x_min = x_mid;
    }

    x_mid = (x_min + x_max) / 2;

  }

      *Mh_s = x_mid;
      *Mh_a = Mh_t_as - *Mh_s;
}

double f_Mh_s_Moore_ASM(double Mh_t, double Mh_s, double R_ms, double Dh, double R_t, double m_mean)
{
  double const a_Moore = 2565;
  double const b_Al = -1.997;
  double const b_Fe = -0.9275;
  double const b_Na = 2.736;
  double const X_Al = 13.7e-2;
  double const X_Fe = 12.4e-2;
  double const X_Na = 2.68e-2;
  double const c_Moore = 1.171;
  double const d_Moore = -14.21;
  double const r_CtoX = 0.309;
  double Pi_Mi = (g * m_mean) / (4 * pi * pow(R_t,2) * (mn_H / N_A) * 2); // ~4e-15
  double Mh_a_Moore = Mh_t - Mh_s * (1 + (R_ms * Dh)); // 1e23

  double P_H2O = Pi_Mi * Mh_a_Moore; // 4e8
  double bXsum = b_Al * X_Al + b_Fe * X_Fe + b_Na * X_Na; // -0.3 
  double f_Mh_s = c_Moore * log(P_H2O * 1e-6) + bXsum / T * (P_H2O + P1 + P3) * 1e-6 - 2 * log((mn_H2O * Mh_s) / (r_CtoX * M_s * 2 * mn_H)) + a_Moore / T + d_Moore; // calculate P1(CO2), P3(N2)
 
  return f_Mh_s;
}

double f_Mh_s_Moore_AS(double Mh_t_as, double Mh_s, double R_t, double m_mean)
{
  double const a_Moore = 2565;
  double const b_Al = -1.997;
  double const b_Fe = -0.9275;
  double const b_Na = 2.736;
  double const X_Al = 13.7e-2;
  double const X_Fe = 12.4e-2;
  double const X_Na = 2.68e-2;
  double const c_Moore = 1.171;
  double const d_Moore = -14.21;
  double const r_CtoX = 0.309;
  double Pi_Mi = (g * m_mean) / (4 * pi * pow(R_t,2) * (mn_H / N_A) * 2); 
  double Mh_a_Moore = Mh_t_as - Mh_s; 
  double P_H2O = Pi_Mi * Mh_a_Moore; 
  double bXsum = b_Al * X_Al + b_Fe * X_Fe + b_Na * X_Na; // -0.3 
  double f_Mh_s = c_Moore * log(P_H2O * 1e-6) + bXsum / T * (P_H2O + P1 + P3) * 1e-6 - 2 * log((mn_H2O * Mh_s) / (r_CtoX * M_s * 2 * mn_H)) + a_Moore / T + d_Moore; // calculate P1(CO2), P3(N2) in this step formerly :-inf
  
  return f_Mh_s;
}
