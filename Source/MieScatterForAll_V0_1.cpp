/* ---------------------------------------------------------------------
 * The MIT License
 *
 * Copyright (c) <2019> <Hongfeng Ma>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * ---------------------------------------------------------------------
 *                  MieScatterForAll Version 0.1
 *
 * Author: Hongfeng Ma, 2015
 *
 *----------------------------------------------------------------------
 * 
 * Example,
 * 
 * MieScatterForAll Mie_test;
 * Mie_test(Wavelength_in_vaccum, Particel_Radius, volume_fraction, n_index_Ag, n_index_medium);
 * ---------------
 * 
 * Input:
 *
 * Wavelength in vaccum      // in unit,[um], e.g., 532nm, here use 0.532
 * Particle Radius           // in unit,[um], e.g., 10nm,  here use 0.010
 * Volume fraction           // in unit,[1],  e.g., equals to 1 if there is only one size of particle
 * Optical index of Particle // in unit,[1], e.g., datatype complex<double>, e.g., n_Ag=0.050+3.021j
 * Optical index of Medium   // in unit,[1], e.g., datatype double, should be real.
 *
 * -----------
 *  
 * Output:
 *
 * Qabs
 *
 *----------------------------------------------------------------------
 *
 *
 */

#include <iostream>
#include <complex>
#include <cmath>
#include <complex_bessel.h>
#include <exception>

#include "MieScatterForAll_V0_1.hpp"




MieScatterForAll::MieScatterForAll()
{
  Nmax = 0;
  pi = 2 * asin(1);    //define PI
  imag_i = sqrt(-1);   //define imagnary number i
}

double MieScatterForAll::Mie_abcd (std::complex<double> m, double x, std::complex<double>* p_an, std::complex<double>* p_bn, std::complex<double>* p_cn, std::complex<double>* p_dn)
{
  //Nmax have already been defined.
	  
  //
  //calculate [xjn(x)]', as div_xjnx

  std::complex<double> *div_xJnx   = new std::complex<double>[Nmax];
  std::complex<double> *div_mxJnmx = new std::complex<double>[Nmax];
  std::complex<double> *div_xHnx   = new std::complex<double>[Nmax];

  double n_loop;
  
  std::complex<double> mx;

  mx = m * x;

  for ( int num_i = 0; num_i < Nmax; num_i++)
  {
	n_loop = num_i + 1;
    div_xJnx[num_i] = x * sp_bessel::sph_besselJ(n_loop-1,x) - n_loop * sp_bessel::sph_besselJ(n_loop,x);
	div_mxJnmx[num_i] = mx * sp_bessel::sph_besselJ(n_loop-1,mx) - n_loop * sp_bessel::sph_besselJ(n_loop,mx);
	div_xHnx[num_i] = x * sp_bessel::sph_hankelH1(n_loop-1,x) - n_loop * sp_bessel::sph_hankelH1(n_loop,x);

	p_an[num_i] = (m * m * sp_bessel::sph_besselJ(n_loop,mx) * div_xJnx[num_i]
			  -sp_bessel::sph_besselJ(n_loop,x) * div_mxJnmx[num_i] )
		      /
		      (m * m * sp_bessel::sph_besselJ(n_loop,mx) * div_xHnx[num_i]
			  -sp_bessel::sph_hankelH1(n_loop,x) * div_mxJnmx[num_i]);

	p_bn[num_i] = (sp_bessel::sph_besselJ(n_loop,mx) * div_xJnx[num_i]
			  -sp_bessel::sph_besselJ(n_loop,x) * div_mxJnmx[num_i] )
		      /
		      (sp_bessel::sph_besselJ(n_loop,mx) * div_xHnx[num_i]
			  -sp_bessel::sph_hankelH1(n_loop,x) * div_mxJnmx[num_i]);

	p_cn[num_i] = (sp_bessel::sph_besselJ(n_loop,x) * div_xHnx[num_i]
			  -sp_bessel::sph_hankelH1(n_loop,x) * div_xJnx[num_i] )
		      /
		      (sp_bessel::sph_besselJ(n_loop,mx) * div_xHnx[num_i]
			  -sp_bessel::sph_hankelH1(n_loop,x) * div_mxJnmx[num_i]);

	p_dn[num_i] = (m * sp_bessel::sph_besselJ(n_loop,x) * div_xHnx[num_i]
			  -m * sp_bessel::sph_hankelH1(n_loop,x) * div_xJnx[num_i] )
		      /
		      (m * m * sp_bessel::sph_besselJ(n_loop,mx) * div_xHnx[num_i]
			  -sp_bessel::sph_hankelH1(n_loop,x) * div_mxJnmx[num_i]);


  }

 //delete arrays

  delete [] div_xJnx;
  delete [] div_mxJnmx;
  delete [] div_xHnx; 



  return 1;
}


// *data_out
// data_out[0] = Qext
// data_out[1] = Qabs
// data_out[2] = Qsca
// polarization = 2*pi*r^3*eps0/q/q/q*sum(j(2l+1)(a_l+b_l));
// alpha_polarization = 2/q/q/q*sum(j(2l+1)(a_l+b_l)) to avoid unit conversion
double MieScatterForAll::MieCal (double lmb, double radius, double fv, std::complex<double> n_particle, double n_medium, double *data_out)
{
  double Vsphere;    //volume of sphere
  double rho;        //concentration of spheres, #um^-3

  //for Mie_abcd calculation
  std::complex<double> m;  //ratio of refraction indices
  double x;  //ratio circumference / wavelength in medium


  Vsphere = 4.0 / 3.0 * pi * radius * radius * radius;
  rho     = fv / Vsphere;
  
  m = n_particle / n_medium;
  x = radius * 2.0 * pi / lmb * n_medium;

  Nmax = round ( 2.0 + x + 4.0*pow(x,1.0/3.0));

  try{

      if (Nmax >=200 || Nmax<0)
      {
//          std::cout << "!____________ERROR___________!"<< std::endl << "________Nmax = " << Nmax << " is too larger or INVALID_________" << std::endl;
          throw "Nmax too large or invalid value";
      }
    
      std::complex<double> *an = new std::complex<double>[Nmax];
      std::complex<double> *bn = new std::complex<double>[Nmax];
      std::complex<double> *cn = new std::complex<double>[Nmax];
      std::complex<double> *dn = new std::complex<double>[Nmax];
      
    
    
      //call Mie_abcd
      //data are stored in an, bn, cn and dn
      Mie_abcd(m, x, an, bn, cn, dn);
    
      double *Dn = new double[Nmax];
      double *En = new double[Nmax];
    
    
    
      
      double real_an,imag_an,real_bn,imag_bn;
    
      for (int num_i=0;num_i<Nmax;num_i++)
      {
    	real_an = std::real(an[num_i]);
    	imag_an = std::imag(an[num_i]);
    	real_bn = std::real(bn[num_i]);
    	imag_bn = std::imag(bn[num_i]);
    
        Dn[num_i] = (2*(num_i+1) +1 ) * ( real_an + real_bn );
    	En[num_i] = (2*(num_i+1) +1 ) * ( real_an * real_an + imag_an * imag_an + real_bn * real_bn + imag_bn * imag_bn);
    
    
      }
    
      double qsca;
      double qext;
      double qabs;
    
    
      qext = 2.0 * sum_func (Dn, Nmax) / x / x;
      qsca = 2.0 * sum_func (En, Nmax) / x / x;
    
      qabs = qext - qsca;
    
    
      //delete arrays
      delete [] an;
      delete [] bn;
      delete [] cn;
      delete [] dn;
      
      delete [] Dn;
      delete [] En;
    
    
      data_out[0] = qext;
      data_out[1] = qabs;
      data_out[2] = qsca;
    
    
    ////////
//      return qext;

  }
  catch(const char *msg)
    {
      std::terminate();
    }

    return 0;

}

  
double MieScatterForAll::sum_func(double* array, int Num)
{
  double sum_comp_result = 0;
  for (int num_i=0;num_i<Num;num_i++)
  {
	sum_comp_result += array[num_i];

  }

  return sum_comp_result;

}



