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


class MieScatterForAll
{
  public:
  MieScatterForAll ();

  double Mie_abcd(std::complex<double> m, double x, std::complex<double>* p_an, std::complex<double>* p_bn, std::complex<double>* p_cn, std::complex<double>* p_dn);
  double MieCal (double lmb, double radius, double fv, std::complex<double> n_particle, double n_medium, double *data_out );

  private:

  double sum_func(double * array, int Num);


  int Nmax;
  double pi;
  std::complex<double> imag_i;
  
};


