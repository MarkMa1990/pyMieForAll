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
 *                  Python Wrapper for MieScatterForAll Version 0.2
 *
 * Author: Hongfeng Ma, 2015
 *
 *----------------------------------------------------------------------
 *
 *
 */

#include <iostream>
#include <complex>
#include "pyMie.hpp"


using namespace std;





PyMieScatterForAll::PyMieScatterForAll()
{
    i_i = -1;

}

PyMieScatterForAll::~PyMieScatterForAll()
{

}

void PyMieScatterForAll::mieCal(double lmb, double radius, double fv, double n_i_real, double n_i_imag, double n_medium, double *data_out, int Nx)
{
    try{

    if (Nx != 3)
    {
          throw "data_out is not the size of 3";

    }
    n_i_cal = n_i_real + n_i_imag * sqrt(i_i);

    Mie_test.MieCal(lmb, radius, fv, n_i_cal, n_medium, data_out);
    }

    catch(const char *msg)
    {
      std::terminate();
    }

}


  
