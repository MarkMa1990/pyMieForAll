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
 *                  Python Wrapper for MieScatterForAll Version 0.1
 *
 * Author: Hongfeng Ma, 2015
 *
 *----------------------------------------------------------------------
 *
 *
 */

#include <iostream>
#include <complex>
#include "../MieScatterForAll_V0_1.hpp"


using namespace std;

class PyMieScatterForAll
{

public:
    PyMieScatterForAll();
    ~PyMieScatterForAll();

    void mieCal(double lmb, double radius, double fv, double n_i_real, double n_i_imag, double n_medium, double *data_out, int Nx);

private:

    MieScatterForAll Mie_test;
    complex<double> i_i;
    complex<double> n_i_cal;

};

