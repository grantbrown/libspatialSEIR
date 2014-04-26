#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

#ifndef RANDOM_NUMBER_PROVIDER_INC
#define RANDOM_NUMBER_PROVIDER_INC


// Include logic about boost availability needs 
// to go here. The library should try to find boost, 
// but fall back on R. Expandin the CMake scripts
// should help with this by using the appropriate
// FIND_BOOST macros and setting appropriate 
// preprocessor directives. 

#include<boost/random/mersenne_twister.hpp>
#include<boost/random/uniform_real.hpp>
#include<boost/random/uniform_int_distribution.hpp>
#include<boost/random/gamma_distribution.hpp>
#include<boost/random/poisson_distribution.hpp>


namespace SpatialSEIR
{
    class RandomNumberProvider
    {
        boost::random::mt19937 generator;
        boost::random::uniform_real_distribution<> unidist;
        boost::random::uniform_int_distribution<> unidist_int;
        boost::random::gamma_distribution<> gammadist;
        public:
            //Random Methods
            RandomNumberProvider(unsigned int seed);
            double uniform();   
            int uniform_int();
            int uniform_int(int a, int b);
            int poisson(int mu);
            double* uniform(int n);
            double* uniform(int n, double* output);
            double gamma();
            double gamma(double a);
            double* gamma(int n);
            double* gamma(int n, double* output);
            double beta(double a, double b);

            //Density Functions
            double dpois(int x, double mu);

            ~RandomNumberProvider();
    };
}

#endif
