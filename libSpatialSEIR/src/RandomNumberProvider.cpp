#include <RandomNumberProvider.hpp>


namespace SpatialSEIR
{

    

    RandomNumberProvider::RandomNumberProvider(unsigned int seed)
    {
        generator.seed(seed);
    }
    double RandomNumberProvider::uniform()
    {
        return(unidist(generator));
    }
    double* RandomNumberProvider::uniform(int n)
    {
        double* output = new double[n];
        int i;
        for (i = 0; i < n; i++)
        {
            output[i] = unidist(generator);
        }
        return(output);
    }
    double* RandomNumberProvider::uniform(int n, double* output)
    {
        int i;
        for (i = 0; i < n; i++)
        {
            output[i] = unidist(generator);
        }
        return(output);
    }
    double RandomNumberProvider::gamma()
    {
        return(gammadist(generator));
    }
    double RandomNumberProvider::gamma(double a)
    {
        boost::random::gamma_distribution<> gdist(a);
        return(gdist(generator));  
    }
    double* RandomNumberProvider::gamma(int n)
    {
        double* output = new double[n];
        int i;
        for (i = 0; i < n; i++)
        {
            output[i] = gammadist(generator);
        }
        return(output);
    }
    double* RandomNumberProvider::gamma(int n, double* output)
    {
        int i;
        for (i = 0; i < n; i++)
        {
            output[i] = gammadist(generator);
        }
        return(output);
    }
    double RandomNumberProvider::beta(double a, double b)
    {
        double v1 = gamma(a); 
        double v2 = gamma(b);
        return(v1/(v1 + v2));
    }
    RandomNumberProvider::~RandomNumberProvider()
    {
        //Nothing special
    }
}
