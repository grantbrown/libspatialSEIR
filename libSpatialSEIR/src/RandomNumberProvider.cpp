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

}
