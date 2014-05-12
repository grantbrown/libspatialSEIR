#include <RandomNumberProvider.hpp>


namespace SpatialSEIR
{

    

    RandomNumberProvider::RandomNumberProvider(unsigned int seed)
    {
        generator.seed(seed);
        int memosize = 1000000;
        logFactorialMemo = new double[memosize];
        maxFactorial = new int;
        *maxFactorial = memosize;
        int i;
        logFactorialMemo[0] = 0;

        for (i = 1; i < memosize; i++)
        {
            logFactorialMemo[i] = std::log(i) + logFactorialMemo[i-1];
        }
    }
    double RandomNumberProvider::factorial(int k)
    {
        if (k <= *maxFactorial)
        {
            return(logFactorialMemo[k]);
        }
        // If not memoized, default to Stirling's approximation
        return(k * std::log(k) - k);
    }

    double RandomNumberProvider::choose(int n, int k)
    {
        // Can we get rid of these bounds checking
        // parts if conditions are satisfied on N?
        // It's only an issue because this stuff ends
        // up in the innermost loop.
        if (n > *maxFactorial)
        {
            // Out of memo
            return(factorial(n) - 
                    factorial(k) - 
                    factorial(n-k));
        } 
        return(logFactorialMemo[n] - 
                logFactorialMemo[k] - 
                logFactorialMemo[n-k]);
    }

    double RandomNumberProvider::choosePartial(int n, int k)
    {
        // Can we get rid of these bounds checking
        // parts if conditions are satisfied on N?
        // It's only an issue because this stuff ends
        // up in the innermost loop.
        if (n > *maxFactorial)
        {
            // Out of memo
            return(factorial(n) - 
                    factorial(n-k));
        } 
        return(logFactorialMemo[n] - 
                logFactorialMemo[n-k]);
    }



    double RandomNumberProvider::uniform()
    {
        return(unidist(generator));
    }
    int RandomNumberProvider::uniform_int()
    {
        return(unidist_int(generator));
    }
    int RandomNumberProvider::uniform_int(int a, int b)
    {
        boost::random::uniform_int_distribution<> udist(a,b);
        return(udist(generator));
    }
    int RandomNumberProvider::poisson(int mu)
    {
        boost::random::poisson_distribution<> pdist(mu);
        return(pdist(generator));
    }
    double RandomNumberProvider::normal(double mu, double sd)
    {
        boost::random::normal_distribution<> ndist(mu, sd);
        return(ndist(generator));
    }
    double RandomNumberProvider::dnorm(double x, double mu, double sd)
    {
        boost::math::normal_distribution<> ndist(mu, sd);
        return(std::log(pdf(ndist, x))); 
        // No way to get the log density directly with boost?
        // I feel like this should be either an optional flag given to 
        // the pdf non-member pdf getter or perhaps an additional non-member
        // getter, as this is a really common need in statistical computing.  
    }
    double RandomNumberProvider::dpois(int x, double mu)
    {
        double out = x*std::log(mu) - mu;
        int i;
        for (i = 1; i <= x; i++)
        {
            out -= std::log(i);
        }
        return(out);
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
        if (a < 0 || b < 0)
        {
            std::cerr << "Invalid (negative) Paramter Values to RandomNumberProvider::beta.\n";
            std::cout << "(a,b): (" << a << ", " << b << ")\n";
            throw(-1);
        }
        double v1 = gamma(a); 
        double v2 = gamma(b);
        return(v1/(v1 + v2));
    }
    RandomNumberProvider::~RandomNumberProvider()
    {
        delete maxFactorial;
        delete[] logFactorialMemo;
    }
}
