/*Copyright 2014, Grant Brown*/

#include <CovariateMatrix.hpp>
#include <iostream>
#include <sstream>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    int CovariateMatrix::genFromText(std::string filename)
    {
        try
        {
            readDataFile(filename.c_str());
        }
        catch(int e)
        {
            cout << "An exception ocurred building CovariateMatrix"
                << " from a text file. Error code " << e << endl;
        }
        return 0;
    }

    int CovariateMatrix::genFromDataStream(double *indata,
                                           int *inrow, 
                                           int *incol)
    {
        double* data = new double[(*incol)*(*inrow)];
        double* nrow; (*nrow) = (*inrow);
        double* ncol; (*ncol) = (*incol);

        int i; 
        for (i = 0; i < (*incol)*(*inrow); i++)
        {
            data[i] = indata[i];
        }
    }

    int CovariateMatrix::readDataFile(const char fn[])
    { 
      /* first row of input file contains two ints: num rows and num cols */
      /* rest of input file is data vals in row major order */
      int status;
      int size;
      int *data;
      int m,n,i;
      FILE* fp;
      fp=fopen(fn,"r");
      if (!fp)
      {
          printf("Error opening file\n");
      }
      status = fseek(fp, 0, SEEK_END);
      if (status != 0)
      {
          printf("Error seeking end of file\n");
          throw(-1);
      }
      size = ftell(fp);
      if (size<0)
      {
          printf("Error getting file position\n");
          throw(-1);
      }
      rewind(fp);

      fscanf(fp,"%d %d",&m, &n);
      printf("Rows: %d, Columns:  %d\n", m, n);
      data = new int[m*n];
      for (i=0;i<m*n;i++)
      {
          fscanf(fp,"%d",&data[i]);
      }
      fclose(fp);

      nrow = new int;
      ncol = new  int;
      varnames = new std::vector<std::string>();

      *nrow = m; *ncol = m; 

      std::stringstream vname;
      for (i = 0; i < *ncol; i++)
      {
          vname.str("V");
          vname << i;
          varnames -> push_back((vname.str()));
      }
      return(1);
    }


    CovariateMatrix::~CovariateMatrix()
    {
        delete[] data;
        delete[] nrow;
        delete[] ncol;
    }

}
