/*Copyright 2014, Grant Brown*/


#include <iostream>
#include <fstream>

#ifndef COMPARTMENTALMODELMATRIX 
#define COMPARTMENTALMODELMATRIX
#include <CompartmentalModelMatrix.hpp>
#endif

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    int CompartmentalModelMatrix::genFromText(std::string filename)
    {
        try
        {
            readDataFile(filename.c_str());
        }
        catch(int e)
        {
            cout << "An exception ocurred building CompartmentalModelMatrix"
                << " from a text file. Error code " << e << endl;
        }
        return 0;
    }

    int CompartmentalModelMatrix::genFromDataStream(int *indata,
                                                    int *inrow, 
                                                    int *incol)
    {
        int numToAlloc = (*incol)*(*inrow);
        this -> data = new int[numToAlloc];
        this -> nrow = new int;
        this -> ncol = new int;
        (*nrow) = (*inrow);
        (*ncol) = (*incol);
        int i; 
        for (i = 0; i < (*incol)*(*inrow); i++)
        {
            data[i] = indata[i];
        }
    }

    int CompartmentalModelMatrix::createEmptyCompartment(int *inrow, int *incol)
    {
        this -> nrow = new int;
        this -> ncol = new int;
        *nrow = *inrow;
        *ncol = *incol;

        data = new int[(*nrow) * (*ncol)];
        int i;
        for (i = 0; i < ((*nrow)*(*ncol)); i ++ )
        {
            data[i] = 0;
        }
    }


    int CompartmentalModelMatrix::readDataFile(const char fn[])
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
      this -> data = new int[m*n];
      for (i=0;i<m*n;i++)
      {
          fscanf(fp,"%d",&data[i]);
      }
      fclose(fp);

      this -> nrow = new int;
      this -> ncol = new  int;
      *nrow = m; *ncol = m; 

      return(1);
    }


    CompartmentalModelMatrix::~CompartmentalModelMatrix()
    {
        delete[] data;
        delete[] nrow;
        delete[] ncol;
    }

}
