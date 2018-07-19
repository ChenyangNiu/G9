#include <iostream>
#include <fstream>


using namespace std;

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

ofstream OutFile("part1c.dat");

int main(int argc, char** argv)
{


int d=1;
double g=1;
double e[100][100] = {0};
int i=0; int j=0;
double x[200]={0};
double y[200]={0};

int n=6,m=6;
for (g=-1; g<=1; g+=0.05)
{
double data[36]={2*d-g, -g/2, -g/2, -g/2, -g/2, 0,
      -g/2, 4*d-g, -g/2, -g/2, 0, -g/2,
      -g/2, -g/2, 6*d-g, 0, -g/2, -g/2,
      -g/2, -g/2, 0, 6*d-g, -g/2, -g/2,
      -g/2, 0, -g/2, -g/2, 8*d-g, -g/2,
       0, -g/2, -g/2, -g/2, -g/2, 10*d-g};

  // allocate data
  char Nchar='N';
  double *eigReal=new double[n];
  double *eigImag=new double[n];
  double *vl,*vr;
  int one=1;
  int lwork=6*n;
  double *work=new double[lwork];
  int info;

  // calculate eigenvalues using the DGEEV subroutine
  dgeev_(&Nchar,&Nchar,&n,data,&n,eigReal,eigImag,
        vl,&one,vr,&one,
        work,&lwork,&info);

  // output eigenvalues to stdout
  cout << "--- Eigenvalues ---" << endl;

  e[i][0]  = eigReal[0];
  e[i][1]  = eigReal[1];
  e[i][2]  = eigReal[2];
  e[i][3]  = eigReal[3];  
  e[i][4]  = eigReal[4];
  e[i][5]  = eigReal[5];

  double min=e[i][0];
  for(j=0;j<6;j++)
  if(e[i][j]<=min)
    min = e[i][j];

  x[i]=g;
  y[i]=min-(2*d-g);

  cout  << min << endl;
  cout<<x[i]<<endl;
  cout<<y[i]<<endl;
  OutFile << x[i]<<"  "<<y[i]<<endl; 


  delete [] eigReal;
  delete [] eigImag;
  delete [] work;

  ++i;
}

 OutFile.close();  
 return 0;

}



