#include <iostream>
#include <fstream>


using namespace std;

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

ofstream OutFile("part1b.dat");

int main(int argc, char** argv)
{


int d=1;
double g=1;
double e[100][100] = {0};
int i=0; int j=0;
double x[200]={0};
double y[200]={0};

int n=2,m=2;
for (g=-1; g<=1; g+=0.1)
{
double data[4]={-g, -g, -g, 2*d-g};

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

  double min=e[i][0];
  for(j=0;j<2;j++)
  if(e[i][j]<=min)
    min = e[i][j];

  //cout  << min << endl;

  x[i]=g;
  y[i]=min+g;

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



