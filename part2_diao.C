#include <iostream>
#include <fstream>
#include<cmath> 

using namespace std;

  int d=1;

//account for E
//ceil 向上取整 floor向下取整
  double sp_f(int p, int f, double g)
  {
  if (p>f)    return floor(p/2.0)*d;

  else        return floor(p/2.0)*d-g/2;
  }

//painring V matrix
  double pairing_v_as(int p, int q, int r,int s, double g)
  {
  if  (p > q)
     return -pairing_v_as(q,p,r,s,g);
  else if(r > s)
     return -pairing_v_as(p,q,s,r,g);    

  else if (floor(p/2.0) == floor(q/2.0) && floor(r/2.0) == floor(s/2.0) && p != q && r != s)
     {
     return -g/2;}
  else
     return 0;
  }
  

int main()
{
  ofstream OutFile("part2_g_Ec.dat");
  int N_Particle = 4;
//  cout<<"The number of the particle:"<<endl;
//  cin>>N_Particle;

  int N_Level = 8;
//  cout<<"The number of the level:"<<endl;
//  cin>>N_Level;
  
  int Particle=N_Particle;
  int Level=N_Level;

  int Fermi_Level_Number; 
  Fermi_Level_Number=N_Particle-1;

//an simple example need to change
  for(double g=-1; g<=1; g+=0.05)
  {
  //cout<<g<<endl;
  double delta=0.;

//用数组实现python里的np.zeros
  double *f_h=new double[Particle];
  for(int i= 0; i<Particle;i++)
  {
     f_h[i] = sp_f(i,Fermi_Level_Number,g);
     //cout<<i<<endl; cout<<Fermi_Level_Number<<endl; cout<<g<<endl;
     //cout<<f_h[i]<<endl; 
  }

  double *f_p=new double[Level-Particle];
  for(int a = 0; a<Level-Particle;a++)
  {
     f_p[a] = sp_f(a + Particle,Fermi_Level_Number,g);
     //cout<<f_p[a]<<endl;
  }

//建立 v_pp_hh
  
  double ****v_pp_hh = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    v_pp_hh[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      v_pp_hh[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        v_pp_hh[a][b][i]= new double [Particle];       
      }
    }
  }
  //cout<<Level-Particle<<endl; cout<<Particle<<endl;

//建立 v_hh_pp
  double ****v_hh_pp = new double ***[Particle];
  for( int i=0; i<Particle; i++)
  {
    v_hh_pp[i]= new double **[Particle];
    for(int j=0;j<Particle; j++)
    {
      v_hh_pp[i][j]= new double *[Level-Particle];
      for(int a=0; a<Level-Particle; a++)
      {
        v_hh_pp[i][j][a]= new double [Level-Particle];
      }
    }
  }


//建立 v_hh_hh
  double ****v_hh_hh = new double ***[Particle];
  for( int i=0; i<Particle; i++)
  {
    v_hh_hh[i]= new double **[Particle];
    for(int j=0;j<Particle; j++)
    {
      v_hh_hh[i][j]= new double *[Particle];
      for(int k=0; k<Particle; k++)
      {
        v_hh_hh[i][j][k]= new double [Particle];
      }
    }
  }

//建立 v_pp_pp
  double ****v_pp_pp = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    v_pp_pp[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      v_pp_pp[a][b]= new double *[Level-Particle];
      for(int c=0; c<Level-Particle; c++)
      {
        v_pp_pp[a][b][c]= new double [Level-Particle];
      }
    }
  }


// python line 58 for matrix in [v_pp_hh, v_hh_hh, v_pp_pp, v_hh_pp]
// v_pp_hh
  int p; int q; int r; int s;
  for(s=0; s<Particle;s++) 
     {for(r=0; r<Particle;r++) 
        {for(q=0; q<Level-Particle;q++)
           {for(p=0; p<Level-Particle; p++)
               {
               //cout<<p<<endl; cout<<q<<endl; cout<<r<<endl; cout<<s<<endl;
               v_pp_hh[p][q][r][s]=pairing_v_as(p,q,r,s,g); }}}}
               //if(v_pp_hh[p][q][r][s]!=0.0) cout<<v_pp_hh[p][q][r][s]<<":"<<p<<","<<q<<","<<r<<","<<s<<endl;}}}}
               

// v_hh_hh
  for(p=0; p<Particle; p++)
    {for(q=0; q<Particle;q++)
      {for(r=0; r<Particle;r++)
         {for(s=0; s<Particle;s++) 
            {v_hh_hh[p][q][r][s]=pairing_v_as(p,q,r,s,g);}}}} 
             //cout<<v_pp_hh[p][q][r][s]<<endl;}}}}
             //if(v_hh_hh[p][q][r][s]!=0.0) cout<<v_hh_hh[p][q][r][s]<<":"<<p<<","<<q<<","<<r<<","<<s<<endl;}}}}
// v_pp_pp
  for(p=0; p<Level-Particle; p++)
    {for(q=0; q<Level-Particle;q++)
      {for(r=0; r<Level-Particle;r++)
       {for(s=0; s<Level-Particle;s++) 
         {v_pp_pp[p][q][r][s]=pairing_v_as(p,q,r,s,g);}}}} 
          //cout<<v_pp_hh[p][q][r][s]<<endl;}}}}
          //if(v_pp_pp[p][q][r][s]!=0.0) cout<<v_pp_pp[p][q][r][s]<<":"<<p<<","<<q<<","<<r<<","<<s<<endl;}}}}

// v_hh_pp
  for(p=0; p<Particle; p++)
    {for(q=0; q<Particle;q++)
      {for(r=0; r<Level-Particle;r++)
       {for(s=0; s<Level-Particle;s++) 
         {v_hh_pp[p][q][r][s]=pairing_v_as(p,q,r,s,g);}}}}
         //cout<<v_hh_pp[p][q][r][s]<<endl;}}}}
         //if(v_hh_pp[p][q][r][s]!=0.0) cout<<v_hh_pp[p][q][r][s]<<":"<<p<<","<<q<<","<<r<<","<<s<<endl;}}}}


//line 65 f_sign_sum = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
//建立 f_sign_sum
  double ****f_sign_sum = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    f_sign_sum[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      f_sign_sum[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        f_sign_sum[a][b][i]= new double [Particle];
      }
    }
  }


  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                {f_sign_sum[a][b][i][j] = -f_p[a] - f_p[b] + f_h[i] + f_h[j];}}}}
                 //cout<<f_sign_sum[a][b][i][j]<<":"<<a<<","<<b<<","<<i<<","<<j<<","<<endl;}}}}


//建立 t_2  
//t_2 = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
//t_2 = v_pp_hh/f_sign_sum

  double ****t_2 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    t_2[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      t_2[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        t_2[a][b][i]= new double [Particle];
      }
    }
  }


//start the cycle do-while


  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                   {t_2[a][b][i][j]=v_pp_hh[a][b][i][j]/f_sign_sum[a][b][i][j];}}}}
                   //if(t_2[a][b][i][j]!=0) 
                   //{
                   //cout<<v_hh_pp[i][j][a][b]<<","<<f_sign_sum[a][b][i][j]<<endl;
     int Itnum=0;                //cout<<t_2[a][b][i][j]<<":"<<a<<","<<b<<","<<i<<","<<j<<endl;}}}}}
   do
  {   

//h_bar = np.zeros(shape=(no_of_states-no_of_prt, no_of_states-no_of_prt, no_of_prt, no_of_prt))
  double ****h_bar = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    h_bar[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      h_bar[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        h_bar[a][b][i]= new double [Particle];
      }
    }
  }


//matrix element construction
//diag0 We already have that
//diag1 P(ab)ft-P(ij)ft
  double ****diag1 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    diag1[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      diag1[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        diag1[a][b][i]= new double [Particle];
      }
    }

  }

  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                   {diag1[a][b][i][j]=-f_sign_sum[a][b][i][j]*t_2[a][b][i][j];}}}}
                    //cout<<diag1[a][b][i][j]<<endl;}}}}



//diag2 1/2*V*t abcd,cdij
  double ****diag2 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    diag2[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      diag2[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        diag2[a][b][i]= new double [Particle];
      }
    }

  }

//初始化 diag2=0;
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                diag2[a][b][i][j]=0;}}}


//diag2 赋值
 for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                 {for(int c=0; c<Level-Particle; c++)
                     {for(int d=0; d<Level-Particle; d++)
                         {diag2[a][b][i][j]+=(1.0/2.0)*v_pp_pp[a][b][c][d]*t_2[c][d][i][j];}}}}}}
                          //if(v_pp_pp[a][b][c][d]*t_2[c][d][i][j]!=0){
                          //cout<<"v_p:"<<v_pp_pp[a][b][c][d]<<":"<<a<<","<<b<<","<<c<<","<<d<<endl;
                          //cout<<"t_2:"<<t_2[c][d][i][j]<<":"<<c<<","<<d<<","<<i<<","<<j<<endl;}}}}}}} 
                          //cout<<diag2[a][b][i][j]<<endl;}}}}}}




//diag3 1/2*V*t klij,klab
  double ****diag3 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    diag3[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      diag3[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        diag3[a][b][i]= new double [Particle];
      }
    }

  }


//初始化 diag3=0;
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                diag3[a][b][i][j]=0;}}}

//diag3 赋值
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                {for(int k=0; k<Particle; k++)
                     {for(int l=0; l<Particle; l++)
                         {diag3[a][b][i][j]+=(1.0/2.0)*v_hh_hh[k][l][i][j]*t_2[a][b][k][l];}}}}}}
                          //if(v_hh_hh[k][l][i][j]*t_2[a][b][k][l]!=0){
                          //cout<<"v_h:"<<v_hh_hh[k][l][i][j]<<":"<<k<<","<<l<<","<<i<<","<<j<<endl;
                          //cout<<"t_2:"<<t_2[a][b][k][l]<<":"<<a<<","<<b<<","<<k<<","<<l<<endl;}}}}}}} 
                          //cout<<"diag3=:"<<diag3[a][b][i][j]<<endl;}}}}
        //cout<<diag3[0][1][0][1]<<endl;



//P(ab)*P(cd)*v*t kbcj,acik don't fit to pairing model, so we just ignore it

//diag4 P(ij)*P(ab)V*t*t ijab,klcd,acik,dblj
  double ****diag4 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    diag4[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      diag4[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        diag4[a][b][i]= new double [Particle];
      }
    }

  }


//初始化 diag4=0;
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                diag4[a][b][i][j]=0;}}}


//diag4 赋值
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                 {for(int k=0; k<Particle; k++)
                     {for(int l=0; l<Particle; l++)
                         {for(int c=0; c<Level-Particle; c++)
                             {for(int d=0; d<Level-Particle; d++)
                                 {diag4[a][b][i][j]+=(1.0/2.0)*v_hh_pp[k][l][c][d]*(t_2[a][c][i][k]*t_2[d][b][l][j]
                                                                                   -t_2[a][c][j][k]*t_2[d][b][l][i]
                                                                                   -t_2[b][c][i][k]*t_2[d][a][l][j]
                                                                                   +t_2[b][c][j][k]*t_2[d][a][l][i]);}}}}}}}}


//diag5 1/2*P(ij)*V*t*t ij,klcd,cdik,ablj
  double ****diag5 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    diag5[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      diag5[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        diag5[a][b][i]= new double [Particle];
      }
    }

  }

//初始化 diag5=0;
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                diag5[a][b][i][j]=0;}}}

//diag5 赋值
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                 {for(int k=0; k<Particle; k++)
                     {for(int l=0; l<Particle; l++)
                         {for(int c=0; c<Level-Particle; c++)
                             {for(int d=0; d<Level-Particle; d++)
                                 {diag5[a][b][i][j]+=(1.0/2.0)*v_hh_pp[k][l][c][d]*(t_2[c][d][i][k]*t_2[a][b][l][j]
                                                                                   -t_2[c][d][j][k]*t_2[a][b][l][i]);}}}}}}}}

//diag6 1/2*P(ab)*V*t*t ab,klcd,ackl,dbij
  double ****diag6 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    diag6[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      diag6[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        diag6[a][b][i]= new double [Particle];
      }
    }

  }

//diag6初始化
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                diag6[a][b][i][j]=0;}}}


//diag6 赋值
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                 {for(int k=0; k<Particle; k++)
                     {for(int l=0; l<Particle; l++)
                         {for(int c=0; c<Level-Particle; c++)
                             {for(int d=0; d<Level-Particle; d++)
                                 {diag6[a][b][i][j]+=(1.0/2.0)*v_hh_pp[k][l][c][d]*(t_2[a][c][k][l]*t_2[d][b][i][j]
                                                                                   -t_2[b][c][k][l]*t_2[d][a][i][j]);}}}}}}}}


//diag7 1/4*V*t*t  klcd,cdij,abkl
  double ****diag7 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    diag7[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      diag7[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        diag7[a][b][i]= new double [Particle];
      }
    }
  }

//diag7初始化
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                diag7[a][b][i][j]=0;}}}

//diag7 赋值
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                 {for(int k=0; k<Particle; k++)
                     {for(int l=0; l<Particle; l++)
                         {for(int c=0; c<Level-Particle; c++)
                             {for(int d=0; d<Level-Particle; d++)
                                 {diag7[a][b][i][j]+=(1.0/4.0)*v_hh_pp[k][l][c][d]*t_2[c][d][i][j]*t_2[a][b][k][l];}}}}}}}}
                                     
//H sum final result
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                  {h_bar[a][b][i][j]=v_pp_hh[a][b][i][j]+diag1[a][b][i][j]+diag2[a][b][i][j]+diag3[a][b][i][j]
                                     +diag4[a][b][i][j]+diag5[a][b][i][j]+diag6[a][b][i][j]+diag7[a][b][i][j];}}}}
                                     //cout<<v_pp_hh[a][b][i][j]<<";"<<diag1[a][b][i][j]<<endl;}}}}

//enter the iteration

//define delta_t_2 = h_bar / f_sign_sum
  double ****delta_t_2 = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    delta_t_2[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      delta_t_2[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        delta_t_2[a][b][i]= new double [Particle];
      }
    }
  }
  //cout<<111<<endl;

//define t_2_new = t_2+delta_t_2
  double ****t_2_new = new double ***[Level-Particle];
  for( int a=0; a<Level-Particle; a++)
  {
    t_2_new[a]= new double **[Level-Particle];
    for(int b=0; b<Level-Particle; b++)
    {
      t_2_new[a][b]= new double *[Particle];
      for(int i=0; i<Particle; i++)
      {
        t_2_new[a][b][i]= new double [Particle];
      }
    }
  }
  //cout<<111<<endl;

// compare of t_2_new and t_2
  delta=0.;

  for(int a=0; a<Level-Particle; a++)
  {
     for(int b=0; b<Level-Particle; b++)
     {
        for(int i=0; i<Particle; i++)
        {
            for(int j=0; j<Particle; j++)
            { 
              delta_t_2[a][b][i][j]=h_bar[a][b][i][j]/f_sign_sum[a][b][i][j];
              delta=delta+delta_t_2[a][b][i][j]*delta_t_2[a][b][i][j];
              t_2_new[a][b][i][j]=delta_t_2[a][b][i][j]+t_2[a][b][i][j];
              //cout<<delta_t_2[a][b][i][j]<<endl;
              //cout<<delta<<endl;
              //cout<<t_2_new[a][b][i][j]<<endl;
              t_2[a][b][i][j]=t_2_new[a][b][i][j]; //cout<<t_2[a][b][i][j]<<endl;


            } 
        } 
     }

  }
  
  Itnum++;
  if(Itnum>1000)
  {cout<<"The cycle is too large"<<endl;
  break;}
  
  //cout<<delta<<endl;
  
  }

  while(delta>0.000001);


//define Ec 
  double Ec=0;

//Ec=1/4*V*t
  for(int a=0; a<Level-Particle; a++)
     {for(int b=0; b<Level-Particle; b++)
        {for(int i=0; i<Particle; i++)
            {for(int j=0; j<Particle; j++)
                {Ec+=(1.0/4.0)*v_pp_hh[a][b][i][j]*t_2[a][b][i][j];}}}}

  OutFile << g<<"  "<<Ec<<endl;
  }      
  OutFile.close();   
  return 0;
}


 

  
