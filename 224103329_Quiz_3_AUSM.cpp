#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float maxi(float lam)
{ float res;
   if(lam>=0)
   {
       res=lam;
   }

   else
   {
       res=0;
   }
return(res);
}

float mini(float lam)
{ float res;
   if(lam>=0)
   {
       res=0;
   }

   else
   {
       res=lam;
   }
return(res);
}

int main()
{float L,u_r,u_l,rho_r,rho_l,p_r,p_l,g;
int N;
float x,y,z;
float dx,dt,T;

printf("enter the value of u_l\n");
scanf("%f",&u_l);
printf("enter the value of u_r\n");
scanf("%f",&u_r);
printf("enter the value of rho_l\n");
scanf("%f",&rho_l);
printf("enter the value of rho_r\n");
scanf("%f",&rho_r);
printf("enter the value of p_l\n");
scanf("%f",&p_l);
printf("enter the value of p_r\n");
scanf("%f",&p_r);
printf("enter the value of gamma\n");
scanf("%f",&g);
printf("enter the domain size\n");
scanf("%f",&L);
printf("enter the number of cells\n");
scanf("%d",&N);
printf("enter the total time\n");
scanf("%f",&T);

dx=L/N;

/*
N=400;
T=0.15;
L=1;


p_r=0.1;
p_l=1;
u_l=0;
u_r=0;
rho_l=1;
rho_r=0.125;
g=1.4;
*/

float rho_e_l,rho_e_r;
rho_e_l=(rho_l*u_l*u_l*0.5)+(p_l/(g-1));      //rho_e terms
rho_e_r=(rho_r*u_r*u_r*0.5)+(p_r/(g-1));

int i,j,k;
float u_1[N+2],u_2[N+2],u_3[N+2];

 for(i=1;i<N+1;i++)     //Initial conditions
 {

    if(i<N/2)
    {
        u_1[i]=rho_l;
        u_2[i]=rho_l*u_l;
        u_3[i]=rho_e_l;

    }

    else
    {
        u_1[i]=rho_r;
        u_2[i]=rho_r*u_r;
        u_3[i]=rho_e_r;
    }


 }

/*float p_o,p_N,rho_o,rho_N,u_o,u_N;

p_o=p_l;                         //Ghost cell data
p_N=p_r;
rho_o=rho_l;

rho_N=rho_r;
u_o=-u_l;
u_N=-u_r;*/

/*u_1[0]=u_1[1];
u_2[0]=-u_2[1];
u_3[0]=u_3[1];

u_1[N+1]=u_1[N];
u_2[N+1]=-u_2[N];
u_3[N+1]=u_3[N];
*/

float t=0;
float F_1_plus[N+2],F_2_plus[N+2],F_3_plus[N+2];
float F_1_minus[N+2],F_2_minus[N+2],F_3_minus[N+2];
float U,A,H,lambda[3];
float P,e,R;
float temp_1=0,temp_2=0,temp_3=0;
float ff_1,ff_2,ff_3;
float fb_1,fb_2,fb_3;
float max_lambda,temp_max;
float M,M_plus,M_minus;
float P_plus,P_minus;

float a,b,c;
 u_1[0]=u_1[1];
 u_2[0]=-u_2[1];
 u_3[0]=u_3[1];

 u_1[N+1]=u_1[N];
 u_2[N+1]=-u_2[N];
 u_3[N+1]=u_3[N];
 while(t<T)
 {
     for(i=0;i<N+1;i++)
     { R=u_1[i];
       P=(u_3[i]-((u_2[i]*u_2[i])/(2*u_1[i])))*(g-1);
       A=(g*P)/R;
       A=sqrt(A);
       U=u_2[i]/u_1[i];
       H=(u_3[i]+P)/u_1[i];
       M=U/A;

       if(M<-1)                           //Finding M_plus
       {
           M_plus=0;
       }

       else if(M<1)
       {
           M_plus=(M+1)*(M+1)*0.25;
       }

       else
       {
           M_plus=M;
       }

       if(M<-1)                           //Finding P_plus
       {
           P_plus=0;
       }

       else if(M<1)
       {
           P_plus=(M+1)*0.5*P;
       }

       else
       {
           P_plus=P;
       }


        F_1_plus[i]=(R*A*M_plus);                                         //Calculation of F_plus
        F_2_plus[i]=(R*U*A*M_plus)+P_plus;
        F_3_plus[i]=(R*H*A*M_plus);




       R=u_1[i+1];
       P=(u_3[i+1]-((u_2[i+1]*u_2[i+1])/(2*u_1[i+1])))*(g-1);
       A=(g*P)/R;
       A=sqrt(A);
       U=u_2[i+1]/u_1[i+1];
       H=(u_3[i+1]+P)/u_1[i+1];
       M=U/A;

       if(M<-1)                           //Finding M_minus
       {
           M_minus=M;
       }

       else if(M<1)
       {
           M_minus=(M-1)*(M-1)*(-0.25);
       }

       else
       {
           M_minus=0;
       }



       if(M<-1)                           //Finding P_minus
       {
           P_minus=P;
       }

       else if(M<1)
       {
           P_minus=(1-M)*0.5*P;
       }

       else
       {
           P_minus=0;
       }



        F_1_minus[i+1]=(R*A*M_minus);                                         //Calculation of F_minus
        F_2_minus[i+1]=(R*U*A*M_minus)+P_minus;
        F_3_minus[i+1]=(R*H*A*M_minus);



       if(i==0)
     { temp_1=F_1_minus[i+1]+F_1_plus[i];
       temp_2=F_2_minus[i+1]+F_2_plus[i];
       temp_3=F_3_minus[i+1]+F_3_plus[i];
     }

      if(i==0)
      {
          continue;
      }

   fb_1=temp_1;
   fb_2=temp_2;
   fb_3=temp_3;

   ff_1=F_1_minus[i+1]+F_1_plus[i];
   ff_2=F_2_minus[i+1]+F_2_plus[i];
   ff_3=F_3_minus[i+1]+F_3_plus[i];


       for(k=0;k<N-1;k++)
     { R=u_1[k];
       P=(u_3[k]-((u_2[k]*u_2[k])/(2*u_1[k])))*(g-1);
       A=(g*P)/R;
       A=sqrt(A);
       U=u_2[k]/u_1[k];



       lambda[0]=U-A;                     //Eigen values
       lambda[1]=U;
       lambda[2]=U+A;

       temp_max=0;
       for(j=0;j<3;j++)
     {
          if(fabs(lambda[j]) > temp_max)
           temp_max = lambda[j];

     }

     if(k==0)
     {
         max_lambda=temp_max;
     }

     if(temp_max>max_lambda)
     {
         max_lambda=temp_max;
     }
    }
      dt=(0.8*dx)/max_lambda;




     if(i!=0)
     {

      u_1[i]=u_1[i]-(((dt/dx))*(ff_1-fb_1));
      u_2[i]=u_2[i]-(((dt/dx))*(ff_2-fb_2));
      u_3[i]=u_3[i]-(((dt/dx))*(ff_3-fb_3));
     }
   temp_1=ff_1;
   temp_2=ff_2;
   temp_3=ff_3;


     }
 u_1[0]=u_1[1];
 u_2[0]=-u_2[1];
 u_3[0]=u_3[1];

 u_1[N+1]=u_1[N];
 u_2[N+1]=-u_2[N];
 u_3[N+1]=u_3[N];
 t=t+dt;
 }

 if(t>T)
 {
        for(i=0;i<N+1;i++)
     { R=u_1[i];
       P=(u_3[i]-((u_2[i]*u_2[i])/(2*u_1[i])))*(g-1);
       A=(g*P)/R;
       A=sqrt(A);
       U=u_2[i]/u_1[i];
       H=(u_3[i]+P)/u_1[i];
       M=U/A;

       if(M<-1)                           //Finding M_plus
       {
           M_plus=0;
       }

       else if(M<1)
       {
           M_plus=(M+1)*(M+1)*0.25;
       }

       else
       {
           M_plus=M;
       }

       if(M<-1)                           //Finding P_plus
       {
           P_plus=0;
       }

       else if(M<1)
       {
           P_plus=(M+1)*0.5*P;
       }

       else
       {
           P_plus=P;
       }


        F_1_plus[i]=(R*A*M_plus);                                         //Calculation of F_plus
        F_2_plus[i]=(R*U*A*M_plus)+P_plus;
        F_3_plus[i]=(R*H*A*M_plus);




       R=u_1[i+1];
       P=(u_3[i+1]-((u_2[i+1]*u_2[i+1])/(2*u_1[i+1])))*(g-1);
       A=(g*P)/R;
       A=sqrt(A);
       U=u_2[i+1]/u_1[i+1];
       H=(u_3[i+1]+P)/u_1[i+1];
       M=U/A;

       if(M<-1)                           //Finding M_minus
       {
           M_minus=M;
       }

       else if(M<1)
       {
           M_minus=(M-1)*(M-1)*(-0.25);
       }

       else
       {
           M_minus=0;
       }



       if(M<-1)                           //Finding P_minus
       {
           P_minus=P;
       }

       else if(M<1)
       {
           P_minus=(1-M)*0.5*P;
       }

       else
       {
           P_minus=0;
       }



        F_1_minus[i+1]=(R*A*M_minus);                                         //Calculation of F_minus
        F_2_minus[i+1]=(R*U*A*M_minus)+P_minus;
        F_3_minus[i+1]=(R*H*A*M_minus);



       if(i==0)
     { temp_1=F_1_minus[i+1]+F_1_plus[i];
       temp_2=F_2_minus[i+1]+F_2_plus[i];
       temp_3=F_3_minus[i+1]+F_3_plus[i];
     }

      if(i==0)
      {
          continue;
      }

   fb_1=temp_1;
   fb_2=temp_2;
   fb_3=temp_3;

   ff_1=F_1_minus[i+1]+F_1_plus[i];
   ff_2=F_2_minus[i+1]+F_2_plus[i];
   ff_3=F_3_minus[i+1]+F_3_plus[i];



      dt=t-T;




     if(i!=0)
     {

      u_1[i]=u_1[i]-(((dt/dx))*(ff_1-fb_1));
      u_2[i]=u_2[i]-(((dt/dx))*(ff_2-fb_2));
      u_3[i]=u_3[i]-(((dt/dx))*(ff_3-fb_3));
     }
   temp_1=ff_1;
   temp_2=ff_2;
   temp_3=ff_3;


     }
 u_1[0]=u_1[1];
 u_2[0]=-u_2[1];
 u_3[0]=u_3[1];

 u_1[N+1]=u_1[N];
 u_2[N+1]=-u_2[N];
 u_3[N+1]=u_3[N];
 t=t+dt;
 }

FILE *fp,*fq,*fr,*fs;
fp=fopen("rho_AUSM_results.dat","w");
fq=fopen("u_AUSM_results.dat","w");
fr=fopen("p_AUSM_results.dat","w");
fs=fopen("ie_AUSM_results.dat","w");


for(i=1;i<N+1;i++)
{ x=((2*i)-1)*(dx*0.5);
  R=u_1[i];
  P=(u_3[i]-((u_2[i]*u_2[i])/(2*u_1[i])))*(g-1);
  e=P/((g-1)*R);
  U=u_2[i]/u_1[i];
    fprintf(fp,"%f\t%f\n",x,R);
    fprintf(fq,"%f\t%f\n",x,U);
    fprintf(fr,"%f\t%f\n",x,P);
    fprintf(fs,"%f\t%f\n",x,e);

}
fclose(fp);
fclose(fq);
fclose(fr);
fclose(fs);
}
