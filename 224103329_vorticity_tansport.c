#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()

{int i,j,m,n;
float beta,dx,dy,Re;

printf("enter the value of m");
scanf("%d",&m);

printf("enter the value of n");
scanf("%d",&n);

printf("enter the value of Re");
scanf("%f",&Re);

dx=pow((m-1),-1);
dy=pow((n-1),-1);
      printf("%f\n",dx);


beta=dx/dy;

float psi[m][n],omega[m][n],psi_old[m][n],omega_old[m][n];

for(i=0;i<m;i++)       //Boundary conditions for psi
        psi[i][0]=0;

for(i=0;i<m;i++)       //Boundary conditions for psi
        psi[i][n-1]=0;

for(j=0;j<n;j++)       //Boundary conditions for psi
        psi[0][j]=0;

for(j=0;j<n;j++)       //Boundary conditions for psi
        psi[m-1][j]=0;


for(i=1;i<m-1;i++)    //Initial guess for psi
{for(j=1;j<n-1;j++)
   psi[i][j]=0;
}





float u[m][n],v[m][n];

for(j=0;j<n;j++)  //left wall
{u[0][j]=0;
v[0][j]=0;
}

for(j=0;j<n;j++)  //right wall
{u[m-1][j]=0;
v[m-1][j]=0;
}

for(i=0;i<m;i++)  //bottom wall
{u[i][0]=0;
v[i][0]=0;
}

for(i=0;i<m;i++)  //top wall
{u[i][n-1]=0;
v[i][n-1]=0;
}


float err1,err2;
int iter=1;





do
{err1=0;
err2=0;

for(j=0;j<n;j++)   //Boundary conditions for omega
    {omega[0][j]=2*(psi[1][j]-psi[0][j]);
      omega[0][j]=-omega[0][j]/pow(dx,2);
    }

for(j=0;j<n;j++)             //Boundary conditions for omega
    {omega[m-1][j]=2*(psi[m-2][j]-psi[m-1][j]);
      omega[m-1][j]=-omega[m-1][j]/pow(dx,2);
    }

for(i=0;i<m;i++)                    //Boundary conditions for omega
    {omega[i][0]=2*(psi[i][1]-psi[i][0]);
      omega[i][0]=-omega[i][0]/pow(dy,2);
    }
for(i=0;i<m;i++)                           //Boundary conditions for omega
    {omega[i][n-1]=2*(psi[i][n-2]-psi[i][n-1]+(1*dy));
      omega[i][n-1]=-omega[i][n-1]/pow(dy,2);
    }



for(i=0;i<m;i++)   // saving old values
{for(j=0;j<n;j++)
  {psi_old[i][j]=psi[i][j];
    omega_old[i][j]=omega[i][j];
}
}

for(i=1;i<m-1;i++)     //calculating psi
{for(j=1;j<n-1;j++)

{   psi[i][j]=((1.0/(2.0*(1+pow(beta,2))))*(psi[i+1][j]+psi[i-1][j]+(pow(beta,2)*(psi[i][j+1]+psi[i][j-1]))+(pow(dx,2)*omega[i][j])));

    err1=err1+(pow((psi[i][j]-psi_old[i][j]),2));
}
}

for(i=1;i<m-1;i++)     //calculating omega
{for(j=1;j<n-1;j++)
{omega[i][j]=((1.0/(2.0*(1+pow(beta,2))))*(((1.0-((psi[i][j+1]-psi[i][j-1])*((beta*Re)/4.0)))*omega[i+1][j])
                             +((1.0+((psi[i][j+1]-psi[i][j-1])*((beta*Re)/4.0)))*omega[i-1][j])
                             +((1.0+((psi[i+1][j]-psi[i-1][j])*(Re/(4.0*beta))))*(pow(beta,2)*omega[i][j+1]))
                             +((1.0-((psi[i+1][j]-psi[i-1][j])*(Re/(4.0*beta))))*(pow(beta,2)*omega[i][j-1]))));
 err2=err2+(pow((omega[i][j]-omega_old[i][j]),2));

}
}

err1=sqrt((err1/((m-2)*(n-2))));
err2=sqrt((err2/((m-2)*(n-2))));

 printf("iteration=%d\t",iter);
printf("error=%.20f\n",err1);
printf("error=%.20f\n",err2);

printf("\n\n");

iter++;
}while(err1>1e-6||err2>1e-6);

FILE *fp,*fq,*fr,*fs,*ft,*fu;

for(j=1;j<(n-1);j++)
{
for(i=1;i<(m-1);i++)
{
u[i][j]=(0.5/dy)*(psi[i][j+1]-psi[i][j-1]);
v[i][j]=(-0.5/dx)*(psi[i+1][j]-psi[i-1][j]);
}
}

fp=fopen("U_graph.plt","w");
fq=fopen("V_graph.plt","w");
fr=fopen("psi_graph.plt","w");
fs=fopen("omega_graph.plt","w");
fu=fopen("comparison_ugraph.plt","w");
ft=fopen("comparison_vgraph.plt","w");

fprintf(fr,"VARIABLES=\"x\",\"y\",\"PSI\"\n");
fprintf(fr,"ZONE T=\"BLOCK1\",i=100,j=100,F=POINT\n\n");
for(i=0;i<m;i++)
{  for(j=0;j<n;j++)
     {
        fprintf(fr,"%f\t%f\t%f\n",i*dx,j*dy,psi[i][j]);
      }
}

fprintf(fp,"VARIABLES=\"x\",\"y\",\"U\"\n");
fprintf(fp,"ZONE T=\"BLOCK1\",i=100,j=100,F=POINT\n\n");
for(i=0;i<m;i++)
{  for(j=0;j<n;j++)
     {
        fprintf(fp,"%f\t%f\t%f\n",i*dx,j*dy,u[i][j]);
      }
}

fprintf(fq,"VARIABLES=\"x\",\"y\",\"V\"\n");
fprintf(fq,"ZONE T=\"BLOCK1\",i=100,j=100,F=POINT\n\n");
for(i=0;i<m;i++)
{  for(j=0;j<n;j++)
     {
        fprintf(fq,"%f\t%f\t%f\n",i*dx,j*dy,v[i][j]);
      }
}



fprintf(fs,"VARIABLES=\"x\",\"y\",\"OMEGA\"\n");
fprintf(fs,"ZONE T=\"BLOCK1\",i=100,j=100,F=POINT\n\n");
for(i=0;i<m;i++)
{  for(j=0;j<n;j++)
     {
        fprintf(fs,"%f\t%f\t%f\n",i*dx,j*dy,omega[i][j]);
      }
}

j=(n/2)-1;
  for(i=0;i<m;i++)
     {
        fprintf(ft,"%f\t%f\n",u[i][j],i*dx);
      }

for(i=0;i<m;i++)
{
     j=(n/2)-1;
        fprintf(fu,"%f\t%f\n",v[i][j],i*dx);

}
fclose(fp);
fclose(fq);
fclose(fr);
fclose(fs);
fclose(ft);
fclose(fu);
}


