#include <stdio.h>
#include <math.h>
#define NPRE 30 
#define NPOST /*1300*//*800*//*1000*//*5000*/3000/*1400 for 33,65; 1700 for 129; 3000 for 257*/
#define ALPHA 0.33
#define NGMAX 15
#define Pi 3.141592653589793
#include "Header1new.h"
#include "Header2new.h"
void mgfas5new( double **u, int n, int maxcyc, double m,int b, int  *qq, double *x, double *y, double R0, double K, double H)
{
	int j,i,jpost,lsw,jsw,ipass; 
	double h,res,v1,v2,v3,v4,v5,v6,v7,v8,v9,h2i,foh2,rjac,omega=/*0.8*//*0.6*/1.0;

	h=1.0/(n-1); 
	h2i=1.0/(h*h);
	foh2=-4.0*h2i;
	
	rjac=/*0.99555*//*0.999*//*0.167080002*//*0.5*/0.9/*0.9 for hyperbolic  K=0~0.6, R=1;0.9 for 33,65; 0.99for 129; 0.9999 for 257*/;
	for (jpost=1;jpost<=NPOST;jpost++)
	{ 
	  for (j=(int)(n*H-H+1)-1;j>=b+2;j--)//inside region
		{
		  for(i=2;i<qq[j+1];i++)/*for(i=2;i<qq[j-1];i++) for waist;  for(i=2;i<qq[j+1];i++) for barrel */
			{
				                v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];
						/*	if (i==qq[j-1]&&j+2<=(int)(n*H-H+1))
						  {
						    v7=2.0*(u[i-1][j-1]-u[i][j-1]);
						   
						    
						  }

						if (i>=qq[j-1]+1&&j+2<=(int)(n*H-H+1))
						  {
						    v2=v1;
						     v1=u[i][j+2]-u[i][j+1];
						    v6=2.0*(u[i+1][j+1]-u[i-1][j+1]);
						    v7=2.0*(u[i-1][j]-u[i+1][j]);
						    }*/
						
						
				        
						res=K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);

				
				v3=2.0*K*cos(2.0*u[i][j])-0.5*(K-1.0)*pow((i-1)*(v1+v2),2)*cos(2.0*u[i][j])+(K-1.0)*(2.0*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)
					+(K-1.0)*(-2.0*(i-1)*(i-1)*sin(2.0*u[i][j]))*0.5*(v1+v2)*(v4+v5)+(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*(1.0*(i-1))*(v4+v5)+(K-1.0)*cos(2.0*u[i][j])*(i-1)*(v1+v2)
					-2.0*(-(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2))*(v4-v5)+4.0*((K*pow(cos(1.0*u[i][j]),2)+1.0*pow(sin(1.0*u[i][j]),2))*pow(1.0*(i-1),2))
					-2.0*(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2)*(v1-v2)+4.0*(K*pow(sin(1.0*u[i][j]),2)+1.0*pow(cos(1.0*u[i][j]),2))*pow(1.0*(i-1),2)+(2.0*(K-1.0)*pow(1.0*(i-1),2)*cos(2.0*u[i][j]))*0.5*(v6+v7);
						u[i][j] -= omega*res/v3;
						if(u[i][j]>=Pi) u[i][j]=Pi;
						  if(u[i][j]<=0.0) u[i][j]=0.0;

	/*if (i>=qq[j-1]) {u[i][j]=u[qq[j-1]][j-1]+(u[qq[j]][j]-u[qq[j-1]][j-1])*(i-qq[j-1])/(qq[j]-qq[j-1]);*//*printf("YES. u[%d][%d]=%f\n", i,j,u[i][j]);}*/
	
	}



		  /*	i=qq[j]-1;
			
				        v1=u[i][j+1]-u[i][j];
						v2=u[i][j+1]-u[i][j];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i+1][j+1]-u[i-1][j+1];
				        
						res=K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);

				
				v3=2.0*K*cos(2.0*u[i][j])-0.5*(K-1.0)*pow((i-1)*(v1+v2),2)*cos(2.0*u[i][j])+(K-1.0)*(2.0*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)
					+(K-1.0)*(-2.0*(i-1)*(i-1)*sin(2.0*u[i][j]))*0.5*(v1+v2)*(v4+v5)+(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*(1.0*(i-1))*(v4+v5)+(K-1.0)*cos(2.0*u[i][j])*(i-1)*(v1+v2)
					-2.0*(-(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2))*(v4-v5)+4.0*((K*pow(cos(1.0*u[i][j]),2)+1.0*pow(sin(1.0*u[i][j]),2))*pow(1.0*(i-1),2))
					-2.0*(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2)*(v1-v2)+4.0*(K*pow(sin(1.0*u[i][j]),2)+1.0*pow(cos(1.0*u[i][j]),2))*pow(1.0*(i-1),2)+(2.0*(K-1.0)*pow(1.0*(i-1),2)*cos(2.0*u[i][j]))*0.5*(v6+v7);
						u[i][j] -= omega*res/v3;
							if(u[i][j]>=Pi) u[i][j]=Pi;
							if(u[i][j]<=0.0) u[i][j]=0.0;*/
	



			
		}








	  

	  

	  for (j=1+b;j>=2;j--)//lattices near the lower boundary
		{
		  for (i=b+2;i<qq[j+1];i++)/* for (i=b+2;i<qq[j-1];i++) for waist; for (i=b+2;i<qq[j+1];i++) for barrel */
			{
				        v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];

						/*if (i==qq[j-1])
						  {
						    v7=2.0*(u[i-1][j-1]-u[i][j-1]);
						    
						  }

						if (i>=qq[j-1]+1)
						  {
						    v2=v1;
						    v1=u[i][j+2]-u[i][j+1];
						    v6=2.0*(u[i+1][j+1]-u[i-1][j+1]);
						    v7=2.0*(u[i-1][j]-u[i+1][j]);
						    }*/
				        
							res=K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
				
				v3=2.0*K*cos(2.0*u[i][j])-0.5*(K-1.0)*pow((i-1)*(v1+v2),2)*cos(2.0*u[i][j])+(K-1.0)*(2.0*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)
					+(K-1.0)*(-2.0*(i-1)*(i-1)*sin(2.0*u[i][j]))*0.5*(v1+v2)*(v4+v5)+(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*(1.0*(i-1))*(v4+v5)+(K-1.0)*cos(2.0*u[i][j])*(i-1)*(v1+v2)
					-2.0*(-(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2))*(v4-v5)+4.0*((K*pow(cos(1.0*u[i][j]),2)+1.0*pow(sin(1.0*u[i][j]),2))*pow(1.0*(i-1),2))
					-2.0*(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2)*(v1-v2)+4.0*(K*pow(sin(1.0*u[i][j]),2)+1.0*pow(cos(1.0*u[i][j]),2))*pow(1.0*(i-1),2)+(2.0*(K-1.0)*pow(1.0*(i-1),2)*cos(2.0*u[i][j]))*0.5*(v6+v7);
						
						u[i][j] -= omega*res/v3;
							if(u[i][j]>=Pi) u[i][j]=Pi;
							if(u[i][j]<=0.0) u[i][j]=0.0;
	/*	if (i>=qq[j-1]) u[i][j]=u[qq[j-1]][j-1]+(u[qq[j]][j]-u[qq[j-1]][j-1])*(i-qq[j-1])/(qq[j]-qq[j-1]);*/
			}
		}


	
		/*	for(j=(int)(n*H-H+1)-b;j<(int)(n*H-H+1);j++)
		{
			for (i=2;i<qq[(int)(n*H-H+1)-b]-b;i++)
			{
				        v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];
				       
							res=K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
				
				v3=2.0*K*cos(2.0*u[i][j])-0.5*(K-1.0)*pow((i-1)*(v1+v2),2)*cos(2.0*u[i][j])+(K-1.0)*(2.0*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)
					+(K-1.0)*(-2.0*(i-1)*(i-1)*sin(2.0*u[i][j]))*0.5*(v1+v2)*(v4+v5)+(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*(1.0*(i-1))*(v4+v5)+(K-1.0)*cos(2.0*u[i][j])*(i-1)*(v1+v2)
					-2.0*(-(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2))*(v4-v5)+4.0*((K*pow(cos(1.0*u[i][j]),2)+1.0*pow(sin(1.0*u[i][j]),2))*pow(1.0*(i-1),2))
					-2.0*(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2)*(v1-v2)+4.0*(K*pow(sin(1.0*u[i][j]),2)+1.0*pow(cos(1.0*u[i][j]),2))*pow(1.0*(i-1),2)+(2.0*(K-1.0)*pow(1.0*(i-1),2)*cos(2.0*u[i][j]))*0.5*(v6+v7);
						
						u[i][j] -= omega*res/v3;
						if(u[i][j]>=Pi) u[i][j]=Pi;
	if(u[i][j]<=0.0) u[i][j]=0.0;
			}
			}*/
			omega=(jpost==1?1.0/(1.0-0.5*rjac*rjac):1.0/(1.0-0.25*rjac*rjac*omega));
	
	}

	
}
