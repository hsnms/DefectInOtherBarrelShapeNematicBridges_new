#include <math.h>
#include <stdio.h>
#include "Header2new.h"
/*#define EPI 1.0e-7*/ 
#define Pi 3.141592653589793
void **nrfunc4new(double **f, double **u, int n, int *qq,int b,double R0, double K, double H)
{
  double h,v1,v2,v4,v5,v6,v7,v8,v9, v;
	int i,j;
	h=1.0/(n-1); 	

	for(j=1;j<=(int)(n*H-H+1);j++)
		for(i=1;i<=(int)(n*R0-R0+1);i++)
		  f[i][j]=0.0;//initially integrand is zero
	
	//inside, only half of the region
	for(j=b+2;j<=(int)((n*H-H+1))/2;j++)/*for(j=b+2;j<(int)((n*H-H+1))-(n-1)/8;j++) for waist*/
	  { v=0.0;
	    for (i=2;i<qq[j-1]-(n-1)/8;i++)/*qq[(int)((n*H-H+1))];i++)*//*(i=qq[j-1]-(n-1)/8;i<qq[j-1]-(n-1)/16;i++)*//*for (i=2;i<qq[j-1]-(n-1)/8;i++) for waist*/
			  {
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
		        
			
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
			
	        
	        	v=v+f[i][j]/((n-1)*(n-1));
			  }
	    /* printf("j=%d\tqq[%d]=%d\tu[qq[%d]][%d]=%f\tv=%f\n",j,j,qq[j],j,j,u[qq[j]][j],v);	*/		  }


	//inside, the other half of the region
	for(j=(int)((n*H-H+1))/2+1;j<(int)((n*H-H+1))-(n-1)/8;j++)/*for(j=b+2;j<(int)((n*H-H+1))-(n-1)/8;j++) for waist*/
	  { v=0.0;
	    for (i=2;/*i<qq[j-1]-(n-1)*R0/4*/i<qq[(int)((n*H-H+1))];i++)/*(i=qq[j-1]-(n-1)/8;i<qq[j-1]-(n-1)/16;i++)*//*for (i=2;i<qq[j-1]-(n-1)/8;i++) for waist*/
			  {
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
		        
			
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
			
	        
	        	v=v+f[i][j]/((n-1)*(n-1));
			  }
	    /* printf("j=%d\tqq[%d]=%d\tu[qq[%d]][%d]=%f\tv=%f\n",j,j,qq[j],j,j,u[qq[j]][j],v);	*/		  }




	    
	//lattices near the lower boundary
	
		for(j=2;j<=1+b;j++)
		{
		  for (i=b+2;i<qq[j]-(n-1)/8;i++) /* for (i=b+2;i<qq[j]-(n-1)/8;i++) for waist*/
			{
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			        
			
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
				
		}
		}	

        

		
		



		for (i=b+2;i<=qq[1]-1-(n-1)/8;i++)/*for (i=b+2;i<=qq[1]-1-(n-1)/8;i++) for waist*/
	{
			v1=u[i][2]-u[i][1];
		    v2=u[i][2]-u[i][1];
			v4=u[i+1][1]-u[i][1];
			v5=u[i][1]-u[i-1][1];
				
			j=1;
		
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
			

				}

		//all the other lines and points
	for (j=b+2;j<=(int)(n*H-H+1)-1;j++)
	{
		    v1=u[1][j+1]-u[1][j];
		    v2=u[1][j]-u[1][j-1];
			v4=u[2][j]-u[1][j];
			v5=u[2][j]-u[1][j];
		        	       
			i=1;
		
			f[i][j]=((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
			
				}




		


	for (j=2;j<=b;j++)
	{
		    v1=u[1+b][j+1]-u[1+b][j];
		    v2=u[1+b][j]-u[1+b][j-1];
			v4=u[b+2][j]-u[1+b][j];
			v5=u[b+2][j]-u[1+b][j];
	        	        
			i=1+b;
		
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
	}


		for (i=2;i<=b;i++)
	{
		    v1=u[i][b+2]-u[i][b+1];
		    v2=u[i][b+2]-u[i][b+1];
			v4=u[i+1][b+1]-u[i][b+1];
			v5=u[i][b+1]-u[i-1][b+1];
	        	        
			j=b+1;
		
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
			        }




         
	  
			v1=u[1][(int)(n*H-H+1)]-u[1][(int)(n*H-H+1)-1];
		    v2=u[1][(int)(n*H-H+1)]-u[1][(int)(n*H-H+1)-1];
			v4=u[2][(int)(n*H-H+1)]-u[1][(int)(n*H-H+1)];
			v5=u[2][(int)(n*H-H+1)]-u[1][(int)(n*H-H+1)];
				       
			i=1;j=(int)(n*H-H+1);
		
			f[i][j]=((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
	        
	
       
	
            
	      
	
           
				v1=u[1+b][2]-u[1+b][1];
		    v2=u[1+b][2]-u[1+b][1];
			v4=u[b+2][1]-u[1+b][1];
			v5=u[b+2][1]-u[1+b][1];
			
			i=1+b;j=1;
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
			
	        
			
	        
			
				i=1;j=b+1;
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i+1][j]-u[i][j];
			f[i][j]=((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
			
			i=1+b;j=b+1;
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow((v1+v2),2)/(4.0*h)-((K-1.0)*sin(2.0*u[i][j])*(i-1))*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j]))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*(v1+v2)/(2.0*h);
	        		
	return 0;
}
