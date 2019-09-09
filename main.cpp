#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <fstream>
using namespace std;
#define  EPSILON    0.0000001
void norm_error(double *error, int n ,double* one_norm,double* two_norm,double* inf_norm)
{
	double sum=0.0;
	for(int i=0;i<n;i++)
	{
		sum = sum + error[i];
	}
	*one_norm=sum;
	sum=0.0;
	for(int i=0;i<n;i++)
	{
		sum = sum + (error[i]*error[i]);
	}
	*two_norm = sqrt(sum);
	double temp=0.0;
	for(int i=0;i<n;i++)
	{
		for(int j=i;j<n;j++)
		{
			if(error[i]>error[j])
			{
				temp=error[i];
				error[i]=error[j];
				error[j]=temp;
			}
		}
	}
	*inf_norm=error[n-1];
}
void initial(double* T,double* s,double** a,double* b,double* oldT,int piece,double h,double w,const double left_bc,const double right_bc,const double up_bc,const double down_bc,const double bc_num,const double source_hot)
{
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			T[i*piece+j]=0.0;
			s[i*piece+j]=0.0;
			oldT[i*piece+j]=0.0;
		}
	}
	for(int i=0;i<piece*piece;i++)
	{
		for(int j=0;j<piece*piece;j++)
		{
			a[i][j]=0.0;
		}
		b[i]=0.0;
	}

	s[(piece/2)*piece+(piece/2)]=source_hot;
	for(int i=0;i<piece;i++)
	{
		T[i]=up_bc;
		T[(piece-1)*piece+i]=down_bc;
	}
	for(int i=0;i<piece;i++)
	{
		T[i*piece]=left_bc;
		T[i*piece+(piece-1)]=right_bc;
	}
	T[0]=25.0,T[piece-1]=25.0,T[(piece-1)*piece]=25.0,T[piece*piece-1]=25.0;
	for(int i=1;i<piece-1;i++)
	{
		for(int j=1;j<piece-1;j++)
		{
			T[i*piece+j]=0.0;
		}
	}
	/*a[0][0]=1.0,a[piece-1][piece-1]=1.0,a[(piece-1)*piece][(piece-1)*piece]=1.0,a[piece*piece-1][piece*piece-1]=1.0;
	for(int i=1;i<piece-1;i++)
	{
		a[i][i]=1.0-w;
		a[i][piece+i]=w;
	}
	for(int i=1;i<piece-1;i++)
	{
		a[i*piece][i*piece]=1.0*w/4;
		a[i*piece+(piece-1)][i*piece+(piece-1)]=1.0*w/4;
		a[(piece-1)*piece+i][(piece-1)*piece+i]=1.0*w/4;
	}*/
	for(int i=1;i<piece-1;i++)
	{
		for(int j=1;j<piece-1;j++)
		{
			a[i*piece+j][(i-1)*piece+j]=1.0*w/4;
			a[i*piece+j][(i+1)*piece+j]=1.0*w/4;
			a[i*piece+j][i*piece+(j-1)]=1.0*w/4;
			a[i*piece+j][i*piece+(j+1)]=1.0*w/4;
			a[i*piece+j][i*piece+j]=-4.0*w/4;
		}
	}
	for(int i=1;i<piece-1;i++)
	{
		b[i]=w*h*bc_num;
	}
	for(int i=1;i<piece-1;i++)
	{
		for(int j=1;j<piece-1;j++)
		{
			b[i*piece+j]=s[i*piece+j]*h*h*(-w/4);
		}
	}

	double one_error=0.0,two_error=0.0,inf_error=0.0;
	double* r=new double[piece*piece];
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			r[i*piece+j]=T[i*piece+j]-oldT[i*piece+j];
		}
	}
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			oldT[i*piece+j] = T[i*piece+j];
		}
	}
	/*printf("T=\n");
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			printf("%lf ",T[i*piece+j]);
		}
		printf("\n");
	}
	printf("a=\n");
	for(int i=0;i<piece*piece;i++)
	{
		for(int j=0;j<piece*piece;j++)
		{
			printf("%lf ",a[i][j]);
		}
		printf("\n");
	}
	printf("b=\n");
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			printf("%lf ",b[i*piece+j]);
		}
		printf("\n");
	}*/
}
void SOR(double* T,double* s,double* oldT,int piece,double h,double w,const double bc_num,double* error)
{
	for(int i=1;i<piece-1;i++)
	{
		T[i]=T[i]+w*(T[piece+i]-T[i]+h*(bc_num));
	}
	for(int i=1;i<piece-1;i++)
	{
		for(int j=1;j<piece-1;j++)
		{
			T[i*piece+j]=(T[(i-1)*piece+j]+T[i*piece+(j-1)]+T[i*piece+(j+1)]+T[(i+1)*piece+j]+(h*h*s[i*piece+j]))*w/4+(T[i*piece+j]-w*T[i*piece+j]);
		}
	}
	double one_error=0.0,two_error=0.0,inf_error=0.0;
	double* r=new double[piece*piece];
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			r[i*piece+j]=fabs(T[i*piece+j]-oldT[i*piece+j]);
		}
	}
	norm_error(r,piece*piece,&one_error,&two_error,&inf_error);
	*error=one_error;
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			oldT[i*piece+j] = T[i*piece+j];
		}
	}
	delete[] r;
}

void conjugate(double **A, double *x,double* oldT, double *b ,int piece, double* error,int* iteration)
{
	double* h=new double[piece*piece];
	double* d=new double[piece*piece];
	double* g=new double[piece*piece];
	double* t=new double[piece*piece];
	double* r=new double[piece*piece];
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			h[i*piece+j]=0.0;
			d[i*piece+j]=0.0;
			g[i*piece+j]=0.0;
			t[i*piece+j]=0.0;
			r[i*piece+j]=0.0;
		}
	}
	for(int i0=0;i0<piece;i0++)// t = Ax.
	{
		for(int i1=0;i1<piece;i1++)
		{
			t[i0*piece+i1] = 0.0;
			for(int j0=0;j0<piece;j0++)
			{
				for(int j1=0;j1<piece;j1++)
				{
					t[i0*piece+i1]=t[i0*piece+i1]+A[i0*piece+i1][j0*piece+j1]*x[j0*piece+j1];
				}
			}
		}
	}
	/*printf("x[]=\n");
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			printf("%lf ",x[i*piece+j]);
		}
		printf("\n");
	}
	system("pause");*/
	/*printf("Ax[]=\n");
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			printf("%lf ",t[i*piece+j]);
		}
		printf("\n");
	}
	system("pause");*/
	for(int i=0;i<piece;i++)// d = b - Ax.//g = -d.
	{
		for(int j=0;j<piece;j++)
		{
			d[i*piece+j]=b[i*piece+j]-t[i*piece+j];
			g[i*piece+j]=d[i*piece+j]*(-1.0);
		}
	}
	/*printf("b[]=\n");
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			printf("%lf ",b[i*piece+j]);
		}
		printf("\n");
	}
	system("pause");
	printf("d[]=\n");
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			printf("%lf ",d[i*piece+j]);
		}
		printf("\n");
	}
	system("pause");
	printf("g[]=\n");
	for(int i=0;i<piece;i++)
	{
		for(int j=0;j<piece;j++)
		{
			printf("%lf ",g[i*piece+j]);
		}
		printf("\n");
	}
	system("pause");*/
	double one_error=0.0,two_error=0.0,inf_error=0.0;
	double newG2=0.0, oldG2=0.0, dh=0.0, alpha=0.0, beta=0.0;
	for((*iteration)=0;(*error)>EPSILON;(*iteration)++)
	{
		newG2=0.0, oldG2=0.0, dh=0.0, alpha=0.0, beta=0.0;
		for(int i0=0;i0<piece;i0++)// h = Ad.
		{
			for(int i1=0;i1<piece;i1++)
			{
				h[i0*piece+i1] = 0.0;
				for(int j0=0;j0<piece;j0++)
				{
					for(int j1=0;j1<piece;j1++)
					{
						h[i0*piece+i1]=h[i0*piece+i1]+A[i0*piece+i1][j0*piece+j1]*d[j0*piece+j1];
					}
				}
			}
		}
		/*printf("h[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",h[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/
		for(int i=0;i<piece*piece;i++)// <d,h>.
		{
			dh=dh+d[i]*h[i];
		}
		/*if(dh==0.0)
		{
			*error=0.0;
			printf("error:%lf\n",*error);
			break;
		}*/
		for(int i=0;i<piece*piece;i++)//oldG2 = <g, g>
		{
			oldG2=oldG2+g[i]*g[i];
		}
		alpha = oldG2/dh;
		/*if(alpha<-1.0||alpha>1.0)
		{
			if(alpha<0.0)
			{
				alpha=-1.0;
			}
			else
			{
				alpha=1.0;
			}
		}*/

		//printf("dh=%lf ,oldG2=%lf ,alpha=%lf\n",dh,oldG2,alpha);

		for(int i=0;i<piece;i++)//t = alpha*d
		{
			for(int j=0;j<piece;j++)
			{
				t[i*piece+j]=d[i*piece+j]*alpha;
			}
		}
		/*printf("alpha*d[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",t[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/
		for(int i=0;i<piece;i++) //x = x + alpha*d
		{
			for(int j=0;j<piece;j++)
			{
				x[i*piece+j]=x[i*piece+j]+t[i*piece+j];
			}
		}
		/*printf("x[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",x[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/

		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				r[i*piece+j]=fabs(x[i*piece+j]-oldT[i*piece+j]);
			}
		}
		/*printf("r[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",r[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/
		norm_error(r,piece*piece,&one_error,&two_error,&inf_error);
		*error=one_error;
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				oldT[i*piece+j] = x[i*piece+j];
			}
		}
		for(int i=0;i<piece;i++)//t = alpha*h = alpha*A*d
		{
			for(int j=0;j<piece;j++)
			{
				t[i*piece+j]=h[i*piece+j]*alpha;
			}
		}
		/*printf("alpha*h[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",t[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/
		for(int i=0;i<piece;i++) //g =g + alpha*A*d
		{
			for(int j=0;j<piece;j++)
			{
				g[i*piece+j]=g[i*piece+j]+t[i*piece+j];
			}
		}
		/*printf("g[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",g[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/
		for(int i=0;i<piece*piece;i++)//newG2 = <g, g>
		{
			newG2=newG2+g[i]*g[i];
		}
		beta = newG2/oldG2;
		//printf("newG2:%lf ,oldG2:%lf ,beta:%lf\n",newG2,oldG2,beta);
		for(int i=0;i<piece;i++)//t = beta*d;
		{
			for(int j=0;j<piece;j++)
			{
				t[i*piece+j]=d[i*piece+j]*beta;
			}
		}
		/*printf("beta*d[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",t[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/
		for(int i=0;i<piece;i++) //d = -g + beta*d
		{
			for(int j=0;j<piece;j++)
			{
				d[i*piece+j]=g[i*piece+j]*(-1)+d[i*piece+j]*beta;
			}
		}
		/*printf("d[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",d[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/

		
		/*printf("iteration=%d\n",iteration);
		printf("error:%lf\n",*error);*/
		/*printf("x[]=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",x[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/
		/*printf("alpha=%lf ,beta=%lf\n",alpha,beta);
		printf("d=\n");
		for(int i=0;i<n;i++)
		{
			printf("%lf ",d[i]);
		}
		printf("\n");
		printf("g=\n");
		for(int i=0;i<n;i++)
		{
			printf("%lf ",g[i]);
		}
		printf("\n");*/
	}
	delete[] h;
	delete[] d;
	delete[] g;
	delete[] t;
	delete[] r;
}
void initial_test(double** a,double* x,double* b,double* oldx,int piece,double h,double w)
{
	for(int i0=0;i0<piece;i0++)
	{
		for(int i1=0;i1<piece;i1++)
		{
			for(int j0=0;j0<piece;j0++)
			{
				for(int j1=0;j1<piece;j1++)
				{
					a[i0*piece+i1][j0*piece+j1]=1.0;
				}
			}
			x[i0*piece+i1]=1.0;
			b[i0*piece+i1]=440.0;
			oldx[i0*piece+i1]=0.0;
		}
	}
}

int main()
{
	int piece=21;double h=(1/(double)piece);
	double* T=new double[piece*piece];
	double* s=new double[piece*piece];
	double** a=new double*[piece*piece];
	double* b=new double[piece*piece];
	double* oldT=new double[piece*piece];
	for(int i=0;i<piece*piece;i++)
	{
		a[i]=new double[piece*piece];
		oldT[i]=0.0;
	}

	int iteration=0,MAX_W=10;
	int* iter=new int[2*MAX_W];
	double error=INT_MAX;
	const double left_bc=30.0,right_bc=30.0,up_bc=0.0,down_bc=20.0,bc_num=0.0,source_hot=1.0;
	initial(T,s,a,b,oldT,piece,h,1.0,left_bc,right_bc,up_bc,down_bc,bc_num,source_hot);
	printf("T=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",T[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");
	ofstream out_file("HW5¤ÀªRµ²ªG.txt");
	printf("SOR:\n");
	for(int i=0;i<MAX_W;i++)
	{
		double w=1.0+0.1*i;
		printf("w=%lf\n",w);
		iteration=0;
		initial(T,s,a,b,oldT,piece,h,w,left_bc,right_bc,up_bc,down_bc,bc_num,source_hot);
		while(error>EPSILON)
		{
			SOR(T,s,oldT,piece,h,w,bc_num,&error);
			iteration++;
		}
		iter[i]=iteration;
		printf("iteration=%d\n",iteration);
		error=INT_MAX;
		/*printf("T=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",T[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");*/
	}
	printf("T=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",T[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");
	printf("---------------------\n");
	printf("conjugate:\n");
	for(int i=0;i<MAX_W;i++)
	{
		double w=1.0+0.1*i;
		printf("w=%lf\n",w);
		iteration=0;
		initial(T,s,a,b,oldT,piece,h,w,left_bc,right_bc,up_bc,down_bc,bc_num,source_hot);
		conjugate(a,T,oldT,b,piece,&error,&iteration);
		iter[MAX_W+i]=iteration;
		printf("iteration=%d\n",iteration);
		error=INT_MAX;
	}
	printf("T=\n");
		for(int i=0;i<piece;i++)
		{
			for(int j=0;j<piece;j++)
			{
				printf("%lf ",T[i*piece+j]);
			}
			printf("\n");
		}
		system("pause");
	/*piece=21;
	double** test_a=new double*[piece*piece];
	for(int i=0;i<piece*piece;i++)
	{
		test_a[i]=new double[piece*piece];
	}
	double* test_b=new double[piece*piece];
	double* test_x=new double[piece*piece];
	double* test_oldx=new double[piece*piece];
	initial_test(test_a,test_x,test_b,test_oldx,piece,h,1.0);
	conjugate(test_a,test_x,test_oldx,test_b,piece,&error);*/
	out_file<<"w=\t\t\t";
	for(int i=0;i<MAX_W;i++)
	{
		out_file<<1.0+0.1*i<<"\t\t";
	}
	out_file<<endl;
	for(int i=0;i<=1;i++)
	{
		if(i==0)
		{
			out_file<<"SQR:      |\t";
		}
		else
		{
			out_file<<"conjugate:|\t";
		}
		for(int j=0;j<MAX_W;j++)
		{
			if(i*MAX_W+j==0||i*MAX_W+j==1)
			{
				out_file<<iter[i*MAX_W+j]<<"\t";
			}
			else
			{
				out_file<<iter[i*MAX_W+j]<<"\t\t";
			}
		}
		out_file<<endl;
	}
	initial(T,s,a,b,oldT,piece,h,1.8,left_bc,right_bc,up_bc,down_bc,bc_num,source_hot);
	iteration=0;
	out_file<<"T=\n";
	while(error>EPSILON)
	{
		SOR(T,s,oldT,piece,h,1.8,bc_num,&error);
		iteration++;
		if((iteration%2==0)&&iteration<30)
		{
			for(int i=0;i<piece;i++)
			{
				for(int j=0;j<piece;j++)
				{
					out_file<<T[i*piece+j]<<" ";
				}
				out_file<<"\n";
			}
		}
	};
	out_file.close();
	system("pause");
	return 0;
}