#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include </home/mi/palic/OMs.h>


//This code takes all uniform oriented matroids that are representatives of reorientation and permutation classes, as found by Finschi, and constructs their lower cones.




int R;			//rank
int N;			//the number of elements

int B;			//the number of bases
int nr_ints;		//the number of integers needed to store the plus (resp. minus) of a chirotope

char **bases;		//the list of bases





int makechirotopes(struct OM M, FILE *out)    //makes the lower cone of a uniform OM M, works only for OMs with at most 64 bases -- it would be too slow otherwise, anyway
{
	int c=0;
	long long int i,j;	
	long long int limit1, limit2; 

	if (B<=32)
	{
		limit1=1l<<B;
		limit2=1;
	}
	else if (B<=64)
	{   
		limit1=1l<<32;
		limit2=1l<<(B&31);
	} 
	else return -1;

	struct OM X;
	makeOM(&X);
	
	for (i=0;i<limit1;i++)			//checks for every subset of the bases of M whether it gives an OM
	{
		X.plus[0]=M.plus[0]&i;
		X.minus[0]=M.minus[0]&i;

		for (j=0;j<limit2;j++)
		{
			if (B>32)
			{
				X.plus[1]=M.plus[1]&j;
				X.minus[1]=M.minus[1]&j;
			}
			if (ischirotope(X))
			{
				writeOM(X,out);
				c++;
			}
		}
	}
		
	removeOM(&X);
	
	return c;
}


int main(int argc, char *argv[])
{
	int step;		//gets the number of a lower cone
	if (argc==2)
	{
		char *ptr;
		step = strtol(argv[1],&ptr,10);
	}
	else 
	{
		printf("An argument expected.\n");
		exit(EXIT_FAILURE);
	}
	R=3;				//The code works only for B<=64!
	N=7;
	makebases();




	FILE *in,*out;
	char text[300];
	struct OM M;
	makeOM(&M);
	int i=0;
	int c;

	sprintf(text,"lower_cones_rank%d_%delements_%d.txt",R,N,step);			//the output - all elements in the lower cone of the step-th uniform OM

	out = fopen(text,"w"); 					
	if(out==NULL)
	{
		fprintf(stderr,"error fopen():  Could not open text file %s.\n",text);
		exit(EXIT_FAILURE);
	}
	fprintf(out,"Elements of lower cone of the %d-th uniform representative (under reorientations and permutations) of rank %d on %d elements:\n",step,R,N);


	if (R>=3)
	{
		sprintf(text,"uniform_representatives_rank%d_%delements.txt",R,N);

		in = fopen(text,"r"); 		 //reads all uniform OMs that are representatives of their classes
		if(in==NULL)
		{
			fprintf(stderr,"error fopen():  Could not open text file %s.txt.\n",text);
			exit(EXIT_FAILURE);
		}

	
	
		fgets(text,300,in);

		for (i=0;i<step;i++)
			readOM(&M,in);

		if (readOM(&M,in)!=0)			//we work only with the step-th OM
		{
			c=makechirotopes(M,out);
			printf("%d chirotopes\n",c);
		}
		else printf("Mistake - the input argument is too large.\n");
		

		fclose(in);
	}
	else if (R==2)					//For R==2, there is exactly one class of uniform oriented matroids
	{
		for (i=0;i<nr_ints-1;i++)
		{
			M.plus[i]=(1l<<32)-1;
			M.minus[i]=0;
		}
		M.plus[nr_ints-1]=(1l<<(B&31))-1;
		M.minus[i]=0;

		showchirotope(M);

		makechirotopes(M,out);
	}



	fclose(out);
	removeOM(&M);
	removebases();

	
	return 0;
}
