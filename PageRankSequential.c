/*
  PIPERIDIS ANESTIS
  AEM : 8689
  PARALLEL AND DISTRIBUTED
  PAGERANK SEQUENTIAL

*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define ITERATIONS 110
#define ALPHA 0.85

struct timeval startwtime, endwtime;
double seq_time;

double **A,*b,*pageranksCurrent,*pageranksPrevious;
int *nodes_vector;
int nodes;
char *nodesFileName,*adjFileName;

void allMemoryAllocation();
void initializeAll();
void createOutput();
void check_arguments(int arguments,char **args);
void create_A_Matrix();
void printA();
void printNodesVector();
void printPageranks();
void gaussSeidelPagerank();
void createLinkMatrix();
void create_Left_Part_Of_Equation();
void create_Right_Part_Of_Equation();
void print_bVector();
double convergence();
void sortingResult();
void freeAll();

int main(int argc,char **argv)
{
  // 2 FILES .txt MUST BE AS INPUT

  check_arguments(argc,argv);
  
  // FIRST INPUT THE .txt FILE WITH THE NODES
  // SECOND INPUT THE .txt FILE WITH THE ADJ MATRIX

  nodesFileName = argv[1];
  adjFileName = argv[2];

  // SEE HOW MANY NODES THE DATASET HAS

  FILE *node_text = fopen(nodesFileName,"r");
  if ( node_text == NULL )
  {
    printf("Cannot open file : %s\n",nodesFileName);
    exit(1);
  }
  int return_value = fscanf(node_text,"%d",&nodes);
  printf("This dataset contains %d nodes\n",nodes);
  fclose(node_text); 

  allMemoryAllocation();

  initializeAll();

  create_A_Matrix();

  // ------------------- START COUNTING TIME -------------------------

  gettimeofday (&startwtime, NULL);

  gaussSeidelPagerank();

  sortingResult();

  gettimeofday (&endwtime, NULL);

  // ------------------- STOP COUNTING TIME --------------------------

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("PageRank Sequential Time = %0.8f\n", seq_time);

  //printA();

  //printNodesVector();

  //printPageranks();

  //print_bVector();

  //createOutput();
  
  freeAll();
  
}

/*

  void gaussSeidelPagerank()

  solve the linear system Ax=b with gauss seidel
  
  we follow these steps:

  calculate the Link Matrix
  calculate the I-a*P

  for each iteration:

  calculate the b vector of the Ax=b
  X(k+1) = 1/Aii*(bi - SUM Aij*Xj(k+1) - SUM Aij*Xj(k))
  calculate the error from the vectors X(k+1), X(k)
 
*/

void gaussSeidelPagerank()
{
  double sum;
  createLinkMatrix();
  create_Left_Part_Of_Equation();

  for(int loops = 0 ;loops<ITERATIONS;loops++)
  {
    create_Right_Part_Of_Equation();
    for(int i=0;i<nodes;i++)
    {
      pageranksPrevious[i] = pageranksCurrent[i];
      sum = 0.0;
      for(int j=0;j<=nodes;j++)
        if(i!=j)
          sum +=  A[i][j]*pageranksCurrent[j]; 
      pageranksCurrent[i] = (b[i] - sum)/(double)(A[i][i]);
    }   
    printf("Iteration %d  , convergence : %0.18f \n",loops+1,convergence());
  }  
}

/*
  void create_Left_Part_Of_Equation()

  we create the matrix  (I-alpha*P)
  where I is the identity matrix and the alpha is 0.85 
  and P is the Link Matrix

  the matrix (I - a*P)
  is the the matrix A of the linear system Ax=b

*/

void create_Left_Part_Of_Equation()
{
  for(int i=0;i<nodes;i++)
    for(int j=0;j<nodes;j++)
      if(i==j)
        A[i][j] = 1 + (-1)*ALPHA*A[i][j];
      else
        A[i][j] =  (-1)*ALPHA*A[i][j];
}

void create_Right_Part_Of_Equation()
{
  // CREATING THE b OF THE EQUATION Ax = b
  // b = (1-alpha)*x

  for(int i=0;i<nodes;i++)  
    b[i] = (1-ALPHA)*pageranksPrevious[i];
}

/*
  double convergence() 

  calculates the :

  ||X(k+1) - X(k)||2

*/

double convergence()
{
  double sum = 0.0;
  for (int i=0;i<nodes;i++) 
    sum += (pageranksCurrent[i] - pageranksPrevious[i])*(pageranksCurrent[i] - pageranksPrevious[i]);
  
  return sqrt(sum);
}

/*
  void createLinkMatrix()

  for each column we count how many 1 we have .. and 
  then we divive each value with the sum.

*/

void createLinkMatrix()
{
  double sum;
  for(int i=0;i<nodes;i++)
  {
    sum = 0.0;
    for(int j=0;j<nodes;j++)
      sum += A[j][i];
    if(sum != 0.0 )
    {
      for(int j=0;j<nodes;j++)
        A[j][i] = A[j][i]/(double)sum;  
    }
  }
}

/*
  void sortingResult()

  sorts the pagerank vector and also the nodes vector
  in order to create the output file

*/
void sortingResult() 
{
  double temp1;
  int temp2;
  for (int i = 0; i < nodes-1; i++)        
    for (int j = 0; j < nodes-i-1; j++) 
      if (pageranksCurrent[j] <= pageranksCurrent[j+1])
      {
        temp1 = pageranksCurrent[j];
        pageranksCurrent[j] = pageranksCurrent[j+1];
        pageranksCurrent[j+1] = temp1;

        temp2 = nodes_vector[j];
        nodes_vector[j] = nodes_vector[j+1];
        nodes_vector[j+1] = temp2;       
      }
}

/*  
  void allMemoryAllocation()

  allocates all the memory we need to do our job.

  nodes_vector : nodes x 1 vector
  pageranksCurrent : nodes x 1 vector
  pageranksPrevious : nodes x 1 vector
  b : nodes x 1 vector
  A : nodes x nodes matrix
  
*/
void allMemoryAllocation()
{
  nodes_vector = (int *)malloc(nodes*sizeof(int));
  if (nodes_vector == NULL )
  {
    printf("Error in memory allocation\n");
    exit(1);
  }

  pageranksCurrent = (double *)malloc(nodes*sizeof(double));
  if (pageranksCurrent == NULL )
  {
    printf("Error in memory allocation\n");
    exit(1);
  }

  pageranksPrevious = (double *)malloc(nodes*sizeof(double));
  if (pageranksPrevious == NULL )
  {
    printf("Error in memory allocation\n");
    exit(1);
  }

  b = (double *)malloc(nodes*sizeof(double));
  if (b == NULL )
  {
    printf("Error in memory allocation\n");
    exit(1);
  }

  A = (double **)malloc(nodes*sizeof(double *));
  if ( A == NULL )
  {
    printf("Error in memory allocation\n");
    exit(1);
  }
  for (int i=0;i<nodes;i++)
  {
    A[i] = (double *)malloc(nodes*sizeof(double));
    if (A[i] == NULL )
    {
      printf("Error in memory allocation\n");
      exit(1);
    }
  }

}

/*
  void create_A_matrix()

  turn an adj list to matrix
  the matrix contains 0 if a node j doesnt have a link to node i
  and has 1 if the node j has a link to node i

*/

void create_A_Matrix()
{
  for(int i=0;i<nodes;i++)
    for(int j=0;j<nodes;j++)
      A[i][j] = 0.0;

  //printf("%s\n",adjFileName);

  int FILE_VALUE;
  int return_value;
  char UNUSED_STRING[20];
  FILE *file = fopen(adjFileName,"r");
  if ( file == NULL )
  {
    printf("Error in opening the file : %s\n",adjFileName);
    exit(1);
  }
  return_value = fscanf(file,"%s",UNUSED_STRING);

  for(int i=0;i<nodes;i++)
  {
    for(;;)
    {
      return_value = fscanf(file,"%d",&FILE_VALUE);
      if( FILE_VALUE != -1 )
      {
        A[FILE_VALUE][i] = 1.0;
      }
      else
      {
        return_value = fscanf(file,"%s",UNUSED_STRING);
        break;
      }
    }
  }
  fclose(file);
}

/*
  void initializeAll()

  initialize :

  nodes_vector[i] = 0,1,....nodes-1
  pageranks = 1/nodes ..... Starting values 
 
*/

void initializeAll()
{
  for(int i=0;i<nodes;i++)
  {
    nodes_vector[i] = i;
    pageranksCurrent[i] = 1/(double)nodes ;
    pageranksPrevious[i] = pageranksCurrent[i];
  }
}

/*
  void freeAll()

  free all the memory used in this program

*/

void freeAll()
{
  free(nodes_vector);
  free(pageranksCurrent);
  free(pageranksPrevious);
  free(b);
  for(int i=0;i<nodes;i++)
    free(A[i]);
  free(A);
}

/*
  void chech_arguments()

  checks if the inputs are 2 and they are .txt files

*/

void check_arguments(int arguments,char **args)
{
  char *file1_extension,*file2_extension;
  if (arguments != 3)
  {
    printf("The arguments must be 3\n");
    printf("First argument (argv[1]): node file\n");
    printf("Second argument (argv[2]): adj list\n");
    exit(1);
  }
  file1_extension = strstr(args[1],".");
  file2_extension = strstr(args[2],".");
  if( strcmp(file1_extension,".txt")!= 0 || strcmp(file2_extension,".txt")!= 0)
  {
    printf("All input files must be .txt\n");
    exit(1);
  }
}

/*
  void createOutput()

  write to .txt file the nodes_vector

*/
 
void createOutput()
{
  char *line = NULL;
  size_t length = 0;
  int loops , return_value ;

  FILE *output = fopen("outputSequential.txt","w");
  if (output == NULL )
  {
    printf("Error opening the files : outputSequential.txt");
    exit(1);
  }
  for(int i=0;i<nodes;i++)
    fprintf(output,"%d\n",nodes_vector[i]);

  fclose(output);
}


void printA()
{
  printf("A MATRIX : \n");
  for(int i=0;i<nodes;i++)
  {
    for(int j=0;j<nodes;j++)
      printf("%f  ",A[i][j]);
    printf("\n");
  }
} 

void printNodesVector()
{
  printf("NODES : \n");
  for(int i=0;i<nodes;i++)
    printf("%d\n",nodes_vector[i]);
}

void print_bVector()
{
  printf("b VECTOR \n");
  for(int i=0;i<nodes;i++)
    printf("%0.15f\n",b[i]);
}

void printPageranks()
{
  printf("PAGERANKS :\n");
  for(int i=0;i<nodes;i++)
    printf("%0.15f\n",pageranksCurrent[i]);
}




