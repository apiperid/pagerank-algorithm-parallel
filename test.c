/*
  PIPERIDIS ANESTIS 
  AEM : 8689
  PARALLEL AND DISTRIBUTED
  PAGERANK - TESTING PROGRAM 
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*

  THIS test.c FILE CHECKS IF 2 .txt FILES ARE THE SAME
  USEFUL IN ORDER TO CHECK IF THE PARALLEL ALGORITHM HAS THE SAME 
  RESULT IN COMPARISON WITH THE RESULT OF THE SEQUENTIAL ALGORITHM
*/

int main(int argc,char **argv)
{
  char *filename1,*filename2;
  FILE *file1,*file2;
  int file1_lines = 0,file2_lines = 0,same_lines = 0;
  char *line = NULL,*line2 = NULL;
  size_t len = 0,len2 = 0;

  // 2 .txt FILES MUST BE AS INPUT

  if(argc!=3)
  {
    printf("Error ..... arguments must be 3\n");
    exit(1);
  }

  filename1 = argv[1];
  filename2 = argv[2];

  if( strcmp(strstr(filename1,"."),".txt")!= 0 || strcmp(strstr(filename2,"."),".txt")!= 0)
  {
    printf("ERROR ..... All input files must be .txt\n");
    exit(1);
  }
  
  file1 = fopen(filename1,"r");
  file2 = fopen(filename2,"r");
  if(file1==NULL || file2==NULL )
  {
    printf("Error in opening the following files:\n");
    printf("%s\n",filename1);
    printf("%s\n",filename2);
    exit(1);
  }

  // CHECK IF THE FILES HAVE THE SAME NUMBER OF NODES

  while ( getline(&line, &len, file1)!= -1 )
    file1_lines ++ ;
  
  while ( getline(&line, &len, file2)!= -1 )
    file2_lines ++ ;

  if(file1_lines != file2_lines )
  {
    printf("The files do not have the same number of lines\n");
    exit(1);
  }
 
  rewind(file1);
  rewind(file2);

  //printf("%d   %d \n",file1_lines,file2_lines);

  //CHECK IF THE RESULT IS THE SAME IN THE 2 FILES

  while ( (getline(&line, &len, file1)!= -1) && (getline(&line2, &len2, file2)!= -1))
    if ( strcmp(line,line2) == 0)
      same_lines++;

  // PRINT THE PERCENTAGE OF THE SUCCESS

  printf("The Success Rate is %0.4f percent\n",((double)same_lines/file1_lines)*100);

  fclose(file1);
  fclose(file2);
}
