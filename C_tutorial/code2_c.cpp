#include <stdio.h>
 

int main(int argc, char** argv) {
   
  //The first argument is always the name of the executable

  printf("There are %d arguments.\n",argc-1);


  if(argc>1)
  {
  	printf("The arguments are:\n");
  	for (int i = 0; i < argc; ++i)
  	{
  		printf("\t%s\n",argv[i]);
  	}

  }

  return 0;
}
