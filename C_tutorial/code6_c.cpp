#include <stdio.h>
#include <stdlib.h> 


struct myArray
{
  int len;
  double * data;
};



void fillArr(myArray * arr)
{
  for (int i = 0; i < arr->len; ++i)
  {
    arr->data[i] = i; //attention: "->" instead of "."
  }

}



int main(int argc, char** argv) 
{
   
  myArray arr;
  arr.len = 10;
  arr.data = (double *)malloc (arr.len * sizeof(double));

  fillArr(&arr);
  
  //how to print?

  return 0;
}







// printArr(&arr);

// void printArr(myArray * arr)
// {
//   for (int i = 0; i < arr->len; ++i)
//   {
//     printf("%f\n", arr->data[i]); 
//   }

// }