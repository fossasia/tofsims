//#include <Rdefines.h>
#include <Rcpp.h>
using namespace Rcpp;

typedef struct {
  float *array;
  size_t used;
  size_t size;
} Array;

int readRawPhi(char* filename, 
               Array *rawdata, 
               float slope, 
               float intercept, 
               int imagepixels);
void initArray(Array *a, size_t initialSize);
void insertArray(Array *a, float element);
void freeArray(Array *a);
unsigned long readCertainBits(unsigned long raw, 
                              int numOfBits, 
                              int readFrom, 
                              int readTo);

// [[Rcpp::export]]
NumericMatrix readRawPhiC(CharacterVector rFilename, 
                          float rSlope, 
                          float rIntercept, 
                          float rImagePixels) {
  char* filename = rFilename[0];
  int x,y;
  x = 0;
  y = 0;
  Array rdata[3];
  int result = readRawPhi(filename, 
                          &rdata[0], rSlope, rIntercept, rImagePixels);
  if(result == 0){
    NumericMatrix Rmatrix(x,y);
    return(Rmatrix);
  }
  x = rdata[0].used;
  y = 2;
  NumericMatrix Rmatrix(x,y);
  int i, j;
  for(i = 0; i < x; i++){
    for(j = 0; j < y; j++){
      Rmatrix[i + x * j] = rdata[j].array[i];
    }
  }
  return Rmatrix;
}

unsigned long long readCertainBits(unsigned long long raw, 
                                   int numOfBits, 
                                   int readFrom, 
                                   int readTo){
  unsigned long long result;
  int shift_number = numOfBits - readTo;
  result = (unsigned long long)pow(2, (readTo - readFrom + 1)) - 1;
  result = result << shift_number;
  result = raw & result;
  result = result >> shift_number;
  return result;
}

int readRawPhi(char* filename, 
               Array *rawdata, 
               float slope, 
               float intercept, 
               int imagepixels)
{
  //  FILE *lf = fopen("log.log", "wb");
  FILE *f = fopen(filename, "rb");
  if(f == NULL){
    //printf("%s", filename);

    perror("Error");
        // printf("\nUnable to open file.");
    return(0);
  }
  
  fseek(f, 0, SEEK_END);
  int fileLen = ftell(f);
  fseek(f, 0, SEEK_SET);
  
  initArray(&rawdata[0], fileLen);
  initArray(&rawdata[1], fileLen);
  initArray(&rawdata[2], fileLen);
  
  char *header;
  header = (char *)malloc(4096 * sizeof(char));
  int readData = fread(header, sizeof(char), 4096, f);
  
  unsigned long long *raw;
  raw = (unsigned long long *) malloc (1 * sizeof(unsigned long long));
  unsigned short *last_in_id;
  last_in_id = (unsigned short *) malloc (1 * sizeof(unsigned short));
  int i = 0;
  Array last_in_ids;
  initArray(&last_in_ids, 8);
  
  while(!feof(f)){
    // Read last in ids
    for(;i < 8; i++){
      readData = readData && fread(last_in_id, 
                                   sizeof(unsigned short), 1, f);
      insertArray(&last_in_ids, *last_in_id);
    }
    // End read last in id
    // Read blocks
    int block_length = (int)(last_in_ids.array[last_in_ids.used - 2] / 8);
    
    i = 0;
    for(; i < block_length; i++){
      readData = fread(raw, sizeof(unsigned long long), 1, f);
      unsigned long long p1 = readCertainBits(*raw, 64, 1, 10);
      unsigned long long p2 = readCertainBits(*raw, 64, 11, 21);
      unsigned long long p3 = readCertainBits(*raw, 64, 22, 32);
      unsigned long long p4 = readCertainBits(*raw, 64, 33, 37);
      unsigned long long p5 = readCertainBits(*raw, 64, 38, 57);
      // unsigned long long p6 = readCertainBits(*raw, 64, 58, 64);
      // fprintf(lf, "%d - %lu,%lu,%lu,%lu,%lu\n", i, p1, p2, p3, p4, p5);
      if(p1 != 0 || p4 == 21)
        continue;
      //       if(p5 == 95268269 || p1 != 0)
      //           continue;
      float mz = pow(((p5-intercept)/slope),2);
      // float mz = p6;
      insertArray(&rawdata[0], p2*imagepixels + p3 + 1);
      insertArray(&rawdata[1], mz);
      insertArray(&rawdata[2], 0);
      //fprintf(lf, "%d,%d,%d,%d,%d,%d\n", p1, p2, p3, p4, p5, p6);
    }
    //    free(raw);
    // End read blocks
  }
  
  //  freeArray(&last_in_ids);
  // fclose(lf);
  fclose(f);
  return(1);
}


void initArray(Array *a, size_t initialSize) {
  a->array = (float *)malloc(initialSize * sizeof(float));
  a->used = 0;
  a->size = initialSize;
}

void insertArray(Array *a, float element) {
  if (a->used == a->size) {
    a->size += 1;
    a->array = (float *)realloc(a->array, a->size * sizeof(float));
  }
  a->array[a->used++] = element;
}

void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}
