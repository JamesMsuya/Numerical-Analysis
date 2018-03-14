#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


typedef struct coordinate_t{
	double x;
	double y;
}coordinate_t;


int dividedMethodInterpolation(coordinate_t *array,int size);

void solveLeastSquares(coordinate_t *matrix,int power, int size);

void freeMemory(double ** memory,int size);

void printMatrix(double ** matrix,int size);

void functionLoader(coordinate_t **array,int *size);

void functionLoader2(coordinate_t **array,int *size);

void solverGesb(double **matrix,int row, int column);

double max(double *matrix, int rows);

void copyMatrix( double ** matrixSource, double ** matrixDest, int rows,int column);

double ** solvedByGESP(double ** augmentedMatrix,int numberOfRows,int numberOfColumns);

double * backwardSub(double ** upperTriangularMatrix,int rows,int columns);

int applyPartialPivoting(double ** matrix,int rows,int columns,int pivoteIndex);






int main(int argc, char const *argv[]){
  	coordinate_t *array;
  	int size;
  
  	functionLoader(&array,&size);
	
	dividedMethodInterpolation(array,size);

	functionLoader2(&array,&size);

	//solveLeastSquares(array,2,size);


	return 0;
}







int dividedMethodInterpolation(coordinate_t *array,int size){
	int lcv=0;
	int i=0;
	int j=0;
	double **matrix = malloc(sizeof(double *)*size);
	for(lcv=0;lcv < size;lcv++ ){
		matrix[lcv]= malloc(sizeof(double)*(size));
		memset(matrix[lcv],0,sizeof(double)*size);
		matrix[lcv][0]= array[lcv].y;
	}
	

	for (i = 0 ; i < size; i++){
		for (j = 0; j < i; j++){
			matrix[i][j+1]= (matrix[i][j] - matrix[i-1][j])/ (array[i].x-array[i-j-1].x);
		}
	}
	printMatrix(matrix,size);
	printf("The polynomial is :\n");

	for (i = 0; i < size; ++i){
		printf("%f ",matrix[i][i]);
		for(j=0;j < i;j++){
			printf("(x- %f)",array[j]);

		}

		if (i!=size-1){
			printf(" + ");
		}
	}

	freeMemory(matrix,size);
	free(array);


}

void functionLoader(coordinate_t **array,int *size){
	int temp;

	temp=5;

	*array= malloc(sizeof(coordinate_t) * temp);

	(*array)[0].x=1.0, (*array)[0].y = 0.7651977;
	(*array)[1].x=1.3, (*array)[1].y = 0.6200860;
	(*array)[2].x=1.6, (*array)[2].y = 0.4554022;
	(*array)[3].x=1.9, (*array)[3].y = 0.2818186;
	(*array)[4].x=2.2, (*array)[4].y = 0.1103623;
	
	
	*size=temp;
}

void functionLoader2(coordinate_t **array,int *size){
	int temp;

	temp=5;

	*array= malloc(sizeof(coordinate_t) * temp);

	(*array)[0].x=0, (*array)[0].y = 1;
	(*array)[1].x=0.25, (*array)[1].y = 1.2840;
	(*array)[2].x=0.5, (*array)[2].y = 1.6487;
	(*array)[3].x=0.75, (*array)[3].y = 2.1170;
	(*array)[4].x=1, (*array)[4].y = 2.7183;
	
	
	*size=temp;
}


void printMatrix(double ** matrix,int size){
	int i = 0;
	int j = 0;
	for(i = 0; i < size; ++i){
		for(j = 0; j < size; ++j)
			fprintf(stdout, "%.5f ",matrix[i][j]);
		
		fprintf(stdout, "%c\n",'\n');
	}
}


void solveLeastSquares(coordinate_t *matrix,int power,int size){
	int i=0;
	int j=0;
	int k=0;
	int n=0;
	double **temp=malloc(size*sizeof(double *));
	for (i = 0; i < size; ++i){
		temp[i]=malloc(sizeof(double)*(size+1));
	}


for(i=0;i<power;i++){

	for (j = 0; j < size; j++){
		
		for(k=0;k<size;k++){

			temp[i][j] += pow(matrix[i].x,j+k);

		}

	}
	for (n = 0; i < size; ++n){
		temp[i][size]= matrix[i].y* pow(matrix[i].x,i);
	}
}


  	solverGesb(temp,size,size+1);













}











void solverGesb(double **matrix,int row, int column){
	int lcv=0;
	double * roots;
	double ** matrixCopy ;


	
	matrixCopy = solvedByGESP(matrix,row,column);


	if(matrixCopy!=NULL){

		roots=backwardSub(matrixCopy,row,column);
	}
	

	if(roots!=NULL){
		for (lcv = 0; lcv < row; lcv++){
			printf(" X%d = %f",lcv+1, roots[lcv]);
		}
		printf("\n");
	}
	freeMemory(matrixCopy,row);
	free(roots);	
}





int applyPartialPivoting(double ** matrix,int rows,int columns,int pivoteIndex){
	int lcv=0;
	int rowWithMaxRatio = pivoteIndex;
	double maxRatio = 0;
	double scaleFactor = 0 ;
	double tempMaxRatio = 0;
	double* swap = NULL;
	if(!(scaleFactor = max(&matrix[pivoteIndex][pivoteIndex],columns - pivoteIndex -1)))
		return 1;
	maxRatio = fabs(matrix[pivoteIndex][pivoteIndex])/ scaleFactor;	

	for(lcv = pivoteIndex+1; lcv < rows; lcv++){

		if(!(scaleFactor = max(&matrix[lcv][pivoteIndex],columns-pivoteIndex-1)))
			return 1;

		if((tempMaxRatio = fabs(matrix[lcv][pivoteIndex])/scaleFactor) > maxRatio)
		{
			maxRatio = tempMaxRatio;
			rowWithMaxRatio = lcv;
		}
	}

	if(rowWithMaxRatio != pivoteIndex){	
		swap = matrix[pivoteIndex];
		matrix[pivoteIndex] = matrix[rowWithMaxRatio];
		matrix[rowWithMaxRatio] = swap;
	}

	return 0 ;
}





double * backwardSub(double ** upperTriangularMatrix,int rows,int columns){
	int lcv=0;
	int lcv2;
	double sum = 0;
	double* roots = (double*) calloc(rows,sizeof(double));
	if(roots == NULL)
		return NULL;

	for(lcv = rows-1; lcv >= 0; --lcv ){
		sum = 0;
		for(lcv2 = lcv+1 ; lcv2 < columns-1; ++lcv2 )
		{
			sum+= upperTriangularMatrix[lcv][lcv2] * roots[lcv2];
		}

		roots[lcv] = (upperTriangularMatrix[lcv][columns-1] - sum) / upperTriangularMatrix[lcv][lcv];
	} 

	return roots;
}


double max(double *matrix,int rows){
	int lcv=0;
	double temp=fabs(matrix[0]);
		
		for (lcv = 0; lcv < rows; lcv++){
		
			if(temp < fabs(matrix[lcv])){
				temp=fabs(matrix[lcv]);
			}
		}

return temp;

}




double ** solvedByGESP(double ** augmentedMatrix,int numberOfRows,int numberOfColumns){
	int i=0;
	int j=0;
	int k=0;
	int error = 0;
	double multiplicand = 0;
	double pivot = 0;
	double ** matrixCopy  = (double**) calloc(numberOfRows,sizeof(double*));
	if(matrixCopy == NULL)
		return NULL;
	for(i = 0; i < numberOfRows; ++i)
	{
		matrixCopy[i] = (double*) calloc(numberOfColumns,sizeof(double));
		if(matrixCopy[i] == NULL)
		{
			freeMemory(matrixCopy,i);
			return NULL;
		}
	}

	copyMatrix(augmentedMatrix,matrixCopy,numberOfRows,numberOfColumns);

	for(i = 0; i < numberOfRows-1; ++i)
	{
		if(applyPartialPivoting(matrixCopy,numberOfRows,numberOfColumns,i))
		{
			freeMemory(matrixCopy,numberOfRows);
			return NULL;
		}
		
		pivot = matrixCopy[i][i];

		if(pivot!=0)
		{
			for(j = i+1; j < numberOfRows; ++j)
			{
				multiplicand = (matrixCopy[j][i]/pivot);

				for(k = i; k < numberOfColumns; ++k)
				{
					matrixCopy[j][k] -= matrixCopy[i][k] * multiplicand;
				}				
			}
		}
	}

	for(i = 0; i < numberOfRows; ++i)
	{
		error = 1;
		for(j = 0; j < numberOfColumns-1; ++j)
		{
			if(matrixCopy[i][j] != 0)
				error = 0;
		}

		if(error)
		{
			freeMemory(matrixCopy,numberOfRows);
			return NULL;
		}
	}

	return matrixCopy;
}


void copyMatrix( double ** matrixSource, double ** matrixDest, int rows,int column){
	int lcv=0;
	int lcv2=0;
	for(lcv = 0; lcv < rows; lcv++){
		for(lcv2 = 0; lcv2 < column; lcv2++){
			matrixDest[lcv][lcv2] = matrixSource[lcv][lcv2];
		}
	}
}


void freeMemory(double ** memory,int size){
	int i=0;
	for(i = 0; i < size; ++i)
	{
		free(memory[i]);
	}

	free(memory);
}