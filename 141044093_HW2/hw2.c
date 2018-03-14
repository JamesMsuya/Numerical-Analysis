#include <stdio.h>
#include <math.h> 
#include <stdlib.h>
#include <string.h>


void solverJCB(double **matrix,int row, int column);

void solverGesb(double **matrix,int row, int column);

double max(double *matrix, int rows);

void copyMatrix( double ** matrixSource, double ** matrixDest, int rows,int column);

double ** solvedByGESP(double ** augmentedMatrix,int numberOfRows,int numberOfColumns);

double * backwardSub(double ** upperTriangularMatrix,int rows,int columns);

int applyPartialPivoting(double ** matrix,int rows,int columns,int pivoteIndex);

void fileReader(const char *file,double ***matrix, int* rows, int* column);

void freeMatrix(double **matrix,int rows);

void freeMemory(double ** memory,int size);

void printMatrix(double ** matrix,int rows,int columns);



int main(int argc, const char **argv){

	int column;
	int rows;
	double **matrix;

 	fileReader(argv[2],&matrix,&rows,&column);


 	if(strcmp(argv[4],"JCB")==0){

		solverJCB(matrix,rows,column);

	}
	if(strcmp(argv[4],"GESP")==0){

		solverGesb(matrix,rows,column);

	}
	else{
		printf("Usage ./slover -i system.txt -m METHOD \n");
	}
	freeMatrix(matrix,rows);
	return 0;
}



void solverJCB(double **matrix,int row, int column){
	int lcv=0;
	int lcv2=0;
	int lcv3=0;
	int count=0;
	double stop=0.1;
	double *matrixSol=malloc(sizeof(double)*column-1);
	double *matrixSol1=malloc(sizeof(double)*column-1);
	double *matrixSol2=malloc(sizeof(double)*column-1);
	memset(matrixSol,0,sizeof(double)*column-1);
	memset(matrixSol1,0,sizeof(double)*column-1);
	memset(matrixSol2,0,sizeof(double)*column-1);

	printMatrix(matrix,row,column);

	if(column-1 > row){

		printf("The System has no unique Solution\n");
		return;
	}

	while(stop > 0.001){
		count++;
		for(lcv2=0;lcv2<row;lcv2++){

			for (lcv = 0; lcv < column-1; lcv++){
				if(lcv2!=lcv){

					matrixSol1[lcv2] += -(matrixSol[lcv] * (matrix[lcv2][lcv])/(matrix[lcv2][lcv2]));
				}

			}

			matrixSol1[lcv2] += matrix[lcv2][lcv]/(matrix[lcv2][lcv2]);;
		}

		for(lcv3=0;lcv3<column-1;lcv3++){
			matrixSol2[lcv3]=matrixSol1[lcv3]-matrixSol[lcv3];
			matrixSol[lcv3]=matrixSol1[lcv3];
		}

		stop=max(matrixSol2,row)/max(matrixSol1,row);
		memset(matrixSol1,0,sizeof(double)*column-1);

		if(count>100){
			printf("The system of Equations does not converge i.e p(A) > 1.\n");
			free(matrixSol);
			free(matrixSol1);
			free(matrixSol2);
			freeMatrix(matrix,row);
			exit(0);
		}
	}
	printf("Solution Set\n");
	for(lcv3=0;lcv3 < column-1;lcv3++){
		printf("x%d = %f ",lcv3+1,matrixSol[lcv3]);
	}
	printf("\n");

free(matrixSol);
free(matrixSol1);
free(matrixSol2);
}



void solverGesb(double **matrix,int row, int column){
	int lcv=0;
	double * roots;
	double ** matrixCopy ;

	if(column < row-1){

		printf("The System has no unique Solution\n");
		return;
	}
	
	printMatrix(matrix,row,column);
	
	printf("\n");

	
	matrixCopy = solvedByGESP(matrix,row,column);

	printf("%6c AUGUMENTED MATRIX \n",' ');

	printf("-----------------------------\n");

	printMatrix(matrixCopy,row,column);

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



void fileReader(const char *file,double ***matrix,int* rows,int* column){
	FILE * input = fopen(file,"r");
	char r='`';
	int lcv=0;
	int lcv1=0;
	int temp=0;
	double num=0;

	if(input==NULL){
		printf("File could not be opened\n");
		exit(0);
	}

	while(!feof(input)){
		r='`';
		while(r!='\n' && !feof(input)){
			temp=fscanf(input, "%lf%*[ ]", &num);
			if(temp!=EOF){
				lcv++;
			}
			fscanf(input,"%c",&r);
			if(temp==EOF){
				lcv1-=1;
				break;
			}
			
		}
		++lcv1;
	}

	*column = lcv/lcv1;
	*rows = lcv1;
	
	*matrix=(double **)malloc(sizeof(double *)*(*rows));
	for(lcv=0;lcv < (*rows);lcv++ ){
		(*matrix)[lcv]=(double *)malloc(sizeof(double)*(*column));

	}
	fclose(input);
	input=fopen(file,"r");

	lcv=0;
	lcv1=0;

	
	while(!feof(input)){
		r='`';
		while(r!='\n' && !feof(input)){
			temp=fscanf(input, "%lf%*[ ]", &num);
			if(temp!=EOF){
				lcv++;
			}
			fscanf(input,"%c",&r);
			
			if(temp==EOF){
				lcv1-=1;
				break;
			}
			(*matrix)[lcv1][(lcv-1)%(*column)]=num;
		}
		lcv1++;
	}
fclose(input);
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



void freeMatrix(double **matrix,int rows){
	int i=0;
	for(i=0; i < rows;i++){
		
		free((matrix)[i]);
	}
	free(matrix);

}
void printMatrix(double ** matrix,int rows,int columns){
	int i = 0;
	int j = 0;
	for(i = 0; i < rows; ++i){
		for(j = 0; j < columns; ++j)
			fprintf(stdout, "%.5f ",matrix[i][j]);
		
		fprintf(stdout, "%c\n",'\n');
	}
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