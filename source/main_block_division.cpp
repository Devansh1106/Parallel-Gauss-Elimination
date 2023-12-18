#include<iostream>
#include<fstream>
#include<cassert>
#include<cmath>
#include"mpi.h"


using namespace std;

int main(int argc,char** argv)
{
	double** A;
	double** LA;
	double* b;
	double* Lb;
	double* pivot_row;
	double pivot_row_rhs;

	double m,time_taken,time_taken_2;
	//MPI variables
	int rank, size,next;

	//Inputting of a matrix
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	int N = 96;

	time_taken = MPI_Wtime();		//record the start time

	LA = new double* [N/size];
	for (int i = 0; i < N/size; i++)
	{
		LA[i] = new double [N];
	}
	Lb = new double [N/size];
	pivot_row = new double [N];
	
	if (rank == 0)
	{
		A = new double* [N];
		for (int i = 0; i < N; i++)
		{
			A[i] = new double [N];
		} 

		ifstream read_file("input96.txt");
		assert(read_file.is_open());

		for(int i = 0; i < N; i++)
		{
			for (int j =0; j < N; j++)
			{
				read_file>>A[i][j];
			}
		}
		read_file.close();
		
		std::cout<<"Input done for matrix!"<<"\n";

		b = new double [N];
		ifstream read_file1("rhs96.txt");
		assert(read_file1.is_open());

		for (int i = 0; i < N; i++)
		{
			read_file1>>b[i];
		}
		read_file1.close();
		std::cout<<"Input done for rhs vector!";

		//Distributing matrix to all rank such that first "size" rows goes to rank 0 ans so on.

		for(int p=size-1; p >= 0; p--)
		{
			for(int i=p*(N/size); i < (N/size)*(p+1); i=i+1)
			{
				for(int j=0; j < N; j++)
				{
					LA[i-(p*(N/size))][j] = A[i][j];
				}
				Lb[i-(p*(N/size))] = b[i];
			}
		
			if(p!=0)
			{
				for(int i=0; i < N/size; i++)
				{
					MPI_Send(&LA[i][0],N,MPI_DOUBLE,p,10,MPI_COMM_WORLD);
					MPI_Send(&Lb[i],1,MPI_DOUBLE,p,20,MPI_COMM_WORLD);
				}
			}
		}
	}
	else
	{
		for(int i=0; i < N/size; i++)
		{
			MPI_Recv(&LA[i][0],N,MPI_DOUBLE,0,10,MPI_COMM_WORLD,&status);
			MPI_Recv(&Lb[i],1,MPI_DOUBLE,0,20,MPI_COMM_WORLD,&status);			
		}
	}

	//******************Forward Elimination starts*******************//

	for (int k=0; k < N; k++)
	{
		next = (k/(N/size));		
		if (rank == next && k < (N/size)*(rank+1))
		{	
			for (int p =next+1; p<size; p++)
			{
				MPI_Send(&LA[k-(rank*(N/size))][0],N,MPI_DOUBLE,p,50,MPI_COMM_WORLD);
				MPI_Send(&Lb[k-(rank*(N/size))],1,MPI_DOUBLE,p,60,MPI_COMM_WORLD);					
			}
		}
		else if (rank > next)
		{
			MPI_Recv(&pivot_row[0],N,MPI_DOUBLE,next,50,MPI_COMM_WORLD,&status);
			MPI_Recv(&pivot_row_rhs,1,MPI_DOUBLE,next,60,MPI_COMM_WORLD,&status);
		}
		if(rank!=next && rank > next)
		{
			for (int i=0; i < N/size; i++)  
			{
				m = LA[i][k]/pivot_row[k];
				for (int j = k; j < N; j++)
				{
					LA[i][j] = LA[i][j] - m*pivot_row[j];
				}
				Lb[i] = Lb[i] - m*pivot_row_rhs;
			}
		}
		else if(rank > next or rank == next)
		{
			for (int i =k+1; i < (N/size)*(rank+1); i++)
			{
				m = LA[i-(rank*(N/size))][k]/LA[k-(rank*(N/size))][k];
				for (int j = k; j < N; j++)
				{
					LA[i-(rank*(N/size))][j] = LA[i-(rank*(N/size))][j] - m*LA[k-(rank*(N/size))][j];
				}
				Lb[i-(rank*(N/size))] = Lb[i-(rank*(N/size))] - m*Lb[k-(rank*(N/size))];
			}
		}
		MPI_Barrier(MPI_COMM_WORLD); //Synchronising all ranks before gathering matrix
	}
	//******************Forward Elimination ends*******************//
	
	//Gathering parts of matrices
	if(rank == 0)
	{
		int i=0;
		while(i < (N/size))
		{
			for(int j=0; j < N; j++)
			{
				A[i][j] = LA[i][j];
			}
			b[i] = Lb[i];
			i++;
		}

		for (int p=1; p < size; p++)
		{
			for(int a=0; a < N/size; a++)
			{
				MPI_Recv(&LA[a][0],N,MPI_DOUBLE,p,30,MPI_COMM_WORLD,&status);
				MPI_Recv(&Lb[a],1,MPI_DOUBLE,p,40,MPI_COMM_WORLD,&status);
			}	
			
			while(i < (N/size)*(p+1))
			{
				for(int w=0; w < N/size; w++)
				{
					for(int j=0; j < N; j++)
					{
						A[i][j] = LA[w][j];
					}
					b[i] = Lb[w];
					i++;
				}
			}
		}
	}
	else
	{
		for(int a=0; a < N/size; a++)
		{
			MPI_Send(&LA[a][0],N,MPI_DOUBLE,0,30,MPI_COMM_WORLD);
			MPI_Send(&Lb[a],1,MPI_DOUBLE,0,40,MPI_COMM_WORLD);
		}		
	}
	MPI_Barrier(MPI_COMM_WORLD); //Synchronising all ranks before calculating solution

	//Calculating solution vector
	
	if(rank == 0)
	{
		double* sol;
		sol = new double [N];
		for ( int i=N-1; i > -1; i--)
		{
			for ( int j=N-1; j > i; j--)
			{
				b[i] = b[i] - sol[j]*A[i][j];
			}
			sol[i] = b[i]/A[i][i];
		}

		//Writing solution vector to a file
		
		ofstream write_output("Solution_block_division.txt");
		assert(write_output.is_open());
		for ( int i=0; i < N; i++)
		{
			write_output<<i<<" "<<sol[i]<<"\n";
		}
		write_output.close();
		std::cout<<"\n"<<"Output done!";
		time_taken_2 = MPI_Wtime();
		std::cout<<"\n"<<"Time taken = "<<time_taken_2 - time_taken; //record end time and print the difference

		//free memory
		delete [] sol;
		delete [] A;
		delete [] b;
	}
	delete [] LA;
	delete [] Lb;
	delete [] pivot_row;
	MPI_Finalize();
}
