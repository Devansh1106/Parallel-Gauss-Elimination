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

		
		// if(rank ==0)
		// {
		// 	cout<<"\n";
		// 	for(int i=0; i < N; i++)
		// 	{
		// 		for(int j=0; j < N; j++)
		// 		{
		// 			cout<<" "<<A[i][j];
		// 		}	
		// 		cout<<"\n";
		// 	}		
		// }	


		b = new double [N];
		ifstream read_file1("rhs96.txt");
		assert(read_file1.is_open());

		for (int i = 0; i < N; i++)
		{
			read_file1>>b[i];
		}
		read_file1.close();
		std::cout<<"Input done for rhs vector!";

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
			/*
			if(p==3)
			{
				for(int i=0;i<N/size; i++)
				{
					for(int j=0; j<N;j++)
					{
						cout<<" "<<LA[i][j];
					}
					cout<<"\n";
				}
			}*/
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

	//************************************************************************//
	// if(rank ==7)
	// {
	// 	//cout<<"\n";
	// 	for(int i=0; i < N/size; i++)
	// 	{
	// 	// 	for(int j=0; j < N; j++)
	// 	// 	{
	// 	// 		cout<<" "<<LA[i][j];
	// 	// 	}
	// 		//cout<<"\n";
	// 		cout<<Lb[i]<<"\n";
	// 	}		
	// }


	//******************Forward Elimination*******************//

	for (int k=0; k < N; k++)
	{
		// if (k ==5)
		// {
		// 	cout<<"2";
		// }
		//for(int i=0; i < N; i++)
		//{
			next = (k/(N/size));
			// if(rank ==1)
			// {
			// 	cout<<k<<"\n";
			// }
			
			if (rank == next && k < (N/size)*(rank+1))
			{	
				// for(int i=0; i<N; i++)
				// {
				// 	pivot_row[i] = LA[k-(rank*size)][i];
				// }
				// pivot_row_rhs = &Lb[k-(rank*size)];
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

		//}
		// if( rank ==3)
		// {
		// 	for(int i=0; i<N;i++)
		// 		{
		// 			cout<<pivot_row[i]<<" ";
		// 		}
		// 	cout<<"\n";
		// }
		// break;
		// if(rank == 6)
		// {
		// 	cout<<"\n"<<pivot_row_rhs;
		// }
		// if (rank ==1)
		// {
		// 	cout<<k<<"\n";
		// }
		if(rank!=next && rank > next)
		{
			// if(rank ==2)
			// {
			// 	cout<<"rank 2";
			// }
			// if (rank == 3)
			// {
			// 	cout<<"rank 3";
			// }
			for (int i=0; i < N/size; i++)  
			{
				// if(pivot_row[k] == 0)
				// {
				// 	cout<<k<<"\n";
				// }
				m = LA[i][k]/pivot_row[k];
				for (int j = k; j < N; j++)
				{
					LA[i][j] = LA[i][j] - m*pivot_row[j];
					// if(rank == 7)
					// {
					// 	cout<<" "<<LA[i][j]<<" ";
					// }
					//cout<<"\n";
				}
				// cout<<"\n";
				Lb[i] = Lb[i] - m*pivot_row_rhs;
			}
			// if(rank == 7)
			// {
			// 	cout<<m<<" "<<LA[i][j]<<"    ";
			// }
		}
		else if(rank > next or rank == next)
		{
			for (int i =k+1; i < (N/size)*(rank+1); i++)
			{
				// if(LA[k-(rank*(N/size))][k]== 0)
				// {
				// 	cout<<"t"<<k<<"\n";
				// }
				m = LA[i-(rank*(N/size))][k]/LA[k-(rank*(N/size))][k];
				for (int j = k; j < N; j++)
				{
					LA[i-(rank*(N/size))][j] = LA[i-(rank*(N/size))][j] - m*LA[k-(rank*(N/size))][j];
					// if (i==15 and k==14)
					// {
					// 	cout<<LA[i][k]<<" ";
					// }
				}
				//cout<<"\n";
				Lb[i-(rank*(N/size))] = Lb[i-(rank*(N/size))] - m*Lb[k-(rank*(N/size))];
				// if (rank ==7)
				// {
				// 	cout<<m<<" "<<Lb[i-(rank*(N/size))]<<"\n";
				// }
			}
			// if(rank == 7)
			// {
			// 	cout<<LA[0][15]<<"\n";
			// }
		}
		// if (rank ==0)
		// {
		// 	cout<<k<<"\n";
		// }
		MPI_Barrier(MPI_COMM_WORLD);
		
	}
	// if(rank ==0)
	// {
	// 	cout<<"rank 0";
	// }
	// if(rank ==1)
	// {
	// 	cout<<"rank 1";
	// }
	// if(rank ==2)
	// {
	// 	cout<<"rank 2";
	// }
	// if(rank ==3)
	// {
	// 	cout<<"rank 3";
	// }


	// if (rank==7)
	// {
	// 	cout<<"\n";
	// 	for(int i=0; i < N/size; i++)
	// 	{
	// 		for(int j=0; j < N; j++)
	// 		{
	// 			cout<<" "<<LA[i][j];
	// 		}	
	// 		// cout<<" "<<b[i];
	// 		cout<<"\n";
	// 	}
	// }
		//******************Forward Elimination End*******************//
			
		// for (int k=0; k < N/size; k++)
		// {
		// 	//forward elimination
		// 	//making elements 0 below lower part of column k
		// 	for (int i=k+1; i < N/size; i++)
		// 	{
		// 		m = LA[i][k]/LA[k][k];
		// 		for (int j = k; j < N; j++)
		// 		{
		// 			LA[i][j] = LA[i][j] - m*LA[k][j];
		// 			//cout<<LA[i][j]<<" ";
		// 		}
		// 		//cout<<"\n";
		// 		Lb[i] = Lb[i] - m*Lb[k];
		// 	}
		// }

	//***********************************************//
	// if(rank ==1)
	// {
	// 	cout<<"\n";
	// 	for(int i=0; i < N/size; i++)
	// 	{
	// 		for(int j=0; j < N; j++)
	// 		{
	// 			cout<<" "<<LA[i][j];
	// 		}	
	// 		cout<<"\n";
	// 	}		
	// }

	//MPI_Barrier(MPI_COMM_WORLD);

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
			// if(p==2)
			// {
			// 	for (int z=0; z < N; z++)
			// 	{

			// 		cout<<LA[1][z]<<" ";
			// 	}
			// }		
			
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
		// if(rank == 7)
		// {
		// 	for(int c=0; c < N/size; c++)
		// 	{
		// 		// for(int d=0; d < N; d++)
		// 		// {
		// 		// 	cout<<LA[c][d]<<" ";
		// 		// }
		// 		cout<<Lb[c]<<"\n";
		// 	}
		// }			
	}
	// if(rank ==0)
	// {
	// 	for(int z=0; z < N;z++)
	// 	{
	// 		for(int n=0;n< N;n++)
	// 		{
	// 			if(fabs(A[z][n])<1e-10)
	// 			{
	// 				A[z][n] = 0.0;
	// 			}
	// 			if(fabs(b[z])<1e-10)
	// 			{
	// 				b[z] = 0.0;
	// 			}
	// 		}
	// 	}
		// for(int z=0; z< N;z++)
		// {
		// 	cout<<b[z]<<"\n";
		// }
		// for(int z=0; z < N;z++)
		// {
		// 	for(int n=0;n< N;n++)
		// 	{
		// 		cout<<A[z][n]<<" ";
		// 	}
		// 	cout<<"\n";
		// }
	// }
	// if(rank ==0)
	// {
	// 	cout<<"\n";
	// 	for(int i=0; i < N; i++)
	// 	{
	// 		// for(int j=0; j < N; j++)
	// 		// {
	// 		// 	cout<<" "<<A[i][j];
	// 		// }	
	// 		cout<<" "<<b[i];
	// 		cout<<"\n";
	// 	}		
	// // //cout<<"\n"<<A[4][2];
	// }
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank == 0)
	{
		// for(int i=0;i<N;i++)
		// {
		// 	cout<<b[i]<<"\n";
		// }
		//cout<<A[95][95];
		double* sol;
		sol = new double [N];
		//sol[N-1] =0;
		for ( int i=N-1; i > -1; i--)
		{
			// if(A[i][i] == 0)
			// {
			// 	cout<<"\n"<<i;
			// }
			for ( int j=N-1; j > i; j--)
			{
				// if (A[i][i] == 0)
				// {
				// 	cout<<i<<";";
				// }
				b[i] = b[i] - sol[j]*A[i][j];
			}
			sol[i] = b[i]/A[i][i];
		}
		// for (int i =0;i<N;i++)
		// {
		// 	cout<<"\n"<<sol[i]<<"\n";
		// }
		//cout<<sol[998]<<"\n";
		ofstream write_output("Solution_block_division.txt");
		assert(write_output.is_open());
		for ( int i=0; i < N; i++)
		{
			write_output<<i<<" "<<sol[i]<<"\n";
		}
		write_output.close();
		std::cout<<"\n"<<"Output done!";
		time_taken_2 = MPI_Wtime();
		std::cout<<"\n"<<"Time taken = "<<time_taken_2 - time_taken;
		delete [] sol;
		delete [] A;
		delete [] b;

	}
	delete [] LA;
	delete [] Lb;
	delete [] pivot_row;

	
	MPI_Finalize();
}