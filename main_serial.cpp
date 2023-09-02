#include<iostream>
#include<fstream>
#include<cmath>
#include<cassert>
#include<chrono>
using namespace std;

int main()
{
	double** A;
	//double** B;
	double* b;
	double* sol;
	double m;
	int p_row;
	double temp,temp1,max;
	int N=5000;
	auto start = chrono::steady_clock::now(); //for execution time 


	//Inputting of a matrix

	A = new double* [N];
	for (int i = 0; i < N; i++)
	{
		A[i] = new double [N];
	} 

	ifstream read_file("input5000.txt");
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

	//Copying the inputted matrix
	// B = new double* [3];
	// for (int i = 0; i < 3; i++)
	// {
	// 	B[i] = new double [3];
	// } 
	// for (int i = 0; i < 3; i++)
	// {
	// 	for (int j =0; j < 3; j++)
	// 	{
	// 		B[i][j] = A[i][j];
	// 	}
	// }	

	//Inputting rhs vector

	b = new double [N];
	ifstream read_file1("rhs5000.txt");
	assert(read_file1.is_open());

	for (int i = 0; i < N; i++)
	{
		read_file1>>b[i];
	}
	read_file1.close();
	std::cout<<"Input done for rhs vector!";

	//forward elimination
	for (int k=0; k < N; k++)
	{
		for (int j=k; j < N; j++)
		{
			max = 0.0;
			if (max < fabs(A[j][k]))
			{
				p_row = j;
			}
		}

		//pivot if required
		if (p_row != k)
		{
			//exchanging p_row with kth row
			for (int j=0; j < N; j++)
			{
				temp1 = A[k][j];
				A[k][j] = A[p_row][j];
				A[p_row][j] = temp1;
			}
			temp = b[k];
			b[k] = b[p_row];
			b[p_row] = temp;
		}
		
		//making elements 0 below lower part of column k
		for (int i=k+1; i < N; i++)
		{
			m = A[i][k]/A[k][k];
			for (int j = k; j < N; j++)
			{
				A[i][j] = A[i][j] - m*A[k][j];
			}
			b[i] = b[i] - m*b[k];
		}
	}

	//Printing the RREF form

	// for (int i=0; i < 3; i++)
	// {
	// 	for (int j=0; j < 3; j++)
	// 	{
	// 		cout<<A[i][j]<<"\t";
	// 	}
	// 	cout<<"\n";
	// }
	// for ( int i=0; i < 3; i++)
	// {
	// 	cout<<b[i]<<"\t";
	// }

	sol = new double [N];
	//sol[N-1] =0;
	for ( int i=N-1; i > -1; i--)
	{
		for ( int j=N-1; j > i; j--)
		{
			b[i] = b[i] - sol[j]*A[i][j];
		}
		sol[i] = b[i]/A[i][i];
	}
	ofstream write_output("Solution_serial.txt");
	assert(write_output.is_open());
	for ( int i=0; i < N; i++)
	{
		write_output<<i<<" "<<sol[i]<<"\n";
	}
	write_output.close();
	std::cout<<"\n"<<"Output done!";
	// cout<<A[15][15]<<"\n";
	// cout<<b[15]<<"\n";
	auto end = chrono::steady_clock::now();  //for execution time
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
	delete [] A;
	delete [] b;
	delete [] sol;
	//delete [] B;
}
