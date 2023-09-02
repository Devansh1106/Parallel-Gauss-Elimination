# Parallel-Gauss-Elimination #
This is code for performing Gauss elimination in two ways:  
1- Using parallel algorithm such as Messsage Passing Interface (MPI).  
2- Using serial algorithm.  

## Structure of the input and output files ##
This code takes only square matrix as input. All files will be in .txt format.  
Structure of the input files will:   
- Input matrix file will contain a square matrix of size NxN where N is the order of the matrix.
- Input right hand side (rhs) file will contain a vector of size Nx1 where N is the order of the matrix.
  
Struture of the output file will be:
- Output file will contain a solution in form of vector of size Nx1 where N is the order of the matrix. 

## Configuration of code for your own input ##
Four changes are required for changing the input and congiure the code according to your own input:  

1- First change the order of the matrix to the new order that you are going to input (N).  
2- Change the name of the input file to the name of the file that you are going input.  
3- Change the name of the right hand side file (rhs) to the new file that you are going to input.  
4- Change the name of the output file in the code to the name that you preferred.  

## Command for compiling code ##
First you need to install g++ and MPICH for compiling `main_serial.cpp` and `main_block_divison.cpp` respectively. Instructions are for using from terminal.  
- For compiling main_serial.cpp file
    - Use `g++ main_serial.cpp`
    - For running `.exe`, use `./<exe_file_name>`
- For compiling main_block_division.cpp file
    - Use `mpic++ main_block_division.cpp`
    - For running `.exe`, use `mpirun -np <Number_of_process> <exe_file_name>`
