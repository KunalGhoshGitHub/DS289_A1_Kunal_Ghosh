#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include<functional>
#include"DS289.h"

using namespace std;

// Using the double as a real number
typedef double real;

void Analytical_Solution(vector<real> &T, vector<real> &r, real S,real C1, real C2)
{
    for(int i = 0;i < r.size();i++)
    {
        T.push_back(-( ((1.0/4.0)*S*r[i]*r[i]) + (C1*log(r[i])) + C2));
    }
}

int main()
{
    // ODE: T'' + (T'/r) + S = 0
    // Interval R0 <= r <= R
    // Interval starts at: R0
    // Interval ends at: R
    // dr = Step size
    // S is the quantity in the ODE: T'' + (T'/r) + S = 0
    real dr,R,R0,S;

    // Vector to read from the input files
    vector<real> input, Values_of_S;

    // Vector to store the names of all the output files generated
    vector<string> Output_file_names;

    // "Input.txt" has all the parameters except the value of the S(s)
    input = input_parameters("Input.txt");

    // "Input_Values_of_S.txt" has the values for different S(s)
    Values_of_S = input_parameters("Input_Values_of_S.txt");


    R0 = input[1];
    R = input[2];

    // num_interval = Number of intervals
    // num_grid_points = Number of grid points
    int num_interval,num_grid_points;

    num_grid_points = input[0];

    // dr = Step size
    dr = real((R-R0)/(num_grid_points-1));

    num_interval = num_grid_points-1;
    // r = Location of the grid points
    vector<real> r(num_grid_points);

    // Function t_interval(address to the vector of the grid points (r), step size(dr), start of interval(R0), end of interval(R))
    // Writes the location of the grid points to the vector of the grid points (r)
    t_interval(r,dr,R0,R);

    // Number of the rows and columns in the matrix
    int rows ,cols;

    // Boundary condition is applied at the first and last point
    // Hence, the size reduced by 2
    rows = num_grid_points - 2;
    cols = num_grid_points - 2;

    // System of equation AX = B
    // Declaring a vector of vectors to use it as a matrix
    vector<vector <real> > Matrix(rows,vector<real> (cols,0.0));

    real dtdr,Tn;

    // Boundary conditions
    // Boundary Condition: T (r = 1) = Tn
    Tn = input[3];
    // Boundary Condition: dT /dr at (r = 0) = dtdr
    dtdr = input[4];

    // This for loop will iterate over different values of the S(s)
    for (int i=0;i<Values_of_S.size();i++)
    {
        // Declaring the vectors X and B
        vector<real> B(num_grid_points-2);
        vector<real> X(num_grid_points-2);

        // Value of S
        S = Values_of_S[i];

        // We are populating the first row and last manually while putting the boundary conditions

        // First wow
        Matrix[0][0] = -2.0;
        Matrix[0][0+1] = 2.0;
        B[0] = -(S*dr*dr) - (dtdr*2*r[1]*(1-(dr/(2*r[1]))));

        // Last Row
        Matrix[rows-1][rows-1] = -2;
        Matrix[rows-1][rows-2] = (1 - (dr/(2*r[rows-1])));
        B[rows-1] = -(S*dr*dr) - (Tn*(1 + (dr/(2*r[rows-1]))));

        // We are assigning the values  of non-zero elements of the matrix except for first and last row
        for (int row = 1;row < rows-1;row++)
        {
            // Only assign values to NON ZERO elements of the Matrix
            Matrix[row][row-1] = (1 - (dr/(2*r[row])));
            Matrix[row][row] = -2;
            Matrix[row][row+1] = (1 + (dr/(2*r[row])));
        }

        // We are assigning the values  of non-zero elements of B except for first and last row
        for (int row = 1;row<rows-1;row++)
        {
            B[row] = -S*dr*dr;
        }

        // Function Thomas_Algorithm_Matix_Solver(A,B,X) will solve the system of equation AX = B
        // where, A is a tridiagonal matrix
        Thomas_Algorithm_Matix_Solver(Matrix,B,X);
        //X[0] = X[2];
        //X[X.size()-1] = Tn;

        X.insert(X.begin(),X[1]);
        X.push_back(Tn);

        // Writes the grid point location in a file
        write_to_file(r,"Q5_a_Grid_Points_Solution_S_"+to_string(S)+".csv");

        Output_file_names.push_back("Q5_a_Grid_Points_Solution_S_"+to_string(S)+".csv");

        // Writes the numerical solution using Thomas Algorithm at different grid point location in a file
        write_to_file(X,"Q5_a_Solution_S_"+to_string(S)+".csv");

        Output_file_names.push_back("Q5_a_Solution_S_"+to_string(S)+".csv");

        // Vector corresponding to analytical solution
        vector<real> Analytical_T;

        real C1,C2;

        C1 = 0.0;

        C2 = -(1.0+(S/4.0));

        Analytical_Solution(Analytical_T,r,S,C1,C2);

        write_to_file(Analytical_T,"Q5_Analytical_Solution_S_"+to_string(S)+".csv");

        Output_file_names.push_back("Q5_Analytical_Solution_S_"+to_string(S)+".csv");

        // Vector corresponding to analytical solution
        vector<real> Error_T;

        Error(Analytical_T,X,Error_T);

        write_to_file(Error_T,"Q5_Error_S_"+to_string(S)+".csv");

        Output_file_names.push_back("Q5_Error_S_"+to_string(S)+".csv");

    }

    write_to_file(Output_file_names,"Output_file_names.csv");

    real T_req,step_s,S_0,S_n;
    vector<real> S_req,T_max,S_optimum,input_optimum;

    // "Input_For_Optimum_S.txt" has all the parameters to optimize the value of S
    input_optimum = input_parameters("Input_For_Optimum_S.txt");
    S_0 = input_optimum[0];
    S_n = input_optimum[1];
    step_s = input_optimum[2];
    T_req = input_optimum[3];

    for (int i = 0; S_0 +(i*step_s) < S_n;i++)
    {
        S_req.push_back(S_0+(i*step_s));
    }

    T_req = 100;

    for (int i=0;i<S_req.size();i++)
    {
        // Declaring the vectors X and B
        vector<real> B(num_grid_points-2);
        vector<real> X(num_grid_points-2);

        // Value of S
        S = S_req[i];

        // We are populating the first row and last manually while putting the boundary conditions

        // First wow
        Matrix[0][0] = -2.0;
        Matrix[0][0+1] = 2.0;
        B[0] = -(S*dr*dr) - (dtdr*2*r[1]*(1-(dr/(2*r[1]))));

        // Last Row
        Matrix[rows-1][rows-1] = -2;
        Matrix[rows-1][rows-2] = (1 - (dr/(2*r[rows-1])));
        B[rows-1] = -(S*dr*dr) - (Tn*(1 + (dr/(2*r[rows-1]))));

        // We are assigning the values  of non-zero elements of the matrix except for first and last row
        for (int row = 1;row < rows-1;row++)
        {
            // Only assign values to NON ZERO elements of the Matrix
            Matrix[row][row-1] = (1 - (dr/(2*r[row])));
            Matrix[row][row] = -2;
            Matrix[row][row+1] = (1 + (dr/(2*r[row])));
        }

        // We are assigning the values  of non-zero elements of B except for first and last row
        for (int row = 1;row<rows-1;row++)
        {
            B[row] = -S*dr*dr;
        }

        // Function Thomas_Algorithm_Matix_Solver(A,B,X) will solve the system of equation AX = B
        // where, A is a tridiagonal matrix
        Thomas_Algorithm_Matix_Solver(Matrix,B,X);
        //X[0] = X[2];
        //X[X.size()-1] = Tn;

        X.insert(X.begin(),X[1]);
        X.push_back(Tn);

        // Getting the address of the maximum element in the vector
        auto int_result = std::max_element(X.begin(),X.end());

        // Getting the value of the maximum element in the vector
        T_max.push_back(*int_result);
    }

    // Setting all the higher values of S to 0 where maximum temperature is greater than T_req.
    for(int ctr = 0;ctr<T_max.size();ctr++)
    {
        if(T_max[ctr] > T_req)
        {
            S_req[ctr] = 0;
        }
    }

    // Getting the address of the maximum element in the vector
    auto int_result = std::max_element(S_req.begin(),S_req.end());

    // Getting the value of the maximum element in the vector
    S_optimum.push_back(*int_result);

    write_to_file(S_optimum,"Q_5_c_Optimum_S_"+to_string(T_req)+".csv");

    Output_file_names.push_back("Q_5_c_Optimum_S_"+to_string(T_req)+".csv");
    return 0;
}
