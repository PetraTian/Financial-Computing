//
//  main.cpp
//  RepeatSquaring
//
//  Created by PC on 2018/11/5.
//  Copyright © 2018 PC. All rights reserved.
//

# include <iostream>
# include <random>
# include <cmath>
# include <fstream>
# include "newmat.h"
using namespace std;


// This function is used to print matrix
void print_matrix(Matrix A, int size)
{
    for (int i=1; i<=size; i++){
        for (int j=1; j<=size; j++)
            cout << A(i,j) << " ";
        cout << endl;
    }
    cout << endl;
}


// This function return a random number drawn from uniform distribution
double get_uniform()
{
    return float(random()) / (pow(2.0, 31.0)-1);
}


// This function use repeated squaring to calculate A^k
Matrix repeated_squaring(Matrix A, int exponent, int no_rows)
{
    if (exponent==0){
        IdentityMatrix I(no_rows);
        return I;
    }
    else if (exponent % 2 == 1)
        return A*repeated_squaring(A*A, (exponent-1)/2, no_rows);
    else
        return repeated_squaring(A*A, exponent/2, no_rows);
}


// This function use brute force method to calculate A^k
Matrix direct_multiplication(Matrix A, int exponent, int no_rows)
{
    Matrix result(no_rows, no_rows);
    result = A;
    for (int k=1; k<exponent; k++)
        result = result * A;
    return result;
}


// This function calculate the running time for RS and BF methods
// of exponent from 1 to 300 for 5*5 matrix
void time_for_different_exponent()
{
    int exponent, size=5;
    // fill matrix with random entries in interval (-5,5)
    Matrix A(size,size);
    for (int i=1; i<=size; i++){
        for (int j=1; j<=size; j++)
            A(i,j) = (get_uniform()-0.5)*10;
    }
    double time_RS[300];
    double time_BF[300];
    for (int e=1; e<=300; e++){
        exponent = e;
        
        // compute A^exponent using repeated squaring
        double time_before = clock();
        repeated_squaring(A, exponent, size);
        double time_after = clock();
        time_RS[e-1] = (float(time_after) - float(time_before)) / CLOCKS_PER_SEC;
        
        // compute A^exponent using brute force method
        time_before = clock();
        direct_multiplication(A, exponent, size);
        time_after = clock();
        time_BF[e-1] = (float(time_after) - float(time_before)) / CLOCKS_PER_SEC;
    }
    
    ofstream outfile("time_for_different_exponent.csv");
    for (int i=0; i<300; i++)
        outfile << time_RS[i] << "," << time_BF[i] << endl;
}

// This function calculate the running time for RS and BF methods
// of exponent 300 for matrix size from 1 to 300
void time_for_different_matrix_size()
{
    int exponent=300, size;
    
    double time_RS[300];
    double time_BF[300];
    for (int s=1; s<=300; s++){
        size = s;
        
        // fill matrix with random entries in interval (-5,5)
        Matrix A(size,size);
        for (int i=1; i<=size; i++){
            for (int j=1; j<=size; j++)
                A(i,j) = (get_uniform()-0.5)*10;
        }
        
        // compute A^exponent using repeated squaring
        double time_before = clock();
        repeated_squaring(A, exponent, size);
        double time_after = clock();
        time_RS[s-1] = (float(time_after) - float(time_before)) / CLOCKS_PER_SEC;
        
        // compute A^exponent using brute force method
        time_before = clock();
        direct_multiplication(A, exponent, size);
        time_after = clock();
        time_BF[s-1] = (float(time_after) - float(time_before)) / CLOCKS_PER_SEC;
    }
    
    ofstream outfile("time_for_different_matrix_size.csv");
    for (int i=0; i<300; i++)
        outfile << time_RS[i] << "," << time_BF[i] << endl;
}


int main(int argc, const char * argv[])
{
    // scan data from command line
    int exponent, size;
    sscanf(argv[1], "%d", &exponent);
    sscanf(argv[2], "%d", &size);
    cout << "The number of rows/columns in the square matrix is: " << size << endl;
    cout << "The exponent is: " << exponent << endl;
    
    // fill matrix with random entries in interval (-5,5)
    Matrix A(size,size);
    for (int i=1; i<=size; i++){
        for (int j=1; j<=size; j++)
            A(i,j) = (get_uniform()-0.5)*10;
    }
    
    // compute A^exponent using repeated squaring
    double time_before = clock();
    Matrix A_to_pow_k_RS = repeated_squaring(A, exponent, size);
    double time_after = clock();
    cout << "Repeatd squaring result: " << endl;
    double diff = (float(time_after) - float(time_before)) / CLOCKS_PER_SEC;
    cout << "It took " << diff << " seconds to complete." << endl;
    
    // compute A^exponent using brute force method
    time_before = clock();
    Matrix A_to_pow_k_BF = direct_multiplication(A, exponent, size);
    time_after = clock();
    cout << "Direct Multiplication result: " << endl;
    //print_matrix(A_to_pow_k_BF, size);
    diff = (float(time_after) - float(time_before)) / CLOCKS_PER_SEC;
    cout << "It took " << diff << " seconds to complete." << endl;
    
    time_for_different_exponent();
    
    time_for_different_matrix_size();
    return 0;
}




