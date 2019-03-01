#include <iostream>
#include <ctime>
#include "card_game.h"

using namespace std;
	
int main (int argc, char * const argv[]) 
{
	int total_number_of_cards;
	sscanf (argv[1], "%d", &total_number_of_cards);
	int size = total_number_of_cards/2+1;

	double **cache = new double*[size];
	for (int i=0; i<size; i++){
		cache[i] = new double[size];
	}
	for (int i=0; i!=size; i++){
		for (int j=0; j!=size; j++)
			cache[i][j] = -1;
	}

	cout << "Total Number of Cards = " << total_number_of_cards << endl;
	cout << "Value of the game = " << value(size-1,size-1,cache) << endl;

    return 0;
}
