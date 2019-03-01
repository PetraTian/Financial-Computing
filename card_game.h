/*
 *  card_game.h
 *  Card Game
 *
 *  Created by Ramavarapu Sreenivas  
*/

#ifndef	CARDGAME_H
#define CARDGAME_H
#include <algorithm>


double value(int r, int b, double **cache)
{
	if (0 == r)
		return ((double) b);
	if (0 == b)
		return (0);
	
	double term1, term2;
	if (cache[r-1][b] == -1){
		double temp_num = value(r-1, b, cache);
		term1 = ((double) r/(r+b)) * temp_num;
		cache[r-1][b] = temp_num;
	}
	else
		term1 = ((double) r/(r+b) * cache[r-1][b]);
	
	if (cache[r][b-1] == -1){
		double temp_num = value(r, b-1, cache);
		term2 = ((double) b/(r+b)) * temp_num;
		cache[r][b-1] = temp_num;
	}
	else
		double term2 = ((double) r/(r+b) * cache[r][b-1]);
	
	return std::max((term1 + term2), (double) (b - r));
}
#endif