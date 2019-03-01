# ifndef K_HANOI_H
# define K_HANOI_H

# include <iostream>
using namespace std;

# include <deque>
using std::deque;

# include <cmath>


class Towers{
	int num_step;
	int num_disk;
	int num_peg;
	deque < deque <int> > towers;

	void initialize(int disks, int pegs){
		num_step=0;
		num_disk=disks;
		num_peg=pegs;

		// construct tower
		for (int i=0; i<num_peg; i++){
			deque <int> peg;
			towers.push_back(peg);
		}
		for (int i=0; i<num_disk; i++){
			towers[0].push_back(i+1);
		}
	}

	// See whether smaller disks are placed below larger disks
	// if yes then return false, else return true
	bool isLegal(){
		for (int i=0; i<num_peg; i++){
			for (int j=0; j<towers[i].size(); j++){
				for (int k=j; k<towers[i].size(); k++){
					if (towers[i][k] < towers[i][j])
						return false;
				}
			}
		}
		return true;
	}

	// calculate what p to chose
	int num_disk_to_move(int num_disk_left, int num_peg_left){
		if (1 == num_peg_left)
			return (num_disk_left-1);
		else
			return (num_disk_left/2);
	}

	void move_top_disk(int source, int dest){
		num_step++;
		cout << "move disk " << towers[source].front() << " from " << source+1 << " to " << dest+1;
		towers[dest].push_front(towers[source].front());
		towers[source].pop_front();
		if (isLegal()==true){
			cout << " (Legal)" << endl;
		}
		else{
			cout << " (Ilegal)" << endl;
		}
	}

	void move(int n, int source, int dest, deque <int> intermd){
		if (intermd.size()>0){
			if (n>1){
				int p=num_disk_to_move(n, intermd.size());

				deque <int> temp=intermd;
				temp.pop_front();
				deque <int> intermd2=temp;
				temp.push_back(dest);
				deque <int> intermd1=temp;
				temp.pop_back();
				temp.push_front(source);
				deque <int> intermd3=temp;

				move(p, source, intermd.front(), intermd1);
				move(n-p, source, dest, intermd2);
				move(p, intermd.front(), dest, intermd3);
			}
			else
				move_top_disk(source, dest);
		}
		else if(0==intermd.size())
			move_top_disk(source, dest);
	}


	void print_peg_states(int n){
		cout << "----------------------------------------" << endl;
		cout << "The state of the " << n+1 << " peg: " << endl;
		for (int j=0; j<towers[n].size(); j++){
			cout << towers[n][j] << " " << endl;
		}
		cout << "Number of step: " << num_step << endl;
		cout << "----------------------------------------" << endl;
	}

public:
	void solve(int num_peg, int num_disk){
		initialize(num_disk, num_peg);
		print_peg_states(0);
		deque <int> intermd; // intermediate pegs
		for (int i=1; i<num_peg-1; i++)
			intermd.push_back(i);
		move(num_disk, 0, num_peg-1, intermd);
		print_peg_states(num_peg-1);
	}
};

# endif