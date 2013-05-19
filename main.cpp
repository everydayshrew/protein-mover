#include <iomanip>
#include "protein_hh.h"

// Get animation to work.
// COMPLETED:  Occupy the area (the 3x3 node surrounding the particle) with a 'can move here' flag (ints)
// COMPLETED: Translate the Protein (can only go to 'can move here' flag distances between 2 and sqrt(10)
// NOTE: sqrt(8) should be AVOIDED.


/*
NOTES: The protein move program currently, with the 0,1,10 flags will excute ~30 times for every
move the protien does to update the flags around it.

The previous method, which checked every protein with the one caused PLENGTH excutions for every
POSSIBLE move of the protein.
*/

void print(){
      int size = 39, counter = 1;
      
      for(int i=0;i<=size;i++){
       for(int j=0;j<=size;j++){
          for(int k=0;k<=size;k++){ 
             if(box[i][j][k].returnField() > 0){
                cout << "(" << i << ", " << j << ", " << k << "): " << box[i][j][k].returnField();
                if(counter%4==0){cout << endl;}else{cout << "  ";}
                counter++;
                }
             }
          }
       }
       cout << endl;
       system("PAUSE");
}

// The art of molecular Dynamics Simulation -Rapport

int main()
{   
    int a, b, c;
    protein A;
    //for(int i=0; i < 30; i++){cout << randMod(rand()%100) << " ";}
    A.printPath();
    cout << endl;
    //cout << endl;
    system("PAUSE");
    
    int trials = 1000000;
    cout << "\n\nNow Calculating: 0%[--------------------]100%\n                    ";
    for(int i=0; i<=trials; i++){
            A.moveProtein();
            if(i%(trials/20)==0){cout << "X";}
    }
    cout << "" << endl << endl;
    A.printDisplacement();
       
       //cin >> a >> b >> c;
       //cout << "State: " << box[a][b][c].returnField();
    
    system("PAUSE");
    return 1;
}
