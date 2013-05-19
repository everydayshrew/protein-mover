#include <cstdlib>
#include <iostream>
#include <cmath>

class point{
      public:
             // Construtor
             point(){
                     x_force = 0; y_force = 0; z_force = 0; occupied = false; value = 0;; Evalue = 0; sign = 0;      
             };
             
             // Check if Point is Occupied
             bool checkStatus(){return occupied;};
             
             // Functions for field
             int returnField(){return value;};
             void modField(int x){value += x;};
             
             // Functions for Electric Field
             int returnElectricField(){return Evalue;};
             void modElectricField(int x){Evalue += x;};
             
             // Occupy Point
             void setStatus(bool a, int b){occupied = a; value = b;};
             void setStatus(bool a, int b, int c){occupied = a; value = b; sign = c;};
             void modStatus(bool a, int b){occupied = a; value += b;};
             
             // Values of Field Force at Point
             double x_force;
             double y_force;
             double z_force;
             
             // Occupy Status
             bool occupied;
             
             // The value indicates protiens and the fields around them for moving
             // proceudres.  See the moveProtein() function in the protein header
             // for more information of this.
             int value;       
             
             // This value is the electric field values used for probabilty crunching
             // See the ProbAlg() function in the protein header for more information.
             int Evalue;
             
             // Sign of that destinct protien at that point
             int sign;
};
             
             
