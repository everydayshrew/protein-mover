#include "point_hh.h"
#include <fstream>

// statement?t:f

const int MAX_SIZE = 80;
// PLENGTH is the maximum length the protein can expand to.
const int PLENGTH = 6;
point box[MAX_SIZE][MAX_SIZE][MAX_SIZE];


// Xorshift RNG Algorithm by George Marsaglia.
uint32_t xor128(void) {
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;
 
  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}

// Function randMod returns -1, 0, or 1
// Function will not return a value that sets it out of bound of box.
int randMod(int a){
    int value;
    switch(xor128()%3){
       case 0: value = 0; break;
       case 1: value = 1; break;
       case 2: value = -1; break;
    }
    if ((value+a)>MAX_SIZE-2 || (value+a)<2){return 0;}else{return value;}
};



// Function modifies the field around a protein
void modSurroundingField(int x, int y, int z, int value){
     for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
           for(int k=-1;k<=1;k++){ 
              if(k!=0 || j!=0 || i!=0){
                  // Possible crash area?  Mods fields OUTSIDE the boundries of the array.
                  if(!((x+i>MAX_SIZE-2||x+i<2)||(y+j>MAX_SIZE-2||y+j<2)||(z+k>MAX_SIZE-2||z+k<2))){
                       box[x+i][y+j][z+k].modField(value);
                       
                       // Electric field around a protien is of 2 magnitude
                       box[x+i][y+j][z+k].modElectricField(box[x][y][z].sign*2*value);
                       if(!((x+2*i>MAX_SIZE-2||x+2*i<2)||(y+2*j>MAX_SIZE-2||y+2*j<2)||(z+2*k>MAX_SIZE-2||z+2*k<2))){
                          // Electric field around a protein 1 unit out is of 1 magnitude.
                          box[x+(2*i)][y+(2*j)][z+(2*k)].modElectricField(box[x][y][z].sign*value);                                                                        
                       }
                  }
                  //box[x+i][y+j][z+k].modField(value);
              }
           }
        }
     }
};

using namespace std;

class protein{
      public:
             // Constructor Class, Choose a random start for the protein.
             protein(){
                       srand(time(NULL));
                       xpath[0] = 8 + rand() % MAX_SIZE-4;
                       ypath[0] = 8 + rand() % MAX_SIZE-4;
                       zpath[0] = 8 + rand() % MAX_SIZE-4;
                       
                       x_start[0] = xpath[0];
                       y_start[0] = ypath[0];
                       z_start[0] = zpath[0];
                       counter = 0;
                       wins = 0;
                       t1 = 0; t2 = 0;  t3 = 0;

                       switch(xor128()%2){
                         case 0: box[xpath[0]][ypath[0]][zpath[0]].sign = -1;
                         case 1: box[xpath[0]][ypath[0]][zpath[0]].sign = 1;
                       };
                        
                       full = false;
                       box[xpath[0]][ypath[0]][zpath[0]].setStatus(true, 10);
                       
                       // Field setup
                       modSurroundingField(xpath[0],ypath[0],zpath[0], 1);

                       while(!full){createProtein();}
                       calculateCoM();
                       centerOfMassStart[0] = centerOfMass[0];
                       centerOfMassStart[1] = centerOfMass[1];
                       centerOfMassStart[2] = centerOfMass[2];
                       
                       cout << "Created protein of " << counter+1 << " length.\n";
                       //outfile.open("data.txt");
                       
             };
             
             bool isFull(){return full;};
             
             // Create the protein.
             void createProtein(){
                  int a=0,b=0,c=0,x=0;
                  // Check if new spot is occupied.  If not, set that as the new location
                  // continue checking spots until one is found.  Will give up after 20 tries.
                  while((box[xpath[counter]+a][ypath[counter]+b][zpath[counter]+c].checkStatus()) && x<=20
                          //&& (xpath[counter]+a < 2 || xpath[counter]+a > MAX_SIZE-2)
                          //&& (ypath[counter]+a < 2 || ypath[counter]+a > MAX_SIZE-2)
                          //&& (zpath[counter]+a < 2 || zpath[counter]+a > MAX_SIZE-2)
                          ){
                        a = 2*randMod(xpath[counter]);
                        b = 2*randMod(ypath[counter]);
                        c = 2*randMod(zpath[counter]);
                        x++;       
                        //cout << x;                                                                     
                  }
                  // Check if it successfully moved.
                  if(x>20){full = true;}
                  
                  // Incremenet counter if not full.
                  if(!full){
                     counter++;
                     // Write new path data
                     
                     xpath[counter] = xpath[counter-1] + a;
                     ypath[counter] = ypath[counter-1] + b;
                     zpath[counter] = zpath[counter-1] + c;
                     
                     //Set Sign
                     int sign = 0;
                     sign = xor128()%10;
                     if(sign<5){sign = -1;}else{sign = 1;}
                     
                     box[xpath[counter]][ypath[counter]][zpath[counter]].sign = sign;
                    
                     x_start[counter] = xpath[counter];
                     y_start[counter] = ypath[counter];
                     z_start[counter] = zpath[counter];
                     box[xpath[counter]][ypath[counter]][zpath[counter]].setStatus(true, 10);
                    
                     // Set Field
                     modSurroundingField(xpath[counter], ypath[counter], zpath[counter], 1);
                     
                     // If full, set the flag to full.
                     if(counter == PLENGTH-1){full = true;}
                     //cout << "COUNTER: " << counter << endl;
                  }
             };
             
             // Returns the distance between two nodes
             double nodeDistance(int a, int b){
                 return sqrt(pow(static_cast<double>(xpath[a]-xpath[b]),2)
                               +pow(static_cast<double>(ypath[a]-ypath[b]),2)
                               +pow(static_cast<double>(zpath[a]-zpath[b]),2));
             };
             
             // Returns the distance between a node and a spot
             double spotDistance(int a, int x, int y, int z){
                  return sqrt(pow(static_cast<double>(xpath[a]-x),2)
                               +pow(static_cast<double>(ypath[a]-y),2)
                               +pow(static_cast<double>(zpath[a]-z),2));  
             };
             
             // Will attempt to move each component of the protien ONCE
             // Currently moves down the chain. (Random choice in the futue?)
             // A protien has a value of 10.  The field around that
             // protein adds 1.  If the value exceeds 1, a protein can no longer
             // move to that location. 
             void moveProtein(){
                  //int x,y,z;
                  bool valid = true;
                  //for(int i=0; i<= counter; i++){
                          // KEEP AN EYE ON THIS.  counter doesn't seem to modify the last protien in the chain
                          // but PLENGTH is innacurate in the case that the protien does not reach full length
                          // (albiet this event has not happened yet through the amount of testing done)
                          int i = xor128()%PLENGTH;  
                  
                          // x, y, and z are the coordinats of the current node we wish to move
                          int x = xpath[i];
                          int y = ypath[i];
                          int z = zpath[i];
                          
                          // Uses the Probablitiy Algorithm to find a new spot:
                          ProbAlg(i, x, y, z);
                          
                          // Check 1: See if the most probable is valid.
                          if(x == xpath[i] && y == ypath[i] && z == zpath[i]){
                               //cout << "Failed (T1)";
                               valid = false;
                               t1++;
                          }
                  
                          //Move to a random position.
                          //x += randMod(x);
                          //y += randMod(y);
                          //z += randMod(z);
                          //cout << "Trying to move Atom " << i << " at [" << xpath[i] << ", " << ypath[i]
                          //     << ", " << zpath[i] << "] to [" << x << ", " << y << ", " << z
                          //     << "]     ";
                  
                          // Check 2: Check if it is a valid spot.
                          if(valid && box[x][y][z].returnField() > 1){
                              //cout << "Failed! (T2)";
                              valid = false;
                              t2++;
                          }
                                              
                          // OUTDATED 
                          // Check 2: If still valid, check if the new spot is within 2 range of another protein.
                          /*if(valid){
                              for(int j=0; j<=counter; j++){
                                  // Do not check self
                                  if(i!=j){
                                      if(spotDistance(j,x,y,z)<2){
                                         valid = false;                        
                                      }
                                  }   
                              }
                          }     */
                          
                          // Check 3: If still valid, check if the adjacent protien nodes are within
                          // the sqrt(10) distance.
                          if(valid){
                              //double ALPHA, BETA;
                              //if(i!=PLENGTH-1){ALPHA = spotDistance(i+1,x,y,z);}else{ALPHA = 0;}
                              //if(i!=0){BETA = spotDistance(i-1,x,y,z);}else{BETA = 0;}
                              //cout << fixed << setprecision(2);
                              //cout << "\n     Distance to node... [" << i-1 << "]: " << BETA << "   [" << i+1 << "]: " << ALPHA << "  ";
                              if(i==0){
                                  if(spotDistance(i+1,x,y,z) > sqrt(10)){valid = false; t3++;}//cout<<"Failed! (T3)";}
                              } else if(i==counter){
                                  if(spotDistance(i-1,x,y,z) > sqrt(10)){valid = false; t3++;}//cout<<"Failed! (T3)";}
                              } else {
                                  if(spotDistance(i+1,x,y,z) > sqrt(10) ||
                                     spotDistance(i-1,x,y,z) > sqrt(10)){valid = false; t3++;}//cout<<"Failed! (T3)";}
                              }
                          }
                          
                          // Result: If still valid after checking all flags, move the protien
                          // to that seleted spot.
                          if(valid){
                              //cout << "Succeeded! ";
                              // First copy the sign.
                              int temp = box[xpath[i]][ypath[i]][zpath[i]].sign;
                              // Mod the field around the existing particle correctly.
                              modSurroundingField(xpath[i], ypath[i], zpath[i], -1);
                              // Remove the proten.
                              box[xpath[i]][ypath[i]][zpath[i]].setStatus(false, 1);      
                                    
                              // Change index path
                              xpath[i] = x;
                              ypath[i] = y;
                              zpath[i] = z;
                              
                              // Set new protein.
                              box[x][y][z].setStatus(true, 10, temp);
                              modSurroundingField(x,y,z, 1);
                              wins++;
                              writePath(static_cast<int>(wins));
                          }              
                          
                          //Reset flag
                          valid = true;
                          //cout << endl;   
                  //}     
             };
             
             void printPath(){
                  for(int i=0; i<=counter; i++){
                       cout << "(" << xpath[i] << ", " << ypath[i] << ", " << zpath[i] << ")[";
                       if(box[xpath[i]][ypath[i]][zpath[i]].sign > 0){cout << "+]   ";}
                       if(box[xpath[i]][ypath[i]][zpath[i]].sign < 0){cout << "-]   ";}
                       if(box[xpath[i]][ypath[i]][zpath[i]].sign == 0){cout << "0]   ";}
                       if((i+1)%3==0 && i!=0){cout << endl;}
                  }   
             };
             
             void writePath(int a){
                  outfile.open("data.txt", std::ios_base::app);
                  outfile << a << "  ";
                  for(int i=0; i<=counter; i++){
                       outfile << "(" << xpath[i] << "," << ypath[i] << "," << zpath[i] << ")   ";
                  }
                  outfile << endl;
                  outfile.close();
             };   
             
             // In an electric field setup of integers representing the strength, for positive proteins
             // it will seek the lowest value, and negitive proteins will seek the highest (which is
             // obtained by inverting the method of searching for lowest).  Using a porportionality
             // function, the probability of the lowest (most likely) is some value, p1.  And a value
             // 1 off from p1 is related by the function 1/c * p1.
             // Thus p1 = 1 / sum(an/c^(n-1))   Where an is the 'amount' of some value present in the field.
             // or pn = [1 / c^(n-1)] * p1
             void ProbAlg(int index, int& x, int& y, int& z){
                  // Random Probability
                  double P = (xor128()%10000000)/10000000.0;
                  x = xpath[index];
                  y = ypath[index];
                  z = zpath[index];
                  
                  int p,q,r;
                  int lowest = 2 * box[x][y][z].sign; // Lowest value to find p1.
                  double c = 26; // Value of c.  True value is unknow and must be figured out.
                  int a[27]; // an
                  double p1 = 0; // p1.
                  int n = 0; //index of a.
                  int sum = 0; // Sum of array a.
                  int a_size = 0; // Size of the array a.
                  
                  for(int i=0; i<=26; i++){a[i] = 0;} // Quick fill of array a.
                  
                  
                  
                  for(int i=-1;i<=1;i++){
                    for(int j=-1;j<=1;j++){
                       for(int k=-1;k<=1;k++){ 
                          if(k!=0 || j!=0 || i!=0){
                              //cout << "[" << x+i << "][" << y+j << "][" << z+k << "]:  " << box[x+i][y+j][z+k].Evalue << "\n";
                              // Possible crash area?  Looks at fields OUTSIDE the boundries of the array.
                              if(!((x+i>MAX_SIZE||x+i<0)||(y+j>MAX_SIZE||y+j<0)||(z+k>MAX_SIZE||z+k<0))){
                                  // Due to the problem with 0, the rules change according to the value of the
                                  // sign changes the if statement below.                                                                                                         
                                   //cout << "  " << box[x+i][y+j][z+k].Evalue << "*" << box[x][y][z].sign << " < " << lowest << "  ";
                                   if(box[xpath[index]][ypath[index]][zpath[index]].sign > 0){
                                       if((box[x+i][y+j][z+k].Evalue < lowest)){
                                           lowest = box[x+i][y+j][z+k].Evalue;   
                                       }                                                                                                
                                   } else {
                                       if((box[x+i][y+j][z+k].Evalue > lowest)){
                                           lowest = box[x+i][y+j][z+k].Evalue;   
                                       }  
                                   }
                                   
                              }
                          }
                       }
                    }
                 }   
                 //system("PAUSE");
                 //cout << "\nSign is: " << box[x][y][z].sign;
                 //cout << "\nLowest is: " << lowest;
                 
                 // Hey this is totally terrible programming!  Calculates the entireity of an.
                 while(sum < 26){
                       for(int i=-1;i<=1;i++){
                        for(int j=-1;j<=1;j++){
                           for(int k=-1;k<=1;k++){ 
                              if(k!=0 || j!=0 || i!=0){
                                 if(box[x+i][y+j][z+k].Evalue == lowest + (box[x][y][z].sign * n)){
                                    a[n] += 1;
                                    sum++;                                                                                     
                                 }
                              }
                           }
                        }
                     }        
                     n++;
                 }   
                 
                 a_size = n;
                 
                 // Calculate p1.
                 p1 = a[0];
                 for(int i = 1; i < a_size; i++){
                     p1 += a[i]/(pow(c,i-1));
                 }
                 p1 = (1/p1);
                 
                 //cout << "\n\np1 is: " << p1;
                 //cout << "\nP  is: " << P;
                 //system("PAUSE");
                 
                 // Now find the spot.
                 int count = 0;
                 bool found = false;
                 // Percentage values.
                 double first = 0;
                 double last = 0;
                 n = 0;
                 while(!found){
                     if (n==0){
                        first = 0;
                        last = a[n]*p1;
                        //cout << "Attempt 0: 0 < " << P << " < " << last << endl;
                        if(first<P && P<last){    
                           while(!found){
                              x += randMod(x);
                              y += randMod(y);
                              z += randMod(z);          
                              
                              if((box[x][y][z].Evalue == lowest) && 
                                 !(x==xpath[index] && y==ypath[index] && z==zpath[index])){
                                 found = true;
                                 //cout << "\nOld: " << xpath[index] << ", " << ypath[index] << ", " << zpath[index];
                                 //cout << "\nNew: " << x << ", " << y << ", " << z;
                              }
                              count++;
                              // Failsafe if it can't find said value.
                              if (count >= 500 && !found){
                                 //cout << "\nThe 80 checks failed!";
                                 found = true;
                              }  
                              
                              // Reset values
                              if(!found){
                                 x = xpath[index];
                                 y = ypath[index];
                                 z = zpath[index];
                              }                
                           }
                        }
                     } else {
                            
                        // Find the next segment of probability.
                        first = last;
                        if (n+1 > a_size){
                           last = 1;
                        } else {
                           if(a[n-1]==0){
                              last = first;
                           }else{
                              last += (a[n]*p1)/(pow(c,n));
                           }
                        }
                        //cout << "Attempt " << n-1 << ": " << first << " < " << P << " < " << last << endl;
                        if(first<P && P<last){
                           bool flag = false;
                           while(!found){
                              x += randMod(x);
                              y += randMod(y);
                              z += randMod(z);          
                              if((box[x][y][z].Evalue == lowest + (box[xpath[index]][ypath[index]][zpath[index]].sign * (n-1))) 
                                                      && !(x==xpath[index] && y==ypath[index] && z==zpath[index])){
                                 found = true;
                              }
                              count++;
                              //cout << box[x][y][z].Evalue << " = " << lowest << " + (" 
                              //     << box[xpath[index]][ypath[index]][zpath[index]].sign << " * " << n-1 << ")\n";
                              // Failsafe if it can't find said value.
                              if (count >= 500 && !found){
                                 found = true;
                                 x = xpath[index];
                                 y = ypath[index];
                                 z = zpath[index];
                                 // Below is debug code.  It should stay commented.
                                 /*cout << a[0]*p1 << "\n";
                                 cout << "(" << x << "," << y << "," << z << ")[" << box[xpath[index]][ypath[index]][zpath[index]].sign << "]\n";
                                 for(int i = 0; i<a_size; i++){
                                     cout << "\na[" << i << "]: " << a[i] << "  (" << lowest + (box[xpath[index]][ypath[index]][zpath[index]].sign * i)
                                      << ")";
                                 }
                                 cout << endl;
                                 cout << "\n\nI was looking for " 
                                      << lowest + (box[xpath[index]][ypath[index]][zpath[index]].sign * (n-1))
                                      << " and couldn't find it!\n";
                                 cout << "My lowest was " << lowest << " and my n was " << n-1 << endl;
                                 for(int i=-1;i<=1;i++){
                                    for(int j=-1;j<=1;j++){
                                       for(int k=-1;k<=1;k++){ 
                                          if(k!=0 || j!=0 || i!=0){
                                          cout << "[" << x+i << "][" << y+j << "][" << z+k << "]:  " << box[x+i][y+j][z+k].Evalue;
                                          cout << " (" << box[x+i][y+j][z+k].value << ") " << "\n";
                                          }
                                       }
                                    }
                                 }  
                                 printPath();
                                 sum = 0; n = 0;
                                  for(int i=0; i<=26; i++){a[i] = 0;}
                                 while(sum < 26 && n < 26){
                                       cout << "\nN: " << n << "\n----\n";
                                       for(int i=-1;i<=1;i++){
                                        for(int j=-1;j<=1;j++){
                                           for(int k=-1;k<=1;k++){ 
                                              if(k!=0 || j!=0 || i!=0){
                                                 cout << box[x+i][y+j][z+k].Evalue << " == " << lowest
                                                 << " +(" << box[x][y][z].sign << " * " << n << ")" << endl;
                                                 if(box[x+i][y+j][z+k].Evalue == lowest + (box[x][y][z].sign * n)){
                                                    a[n] += 1;
                                                    cout << a[n] << " ";
                                                    sum++;                                                                                     
                                                 }
                                              }
                                           }
                                        }
                                     }        
                                     n++;
                                     cout << endl;
                                 }   
                                 system("PAUSE");*/
                              }  
                              
                              // Reset values
                              if(!found){
                                x = xpath[index];
                                y = ypath[index];
                                z = zpath[index];
                              }
                           }
                        }    
                     }
                     n++;
                     //cout << "N is now:   " << n << endl;
                 }   
             };
             
             void printDisplacement(){
                  calculateCoM();
                  
                  for(int i=0; i<=counter; i++){
                       cout << i+1 << ": " << "(" << x_start[i] << ", " << y_start[i] << ", " << z_start[i] << ")" << " -> "
                            << "(" << xpath[i] << ", " << ypath[i] << ", " << zpath[i] << ")" << "       Displacemnet: "
                            << "(" << xpath[i]-x_start[i] << ", " << ypath[i]-y_start[i] << ", " << zpath[i]-z_start[i] 
                            << ")\n";
                  }
                  cout << "Center of Mass:    (" << centerOfMass[0] << ", " << centerOfMass[1] << ", " << centerOfMass[2]
                       << ") -> (" << centerOfMassStart[0] << ", " << centerOfMassStart[1] << ", "
                       << centerOfMassStart[2] << ")     \n                           Displacement: ("
                       << centerOfMass[0] - centerOfMassStart[0] << ", " << centerOfMass[1] - centerOfMassStart[1]
                       << ", " << centerOfMass[2] - centerOfMassStart[2] << ")\n";
                  cout << "Avg Velocity:\n" << fixed << setprecision(5) << "        (" 
                       << (centerOfMass[0]-centerOfMassStart[0])/wins << ", "
                       << (centerOfMass[1]-centerOfMassStart[1])/wins << ", "
                       << (centerOfMass[2]-centerOfMassStart[2])/wins << ") per tick\n";
                  cout << "Number of successful moves: " << wins << "\n";
                  cout << fixed << setprecision(2);
                  cout << "Failures:   [T1] " << t1 << "    [T2] " << t2/1000000 << "x10^6    [T3] " << t3/1000000 << "x10^6 \n";
                  //outfile.close();
             };
             
             void calculateCoM(){
                  double x=0,y=0,z=0;
                  for(int i=0; i<=counter; i++){
                       x += xpath[i]; 
                       y += ypath[i];
                       z += zpath[i];
                  }
                  
                  centerOfMass[0] = x/(counter+1);
                  centerOfMass[1] = y/(counter+1);
                  centerOfMass[2] = z/(counter+1);
                   
             };
             
             
             // Center of mass variables
             double centerOfMass[3];
             double centerOfMassStart[3];
             
             // Stores data of chosen path of protein
             int xpath[15];
             int ypath[15];
             int zpath[15];
             
             // Copy of the initial position of the protiens, used for overal displacement
             int x_start[15];
             int y_start[15];
             int z_start[15];
             int counter;
             bool full;
             
             // Amount of times it's successfully moved a protein.
             double long wins;
             
             // Amount of times it failed.
             double long t1, t2, t3;
             
             // Output file
             ofstream outfile;
             
};
