#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void print_vector_s(vector<string> const &ptr) //function to print a vector of strings
{
 for(unsigned int i=0; i<ptr.size(); i++)
 {
  cout << ptr.at(i) << endl;
 }
 cout << endl;
}

void Read_spew()
{
 ifstream file;
 file.open("45_0210_1.txt");
 
 char output[500];
 vector<string> output_list;

 if (file.is_open())
 {
  while(!file.eof())
  {
   file >> output; 
   output_list.push_back(output);
  }
 }
 
 file.close();
//print_vector_s(output_list);

int i=0; //counting variable for how many photons are produced in killMedia
/*for(unsigned int j=0; j<output_list.size(); j++)
 {
  if(output_list.at(j)=="E_PbWO4" && output_list.at(j-4)=="1" && output_list.at(j+21)=="1"
  //Number of photons produced exactly like in the SDAction
  )
  {
   i+=1;
   //cout << output_list.at(j-4) << endl;
  }
 }*/

/*for(unsigned int j=0; j<output_list.size(); j++)
{
 if (output_list.at(j)=="Air" && output_list.at(j+61)=="1" //parent id of next step
 && output_list.at(j+65)=="E_PbWO4" && output_list.at(j+86)=="1"//just produced
 )
 {
  i+=1;
 }
}*/

for(unsigned int j=0; j<output_list.size(); j++)
{
 if (output_list.at(j)=="GG475_filter" //&& output_list.at(j-14)=="166" //only left kill media
 //&& output_list.at(j-4)=="1"//Parent ID is not 1
 )
 {
  cout << output_list.at(j-17)<< endl;
  //cout << output_list.at(j-19)<< endl;
  i+=1;
 }
}

/*for(unsigned int j=0; j<output_list.size(); j++)
{
 if (output_list.at(j)=="escaped" && output_list.at(j+1)=="scintillator" //only left kill media
 && output_list.at(j+2)=="photon"
 )
 {
  i+=1;
  cout << output_list.at(j-86)<<" "<< output_list.at(j-81) << endl; //Looking at pre and post position in crystal PbWO4
 }
}*/

 cout << i << endl;
}
