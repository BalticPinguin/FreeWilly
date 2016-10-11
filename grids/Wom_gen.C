#include "Wom.h"

// see http://web.maths.unsw.edu.au/~rsw/Sphere/
int Wom_precision_table( unsigned int rule )
{
   unsigned int rule_max=10;
   int table[10]= {4, 9, 16, 25, 36, 49, 64, 81, 100, 121};

   if ( rule < 1 )
   {
      std::cerr << "\n";
      std::cerr << "PRECISION_TABLE - Fatal error!\n";
      std::cerr << "  RULE < 1.\n";
      exit ( 1 );
   }
   else if ( rule_max < rule )
   {
      std::cerr << "\n";
      std::cerr << "PRECISION_TABLE - Fatal error!\n";
      std::cerr << "  RULE_MAX < RULE.\n";
      exit ( 1 );
   }

   return table[rule-1];
}

void Wom_points (int rule, unsigned int order, double* x, double* y, double* z, double* w)
{
   switch (rule)
   {
   case 1:{  // ev
      switch (order){
      case 1:
         ev01_0004 (x,y,z,w);
      case 2:
         ev02_0009 (x,y,z,w);
      case 3:
         ev03_0016 (x,y,z,w);
      case 4:
         ev04_0025 (x,y,z,w);
      case 5:
         ev05_0036 (x,y,z,w);
      case 6:
         ev06_0049 (x,y,z,w);
      case 7:
         ev07_0064 (x,y,z,w);
      case 8:
         ev08_0081 (x,y,z,w);
      case 9:
         ev09_0100 (x,y,z,w);
      case 10:
         ev10_0121 (x,y,z,w);
      default:{
         std::cerr << "\n";
         std::cerr << "WOB_POINTS - Fatal error!\n";
         std::cerr << "- Order unknown.\n";
         exit ( 1 );
      }}
   }
   case 2:{ //md
      switch (order){
      case 1:
         md01_004 (x,y,z,w);
      case 2:
         md02_009 (x,y,z,w);
      case 3:
         md03_016 (x,y,z,w);
      case 4:
         md04_025 (x,y,z,w);
      case 5:
         md05_036 (x,y,z,w);
      case 6:
         md06_049 (x,y,z,w);
      case 7:
         md07_064 (x,y,z,w);
      case 8:
         md08_081 (x,y,z,w);
      case 9:
         md09_100 (x,y,z,w);
      case 10:
         md10_121 (x,y,z,w);
      default:{
         std::cerr << "\n";
         std::cerr << "WOB_POINTS - Fatal error!\n";
         std::cerr << "- Order unknown.\n";
         exit ( 1 );
      }}
   }
   case 3:{ //mn
      switch (order){
      case 1:
         md01_004 (x,y,z,w);
      case 2 :
         md02_009 (x,y,z,w);
      case 3 :
         md03_016 (x,y,z,w);
      case 4:
         md04_025 (x,y,z,w);
      case 5:
         md05_036 (x,y,z,w);
      case 6:
         md06_049 (x,y,z,w);
      case 7:
         md07_064 (x,y,z,w);
      case 8:
         md08_081 (x,y,z,w);
      case 9:
         md09_100 (x,y,z,w);
      case 10:
         md10_121 (x,y,z,w);
      default:{
         std::cerr << "\n";
         std::cerr << "WOB_POINTS - Fatal error!\n";
         std::cerr << "- Order unknown.\n";
         exit ( 1 );
      }}
   }
   case 4:{ //me
      switch (order){
      case 1:
         me01_004 (x,y,z,w);
      case 2 :
         me02_009 (x,y,z,w);
      case 3 :
         me03_016 (x,y,z,w);
      case 4:
         me04_025 (x,y,z,w);
      case 5:
         me05_036 (x,y,z,w);
      case 6:
         me06_049 (x,y,z,w);
      case 7:
         me07_064 (x,y,z,w);
      case 8:
         me08_081 (x,y,z,w);
      case 9:
         me09_100 (x,y,z,w);
      case 10:
         me10_121 (x,y,z,w);
      default:{
         std::cerr << "\n";
         std::cerr << "WOB_POINTS - Fatal error!\n";
         std::cerr << "- Order unknown.\n";
         exit ( 1 );
      }}
   }
   default:{
         std::cerr << "\n";
         std::cerr << "WOB_POINTS - Fatal error!\n";
         std::cerr << "- Rule unknown.\n";
         exit ( 1 );
   }
   }
}

bool unavailable(int rule){
   if(rule<1)
      return true;
   if(rule>10)
      return true;
   return false;
}

int max_avail(){
   return 10;
}
