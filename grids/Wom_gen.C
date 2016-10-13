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

void Wom_points (int rule, unsigned int num_pts, double* x, double* y, double* z, double* w)
{
   switch (rule)
   {
   case 1:{  // ev
      switch (num_pts){
      case 4:
         ev01_0004 (x,y,z,w);
         return;
      case 9:
         ev02_0009 (x,y,z,w);
         return;
      case 16:
         ev03_0016 (x,y,z,w);
         return;
      case 25:
         ev04_0025 (x,y,z,w);
         return;
      case 36:
         ev05_0036 (x,y,z,w);
         return;
      case 49:
         ev06_0049 (x,y,z,w);
         return;
      case 64:
         ev07_0064 (x,y,z,w);
         return;
      case 81:
         ev08_0081 (x,y,z,w);
         return;
      case 100:
         ev09_0100 (x,y,z,w);
         return;
      case 121:
         ev10_0121 (x,y,z,w);
         return;
      default:{
         std::cerr << "\n";
         std::cerr << "WOB_POINTS - Fatal error!\n";
         std::cerr << "- Order unknown.\n";
         exit ( 1 );
      }}
   }
   case 2:{ //md
      switch (num_pts){
      case 4:
         md01_004 (x,y,z,w);
         return;
      case 9:
         md02_009 (x,y,z,w);
         return;
      case 16:
         md03_016 (x,y,z,w);
         return;
      case 25:
         md04_025 (x,y,z,w);
         return;
      case 36:
         md05_036 (x,y,z,w);
         return;
      case 49:
         md06_049 (x,y,z,w);
         return;
      case 64:
         md07_064 (x,y,z,w);
         return;
      case 81:
         md08_081 (x,y,z,w);
         return;
      case 100:
         md09_100 (x,y,z,w);
         return;
      case 121:
         md10_121 (x,y,z,w);
         return;
      default:{
         std::cerr << "\n";
         std::cerr << "WOB_POINTS - Fatal error!\n";
         std::cerr << "- Order unknown.\n";
         exit ( 1 );
      }}
   }
   case 3:{ //mn
      switch (num_pts){
      case 4:
         md01_004 (x,y,z,w);
         return;
      case 9 :
         md02_009 (x,y,z,w);
         return;
      case 16 :
         md03_016 (x,y,z,w);
         return;
      case 25:
         md04_025 (x,y,z,w);
         return;
      case 36:
         md05_036 (x,y,z,w);
         return;
      case 49:
         md06_049 (x,y,z,w);
         return;
      case 64:
         md07_064 (x,y,z,w);
         return;
      case 81:
         md08_081 (x,y,z,w);
         return;
      case 100:
         md09_100 (x,y,z,w);
         return;
      case 121:
         md10_121 (x,y,z,w);
         return;
      default:{
         std::cerr << "\n";
         std::cerr << "WOB_POINTS - Fatal error!\n";
         std::cerr << "- Order unknown.\n";
         exit ( 1 );
      }}
   }
   case 4:{ //me
      switch (num_pts){
      case 4:
         me01_004 (x,y,z,w);
         return;
      case 9 :
         me02_009 (x,y,z,w);
         return;
      case 16 :
         me03_016 (x,y,z,w);
         return;
      case 25:
         me04_025 (x,y,z,w);
         return;
      case 36:
         me05_036 (x,y,z,w);
         return;
      case 49:
         me06_049 (x,y,z,w);
         return;
      case 64:
         me07_064 (x,y,z,w);
         return;
      case 81:
         me08_081 (x,y,z,w);
         return;
      case 100:
         me09_100 (x,y,z,w);
         return;
      case 121:
         me10_121 (x,y,z,w);
         return;
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
