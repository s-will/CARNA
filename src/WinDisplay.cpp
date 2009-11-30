#include <CImg.h>
#include <iostream>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>

using namespace cimg_library;
using namespace std;
using namespace Gecode;

class WinDisplay {

  CImg<unsigned char> imageOut; 
  CImgDisplay main_disp;

  int nvariables;
  int domain_size;
  int scalex;
  int scaley;
  int img_x;
  int img_y;  

public: 
  WinDisplay(int var, int dom){   
    int min_win_size=200;

    nvariables=var;
    domain_size=dom;
    scalex=1;
    scaley=1;
    img_x=0;
    img_y=0;    
    if (var<min_win_size){
      scalex=min_win_size/var;
    }
    if (dom<min_win_size){
      scaley=min_win_size/dom;
    }
    img_x=var*scalex;
    img_y=dom*scaley;

    imageOut = CImg<>(img_x,img_y,1,3,0);
    cout << "allocating image" << img_x << " " << img_y << "\n";
    
    for (int i=0;i<img_x;i++)
      for(int j=0;j<img_y;j++)
	imageOut(i,j,1)=200;   //uninitialized
        
    main_disp=CImgDisplay();
   }

  ~WinDisplay(){
    cout << "calling destructor\n";
    //delete imageOut;
    //delete main_disp;
  }

  void display(){
    cout << "display\n";
    if (main_disp.is_empty() || main_disp.is_closed())
      main_disp=CImgDisplay(imageOut,"Bitmap",0);
    main_disp.display(imageOut);  
  }

  void update(IntVarArray q){
    cout << "update\n";    

    // clean image
    for (int i=0;i<img_x;i++){
      for(int j=0;j<img_y;j++){
	imageOut(i,j,1)=0;
      }
    }

    // for each domain value paint block if it is in the current domain
    for (int i=0;i<q.size();i++){
      for (int j=0;j<domain_size;j++){

	if (q[i].in(j)){
	  for (int dx=0;dx<scalex;dx++)
	    for (int dy=0;dy<scaley;dy++)
	      imageOut(i*scalex+dx,j*scaley+dy,1)=200;
	}
      }

    }

    cout <<"ok update\n";
    display();    
  }

  void close(){
    cout << "close (wait 3 sec)\n";
    if (!main_disp.is_empty() && !main_disp.is_closed()) {
      for (int i=0;i<3;i++){
	cout << 3 - i << endl;
	main_disp.wait(1000);
      }
      main_disp.close();
    }
  }


 }; 
