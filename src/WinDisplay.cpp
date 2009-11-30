#include <CImg.h>
#include <iostream>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include <gecode/set.hh>

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

  char* title;

public: 
  WinDisplay(int var, int dom){   
    WinDisplay(var,dom,"Untitled");
  }

  WinDisplay(int var, int dom, char* t){   
    int min_win_size=200;

    title=t;
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
    cout << "allocating image: " << img_x << " " << img_y << "\n";
    
    for (int i=0;i<img_x;i++)
      for(int j=0;j<img_y;j++)
	imageOut(i,j,1)=200;   //uninitialized
        
        main_disp=CImgDisplay();
	//main_disp.close();
   }

  ~WinDisplay(){
    //    cout << "calling destructor\n";
    //delete imageOut;
    //delete main_disp;
  }

  void display(){

    if (main_disp.is_empty() || main_disp.is_closed())
      main_disp=CImgDisplay(imageOut,title,0);
    main_disp.display(imageOut);  
  }

  void update(IntVarArray q){

    // clean image
    for (int i=0;i<img_x;i++){
      for(int j=0;j<img_y;j++){
	rgb(i,j,0,0,0);
      }
    }

    // for each domain value paint green block if it is in the current domain
    for (int i=0;i<q.size();i++){
      for (int j=0;j<domain_size;j++){

	if (q[i].in(j)){
	  for (int dx=0;dx<scalex;dx++)
	    for (int dy=0;dy<scaley;dy++)

	      rgb(j*scalex+dx,i*scaley+dy,0,200,0);
	}
      }

    }


    // idea of drawing edges for variable M
    for (int i=0;i<q.size();i++){
      for (int j=0;j<domain_size;j++){
	if (q[i].in(j)){
	  for (int dx=0;dx<scalex;dx++){	    
	    for (int k=-1;k<2;k++){
	      int x,y;
	      x=j*scalex-dx+scalex/2+k;
	      y=i*scaley-dx+scalex/2;
	      if (x>0 && y>0 && x<img_x && y<img_y){
		rgb(x,y,0,0,255);
	      }
	    }	    
	  }
	}
      }
      
    }

    // variables G
    for (int i=0;i<q.size();i++){
      for (int j=0;j<domain_size;j++){
	if (q[i].in(j)){
	  for (int dx=0;dx<scalex;dx++){	    
	    for (int k=-1;k<2;k++){
	      int x,y;
	      x=j*scalex   +scalex/2+k;
	      y=i*scaley-dx+scalex/2;
	      if (x>0 && y>0 && x<img_x && y<img_y){
		rgb(x,y,255,0,0);
	      }
	    }	    
	  }
	}
      }      
    }

    // variables H
    for (int i=0;i<q.size();i++){
      for (int j=0;j<domain_size;j++){
	if (q[i].in(j)){
	  for (int dx=0;dx<scalex;dx++){	    
	    for (int k=-1;k<2;k++){
	      int x,y;
	      x=j*scalex-dx+scalex/2;
	      y=i*scaley+   scalex/2+k;
	      if (x>0 && y>0 && x<img_x && y<img_y){
		rgb(x,y,100,100,100);
	      }
	    }	    
	  }
	}
      }      
    }



    display();    
  }

  void update(SetVarArray H){
    const unsigned char green[] = { 64,255,32 }, blue[] = { 128,200,255}, red[] = { 255,0,0 }, white[] = { 255,255,255 };
    // clean image
    for (int i=0;i<img_x;i++){
      for(int j=0;j<img_y;j++){
	rgb(i,j,0,0,0);
      }
    }

    // set variables (red=outside set, grey in the lattice, green in the set
    for (int i=0;i<H.size();i++){
      for (int j=0;j<domain_size;j++){
	int r=100,g=100,b=100;
	
	if (H[i].contains(j)) {g=200;r=0;b=40;}
	if (H[i].notContains(j)) {g=10;r=200;b=10;}
       	
	for (int dx=0;dx<scalex;dx++){	    
	  for (int dy=0;dy<scaley;dy++){	    
	    //	  for (int k=-1;k<2;k++){
	    int x,y;
	    x=j*scalex+dx;
	    y=i*scaley+dy;
	    if (x>0 && y>0 && x<img_x && y<img_y){
	      rgb(x,y,r,g,b);
	    }	    
	  }
	}

	imageOut.draw_text(j*scalex,
			   i*scaley,
			   "i=%d",
			   white,
			   0,1,4,i);

      }      
    }
    display();    
  }
    


  void close(){
    cout << "close (wait 1 sec)\n";
    if (!main_disp.is_empty() && !main_disp.is_closed()) {
      for (int i=0;i<10;i++){
	main_disp.wait(100);
      }
      main_disp.close();
    }
  }


  /**
   * Draws a pixel at x,y with color r g b
   */
  void rgb(int x, int y, int r, int g, int b){
    imageOut(x,y,0)=r;
    imageOut(x,y,1)=g;
    imageOut(x,y,2)=b;
  }

 }; 