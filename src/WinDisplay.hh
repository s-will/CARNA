#include <iostream>

#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include <gecode/set.hh>
#include <CImg.h>


using namespace cimg_library;
using namespace std;

class WinDisplay {

  CImg<unsigned char> imageOut; 
  CImgDisplay main_disp;

  int nRows;
  int nCols;
  int scale;
  int img_x;
  int img_y;  

  int undef;

  const char* title;

  


public: 
  WinDisplay(int Rows, int Cols){   
    WinDisplay(Rows,Cols,"Untitled");
  }

  WinDisplay(int Rows, int Cols, string t){   
    int min_win_size=1000;

    undef= Rows;

    title=t.c_str();
    nRows=Rows;
    nCols=Cols;
    int scalex=1;
    int scaley=1;
    img_x=0;
    img_y=0;    
    if (Cols<min_win_size){
      scalex=min_win_size/Cols;
    }
    if (Rows<min_win_size){
      scaley=min_win_size/Rows;
    }
    
    if (scalex<scaley) scale=scalex;
    else
      scale=scaley;
    img_x=Cols*scale;
    img_y=Rows*scale;

    imageOut = CImg<>(img_x,img_y,1,3,0);
    cout << "allocating image: " << img_x << " " << img_y << "\n";
    
    for (int i=0;i<img_x;i++)
      for(int j=0;j<img_y;j++)
	imageOut(i,j,1)=200;   //uninitialized
        
        main_disp=CImgDisplay();
	//main_disp.close();
   }

  ~WinDisplay(){
    cout << "calling destructor\n";
    close();
    //delete imageOut;
    //delete main_disp;
  }

  void display(){

    if (main_disp.is_empty() || main_disp.is_closed())
      //if (main_disp.is_closed) // the above line is not compiling (change in current CImg version?. What is correct, Ale?
      main_disp=CImgDisplay(imageOut,title,0);
    main_disp.display(imageOut);  
  }

    void update(const Gecode::IntVarArray &MD, const Gecode::BoolVarArray &M){
	//const unsigned char green[] = { 64,255,32 }, blue[] = { 128,200,255}, red[] = { 255,0,0 }, white[] = { 255,255,255 };
        printf("OK: imgx %d, imgy %d, col %d, row %d, scale %d\n",img_x,img_y,nCols,nRows,scale);
    // clean image
    for (int i=0;i<img_x;i++){
      for(int j=0;j<img_y;j++){
	rgb(i,j,0,0,0);
      }
    }

    /*    // for each domain value paint green block if it is in the current domain
    for (int i=0;i<q.size();i++){
      for (int j=0;j<domain_size;j++){

	if (q[i].in(j)){
	  for (int dx=0;dx<scalex;dx++)
	    for (int dy=0;dy<scaley;dy++)

	      rgb(j*scalex+dx,i*scaley+dy,0,200,0);
	}
      }

    }
    */


    // idea of drawing edges for matches
    for (int i=0;i<nRows;i++){
      for (int j=0;j<nCols;j++){
	  if (MD[i].in(j) && M[i].in(1)){
	  for (int dx=0;dx<scale;dx++){	    
	    for (int k=-1;k<2;k++){
	      int x,y;
	      x=j*scale-dx+scale/2+k;
	      y=i*scale-dx+scale/2;
	      if (x>0 && y>0 && x<img_x && y<img_y){
		  if (!M[i].assigned())
		      rgb(x,y,0,128,255); // if undecided -> cyan color
		  else
		      rgb(x,y,0,0,255);   // blue if chosen
	      }
	    }	    
	  }
	}
      }
      
    }

    // deletions
    for (int i=0;i<nRows;i++){
      for (int j=0;j<nCols;j++){
	if (MD[i].in(j) && M[i].in(0)){
	  for (int dx=0;dx<scale;dx++){	    
	    for (int k=-1;k<2;k++){
	      int x,y;
	      x=j*scale   +scale/2+k;
	      y=i*scale-dx+scale/2;
	      if (x>0 && y>0 && x<img_x && y<img_y){
		  if (!M[i].assigned())
		      rgb(x,y,255,128,0); // if undecided -> orange color
		  else
		      rgb(x,y,255,0,0);  // if can not be undef -> red
	      }
	    }	    
	  }
	}
      }      
    }

    // insertions
    for (int i=0;i<nRows;i++){
      for (int j=0;j<nCols;j++){
	int r=100,g=100,b=100;
	
	int minj=nCols-1;
	if (i+1<nRows) {
	    minj = MD[i+1].min()-(!M[i+1].in(0)?1:0);
	}
	
	int maxj=nCols-1;
	if (i+1<nRows) {
	    maxj = MD[i+1].max()-(!M[i+1].in(0)?1:0);
	}
	
	if (MD[i].max()<j && j<=minj) {g=200;r=0;b=40;}
	if (MD[i].min()<j && j<=maxj) { //{g=10;r=200;b=10;}
	  for (int dx=0;dx<scale;dx++){	    
	    for (int k=-1;k<2;k++){
	      int x,y;
	      x=j*scale-dx+scale/2;
	      y=i*scale+   scale/2+k;
	      if (x>0 && y>0 && x<img_x && y<img_y){
		rgb(x,y,r,g,b);
	      }
	    }	    
	  }
	}
      }      
    }


    /*
      // keep for future use: text drawing in the picture
    for (int i=0;i<nRows;i++){
      imageOut.draw_text(0*scalex,
			 i*scaley,
			 "i=%d",
			 white,
			 0,1,4,i);
    }
    */

    display();    
  }
  
  void close(){
    cout << "close (wait 1 sec)\n";
    if (!main_disp.is_empty() && !main_disp.is_closed()) {
      //if (! main_disp.is_closed) // the above line is not compiling (change in current CImg version?. What is correct, Ale? 
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

