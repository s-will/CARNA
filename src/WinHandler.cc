#include "config.h"

#ifdef HAVE_GIST

#include "WinHandler.hh"

WinHandler::WinHandler(int Rows, int Cols, std::string t)
    : m_stoprequested(false), m_running(false)
{
    nRows=Rows;
    nCols=Cols;

    title=t.c_str();

    img_x=Cols;
    img_y=Rows;

    
    imageOut = cimg_library::CImg<>(img_x,img_y,1,3,0);
    //cout << "allocating image: " << img_x << " " << img_y << "\n";
    
    for (int i=0;i<img_x;i++)
	for(int j=0;j<img_y;j++)
	    imageOut(i,j,1)=128;   //uninitialized

    go();
}

void WinHandler::update(const Gecode::IntVarArray &MD, const Gecode::BoolVarArray &M) {
    //printf("Update request\n");
    //printf("OK: imgx %d, imgy %d, col %d, row %d\n",img_x,img_y,nCols,nRows);
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
	    if (MD[i].in(j) && M[i].in(1)
		&& (i==0 || j>MD[i-1].min())
		){
		int x,y;
		x=j;//*scale-dx+scale/2+k;
		y=i;//*scale-dx+scale/2;
		if (x>0 && y>0 && x<img_x && y<img_y){
		    rgbAdd(x,y,0,0,255);   // blue
		}
	    }
	}
      
    }
    // deletions
    for (int i=0;i<nRows;i++){
	for (int j=0;j<nCols;j++){
	    if (MD[i].in(j) && M[i].in(0)){
		int x,y;
		x=j;//*scale   +scale/2+k;
		y=i;//*scale-dx+scale/2;
		if (x>0 && y>0 && x<img_x && y<img_y){
		    rgbAdd(x,y,255,0,0);  // red
		}
	    }
	}
    }
    // insertions
    for (int i=0;i<nRows;i++){
	for (int j=0;j<nCols;j++){
	
	    int maxj=nCols-1;
	    if (i+1<nRows) {
		maxj = MD[i+1].max()-(!M[i+1].in(0)?1:0);
	    }
	
	    if (MD[i].min()<j && j<=maxj) {
		int x,y;
		x=j;//*scale-dx+scale/2;
		y=i;//*scale+   scale/2+k;
		if (x>0 && y>0 && x<img_x && y<img_y){
		    rgbAdd(x,y,0,255,0); //green
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


    sem_post(&s1);   
}


void WinHandler::window_thread(){
    struct timespec   ts = {0, 0};

    ts.tv_sec  = 0;
    ts.tv_nsec = 10000000;

    while (m_stoprequested==false){
	//      sem_wait(&lock);
	if (!main_disp.is_closed()){
	    //	sem_post(&lock);
	    if (main_disp.is_resized()){
		//printf("resize\n");
		main_disp.resize().display(imageOut);
	    }
	    //	printf("waiting for a window event\n");
	    nanosleep(&ts,NULL);
	    //CImgDisplay::wait(main_disp);
	    //printf("window event received\n");
	}
	//      else
	//	sem_post(&lock);

      
    }
}


void WinHandler::do_work()
{
    //printf("Working\n");
    main_disp=cimg_library::CImgDisplay(imageOut,title,0);

    
    // --------------------
    // here we can set the initial display size
    // set image dimensions in order to fit minimal and maximal 
    // dimensions
    
    // set dimensions such that each is at least 300 and at least
    // twice sequence length+1
    // except this exceeds some maximal dimensions

    size_t wdim_x=img_x;
    size_t wdim_y=img_y;
    
    const size_t mindim=std::max(300,std::min(2*img_x,2*img_y));
    const size_t maxdim_x=1000;
    const size_t maxdim_y=750;
    
    double ratio=1;
    
    if ((size_t)std::min(wdim_x,wdim_y) < mindim) {
	ratio = ((double)mindim/std::min(wdim_x,wdim_y));
    }
    
    wdim_x *= ratio;
    wdim_y *= ratio;
    
    ratio=1;

    if ((size_t)wdim_x>maxdim_x) {
	ratio = ((double)maxdim_x/wdim_x);
    }
    if ((size_t)wdim_y*ratio>maxdim_y) {
	ratio = ((double)maxdim_y/wdim_y);
    }
    wdim_x *= ratio;
    wdim_y *= ratio;

    main_disp.resize(wdim_x,wdim_y,true);
    // end resizing
    // --------------------


    main_disp.display(imageOut);
    
    
    
    
    while (m_stoprequested==false){
	// wait for something to do
	sem_wait(&s1);
	//printf("update requested\n");
	//      sem_wait(&lock);
	//      if (main_disp.is_empty()){
	//	printf("Empty display -> link image\n");
	main_disp.display(imageOut);  
	//      }
	//sem_post(&lock);
      
	// sem_wait(&lock);
	if (main_disp.is_closed()){
	    //printf("Window closed -> show\n");
	    main_disp.show();
	    main_disp.resize(wdim_x,wdim_y,true);
	}
	//sem_post(&lock);
      
	main_disp.paint();
      
	//      printf("win closed: %d\n",main_disp.is_closed());

    }

}                    

#endif //HAVE_GIST
