#include <iostream>

#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include <gecode/set.hh>
#include <CImg.h>

#include <semaphore.h>
#include <pthread.h>

using namespace cimg_library;
using namespace std;

class WinHandler
{

  pthread_t m_thread,w_thread;
  sem_t s1,lock;  /// semaphores

  CImg<unsigned char> imageOut; 
  CImgDisplay main_disp;
  const char* title;

  int nRows;
  int nCols;

  int img_x;
  int img_y;  

public:
  WinHandler(int Rows, int Cols, string t)
    : m_stoprequested(false), m_running(false)
  {
    nRows=Rows;
    nCols=Cols;

    title=t.c_str();

    img_x=Cols;
    img_y=Rows;

    imageOut = CImg<>(img_x,img_y,1,3,0);
    cout << "allocating image: " << img_x << " " << img_y << "\n";
    
    for (int i=0;i<img_x;i++)
      for(int j=0;j<img_y;j++)
	imageOut(i,j,1)=128;   //uninitialized

    go();
  }
  
  ~WinHandler()
  {

  }
  
  // Create the thread and start work
  void go() 
  {

    sem_init(&s1,0,1); // initialize at 1
    sem_init(&lock,0,1); // initialize at 1

    m_running = true;
    //    printf("create thread (dispatcher)\n");
    pthread_create(&m_thread, 0, &WinHandler::start_dispatcher_thread, this);
    //    printf("create thread (window)\n");
    pthread_create(&w_thread, 0, &WinHandler::start_window_thread, this);
  }
  
  void stop() // Note 2
  {
    m_running = false;
    m_stoprequested = true;
    pthread_join(m_thread, 0);
    pthread_join(w_thread, 0);
  }
  
  void update(const Gecode::IntVarArray &MD, const Gecode::BoolVarArray &M){
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
	if (MD[i].in(j) && M[i].in(1)){
	  int x,y;
	  x=j;//*scale-dx+scale/2+k;
	  y=i;//*scale-dx+scale/2;
	  if (x>0 && y>0 && x<img_x && y<img_y){
	    if (!M[i].assigned())
	      rgbAdd(x,y,0,0,255); // if undecided -> cyan color
	    else
	      rgbAdd(x,y,0,0,255);   // blue if chosen
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
	    if (!M[i].assigned())
	      rgbAdd(x,y,255,0,0); // if undecided -> orange color
	    else
	      rgbAdd(x,y,255,0,0);  // if can not be undef -> red
	  }
	}
      }      
    }
    // insertions
    for (int i=0;i<nRows;i++){
      for (int j=0;j<nCols;j++){
	int r=0,g=255,b=0;
	
	int minj=nCols-1;
	if (i+1<nRows) {
	  minj = MD[i+1].min()-(!M[i+1].in(0)?1:0);
	}
	
	int maxj=nCols-1;
	if (i+1<nRows) {
	  maxj = MD[i+1].max()-(!M[i+1].in(0)?1:0);
	}
	
	if (MD[i].max()<j && j<=minj) {g=255;r=00;b=0;}
	if (MD[i].min()<j && j<=maxj) { //{g=10;r=200;b=10;}
	  int x,y;
	  x=j;//*scale-dx+scale/2;
	  y=i;//*scale+   scale/2+k;
	  if (x>0 && y>0 && x<img_x && y<img_y){
	    rgbAdd(x,y,r,g,b);
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

  void rgb(int x, int y, int r, int g, int b){
    imageOut(x,y,0)=r;
    imageOut(x,y,1)=g;
    imageOut(x,y,2)=b;
  }

  void rgbAdd(int x, int y, int r, int g, int b){
    int t;
    t=imageOut(x,y,0)+r;if (t>255) t=255;
    imageOut(x,y,0)=t;

    t=imageOut(x,y,1)+g;if (t>255) t=255;
    imageOut(x,y,1)=t;
    t=imageOut(x,y,2)+b;if (t>255) t=255;
    imageOut(x,y,2)=t;
  }

  
private:
  volatile bool m_stoprequested; // Note 5
  volatile bool m_running;
 
  // This is the static class function that serves as a C style function pointer
  // for the pthread_create call
  static void* start_dispatcher_thread(void *obj)
  {
    //All we do here is call the do_work() function
    reinterpret_cast<WinHandler *>(obj)->do_work();
    return obj;
  }

  static void* start_window_thread(void *obj)
  {
    //All we do here is call the do_work() function
    reinterpret_cast<WinHandler *>(obj)->window_thread();
    return obj; 
  }

  void window_thread(){
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
 
  void do_work()
  {
    //printf("Working\n");
    main_disp=CImgDisplay(imageOut,title,0);
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
      }
      //sem_post(&lock);
      
      main_disp.paint();
      
      //      printf("win closed: %d\n",main_disp.is_closed());

    }

  }                    
};
