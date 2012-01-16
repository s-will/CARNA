#ifndef WIN_HANDLER_HH
#define WIN_HANDLER_HH

#include <iostream>

#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include <gecode/set.hh>
#include <CImg.h>

#include <semaphore.h>
#include <pthread.h>


class WinHandler
{

    pthread_t m_thread,w_thread;
    sem_t s1,lock;  /// semaphores

    cimg_library::CImg<unsigned char> imageOut; 
    cimg_library::CImgDisplay main_disp;
    const char* title;

    int nRows;
    int nCols;

    int img_x;
    int img_y;  

public:
    WinHandler(int Rows, int Cols, std::string t);
  
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
    
    
    void update(const Gecode::IntVarArray &MD, const Gecode::BoolVarArray &M);
    
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

    void window_thread();

    void do_work();
 
};

#endif // WIN_HANDLER_HH
