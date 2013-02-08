/* 
 * (c) Marek Waligorski 2004, www.icpnet.pl/~emwal
 * You are free to use, distribute and/or modify this code as long as
 * you do it free of charge and mention the original author's name.
 * This code comes WITHOUT WARRANTY of any kind.
 */

/** Log Class **/

#ifndef _MW_LOG_CLASS_
#define _MW_LOG_CLASS_

#include <fstream>
#include <time.h>
#include <stdarg.h>
#include <cstring>

//a nice shortcut (macro)
#define LOG Log.Add	

#ifdef _WIN32_ 	//define the message procedure
 #include "E_aux_func.hpp"
 #define Notice MBprint
#else
 #define Notice printf
#endif

#define NO_DEBUG 	0
#define LIGHT_DEBUG 	1
#define HEAVY_DEBUG 	2
#define EXTREME_DEBUG 	3

using namespace std;

class C_log {
 private:
  char filename [20];
  FILE* logFile;
  static int exists;
  int num_errors;
  int num_warnings;

 public:
  C_log (const char * nazwa);
  ~C_log ();  
  bool Open ();
  void Close ();
  void Add (const char * fmt, ...);

  bool auto_flush;
  bool warning_notify;
  bool error_notify;
  unsigned int debug; //debug level (0- no debug, 3- extreme debug)
};

extern C_log Log;

#endif //_MW_LOG_CLASS_
