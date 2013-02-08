/* 
 * (c) Marek Waligorski 2004, www.icpnet.pl/~emwal
 * You are free to use, distribute and/or modify this code as long as
 * you do it free of charge and mention the original author's name.
 * This code comes WITHOUT WARRANTY of any kind. See LICENSE file for details.
 */

#include "log.h"

//#define DEBUGGING_OFF  		//enable this in the final version for better performance

int C_log::exists = 0; //false

/*------------- Konstruktor & Destruktor ----------------*/
C_log::C_log (const char * nazwa)
{
 this->auto_flush = TRUE;
 this->logFile = NULL; //open=FALSE
 this->error_notify = TRUE;

 this->warning_notify = TRUE;
 #ifdef _WIN32
  this->warning_notify = FALSE;	//don't pop a box on every warning
 #endif

 this->debug = 0;
 num_errors = num_warnings = 0;
   
 if (strlen (nazwa)>30) Notice("Filename too long!"); 
 else
 {
  strcpy (this->filename,nazwa);
  if (exists==1) {Notice ("One log already in use!"); exists=-1;} //error!
  else if (!this->Open ()) Notice("ERROR! can't log!");
  else exists=1; //true
 } 
}

C_log::~C_log ()
{
 this->Close ();
}  

/***************** Open & Close ********************************/
bool C_log::Open ()
{
 if (logFile!=NULL) return FALSE; //prevent from double opening

 if ((logFile = fopen (filename,"wt"))==NULL) return FALSE; //open the file
  
 time_t rawtime; time (&rawtime); //get current time
 struct tm * timeinfo = localtime (&rawtime);
 
 Add ("Starting log: %s",asctime (timeinfo)); 
 Add ("Log class by Marek Waligorski");
 Add ("------------------------------\n");
 
 return TRUE; //log successfully opened
}

void C_log::Close ()
{
 if (logFile==NULL) return; //don't close what ain't open

 fprintf (logFile,"\n----------- END OF FILE ---------------");
 if (exists==-1) fprintf (logFile,"\nCRITICAL LOG ERROR!");

 if (num_errors || num_warnings) fprintf (logFile,"\n%d error(s) and %d warning(s) encountered.",num_errors, num_warnings);

 fclose (logFile);
 this->logFile = NULL; //just to make sure
}

/****************** Add info/error to log *****************************/
void C_log::Add (const char *fmt, ...)
{
 bool isItError = FALSE;
 bool isItWarning = FALSE;
 
 if (exists==-1) return;
 if (logFile==NULL) return;	// do not write to an inexistent file
 if (fmt == NULL) return; 	// If There's No Text, Do Nothing
 
 if (fmt[0]=='!')                         //if it's an error
 {isItError = TRUE; num_errors++;
  fmt++; fprintf(logFile,"<!> ERROR: ");} //write error notification

#ifndef DEBUGGING_OFF 
 else if (fmt[0]=='@') //if it's debug mode
 {if (this->debug==0) return; fmt++; fprintf (logFile,"- ");} //write a hyphen
 
 else if (fmt[0]=='#') //if it's heavy debug information
 {if (this->debug<=1) return; fmt++; fprintf (logFile,"\t");}

 else if (fmt[0]=='$') //if it's extreme debug information (3 or more)
 {if (!this->debug<=2) return; fmt++;fprintf (logFile,"\t\t");}

 else if (fmt[0]=='*')                     		//if it's an warning
  {fmt++; isItWarning=TRUE; fprintf(logFile,"<*> WARNING: ");
   num_warnings++;} 	

#endif //DEBUGGING_OFF

 else fprintf (logFile, "* "); //write a star for important messages
 
 char		text[256];			// Holds Our String
 va_list		ap;			// Pointer To List Of Arguments
 va_start(ap, fmt);				// Parses The String For Variables
   vsprintf(text, fmt, ap);	 		// And Converts Symbols To Actual Numbers
 va_end(ap);					// Results Are Stored In Text
 
 fprintf (logFile, "%s\n", text);
 
 if (auto_flush) fflush (this->logFile);
 
 if (error_notify && isItError) Notice ("ERROR: %s\n",text);  //print error message on the screen
 if (warning_notify && isItWarning) Notice ("Warning: %s\n",text);
}
