#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdlib>
#include <cstdio>
#include <cstring> // for memset
#include <unistd.h>
#ifndef CH_DISABLE_SIGNALS
#include <csignal>
#endif


#ifdef CH_Darwin
# define _Atomic
#include <libproc.h>
#endif

#include "SPMD.H"
#include "CH_Attach.H"
#include "parstream.H"
#include "BaseNamespaceHeader.H"

int pfds[2];

bool alreadyAttached = false;
using std::endl;
void DebugCont()
{
// #ifndef CH_DISABLE_SIGNALS
//   signal(SIGSEGV, SIG_DFL);
//   signal(SIGABRT, SIG_DFL);
// #endif
  ssize_t ret;
  ret = write(pfds[1], "OK", 3);
}

#ifndef CH_NTIMER
extern const char* currentTimer();
#else
const char* currentTimer(){ static const char* rtn="empty"; return rtn;}
#endif

void AttachDebugger(int a_sig)
{
#ifndef CH_NTIMER
  pout()<<"current Timer " <<currentTimer()<<"\n";
#endif

#ifndef CH_DISABLE_SIGNALS
  if (alreadyAttached) return;

  alreadyAttached = true;
  char buf[1200];
  char* display;
  FILE *f;
  int proc = getpid();
  char title[1024];

  char binaryName[1025];
  memset(binaryName, 0, 201); // readlink doesn't insert NULL at end

  char linkName[1024];

#ifdef CH_Linux   // note that a machine can be Linux, but not support /proc (Catamount)
  // voodoo specific to linux /proc system to get the full name of the binary file
  sprintf(linkName,"/proc/%d/exe",proc);
  ssize_t s = readlink(linkName, binaryName, 1024);
#else
  ssize_t s = proc_pidpath (proc, binaryName, 1024);
  std::cout<<"s: "<<s<<"\n";
  std::cout<<binaryName<<std::endl;
  
#endif

  int a = 50; // used for title of emacs window

  if (s < a)
  {
    a = s;
  }

  int ret;
  ret = pipe(pfds);

  //System V vs. POSIX naming differences for child signal handler.
#ifdef SIGCHLD
#ifndef SIGCLD
#define SIGCLD SIGCHLD
#endif
#endif

  signal(SIGCLD, SIG_IGN);
  if (!fork())
  {
    // try to find a display to send this debugger to.
    display = getenv("DISPLAY");
    printf("%s\n%s\n",linkName, display);

    if (display != NULL)
      {
        sprintf(title, "'%s %d'", binaryName+s-a, procID());

#define useemacs 1
#if useemacs
        //  sprintf(buf,"emacs -geometry 80x27 -title %s -display %s --eval '(progn (gdb \" gdb -i=mi -q -f %s --pid %d\") )'",
        //          title, display, binaryName, proc);
        sprintf(buf,"emacs -geometry 80x27 -title %s -display %s --eval '(progn (lldb \" lldb -p %d\") )'",
                title, display, proc);
#else
        sprintf(buf,"xterm -e \"gdb -q -f %s --pid %d\" ", binaryName, proc);
#endif
        printf("%s\n",buf);
        f = popen(buf,"w");
      }
    else
      {
        sprintf(buf,"gdb -q  %s %d", binaryName, proc);
        f=popen(buf, "w");
        fputs("where\ncall DebugCont()\ndetach\nquit\n", f);
        fflush(f);
        sleep(2);
      }

    pclose(f); // child process stalls
    exit(0);
  }
  else
  {
    char rbuf[3];
    ssize_t ret;
    ret = read(pfds[0], rbuf, 3);
  }
#endif
}

#ifdef CH_MPI
#include "SPMD.H"
#endif

// when running the code inside an actual debugger, the debugger signal handler
// will come in before this one, so the debugger should function correctly.

extern bool maydayabortFlag;

int registerDebugger()
{
  int rtn = 0;
  maydayabortFlag = true;

#ifndef CH_DISABLE_SIGNALS
  signal(SIGSEGV, AttachDebugger);
  signal(SIGABRT, AttachDebugger);

#ifdef CH_MPI
  int r = setChomboMPIErrorHandler();
  rtn+=r;
#endif
#else
  rtn = 2;
#endif
  return rtn;
}

#ifdef CH_MPI

void MPI_Shutdown()
{
  MPI_Abort(MPI_COMM_WORLD, -4);
}

void mpierrorfunction(MPI_Comm* a_comm, int* a_error, ...)
{
  char string[MPI_MAX_ERROR_STRING];
  int resultlen;
#ifndef CH_NTIMER
  pout() << "current Timer " << currentTimer() << endl;
#endif

  MPI_Error_string(*a_error, string, &resultlen);
  pout() << string << endl;
  //printf("%s\n",string);

  MayDay::Error("Chombo MPI Error handler called \n");

}

int setChomboMPIErrorHandler()
{
  int rtn=1;

  MPI_Errhandler handler;

  // Some MPI implementations are having trouble with this.  In particular,
  // halem at GSFC.  For now, i'm just going to not do this for machines
  // running OSF (which currently, is only halem for us) (ndk)
#ifndef CH_OSF1
  // deprecated API
  //MPI_Errhandler_create(mpierrorfunction, &handler);
  //rtn = MPI_Errhandler_set(Chombo_MPI::comm, handler);
  MPI_Comm_create_errhandler(mpierrorfunction, &handler);
  rtn = MPI_Comm_set_errhandler(Chombo_MPI::comm, handler);
  if (rtn == MPI_SUCCESS) rtn = 0;
#else
  rtn = 2;
#endif

  return rtn;
}
#endif

#include "BaseNamespaceFooter.H"
