#include "interface/main.hpp"
#include "interface/plugin.hpp"
#include <signal.h>

using namespace theta;

REGISTER_PLUGIN_BASETYPE(Main);

bool theta::stop_execution = false;

static void sigint_handler(int sig){
   if(stop_execution){
      throw ExitException("user insisted with several SIGINT");
   }
   stop_execution = true;
}

void theta::install_sigint_handler(){
    struct sigaction siga;
    siga.sa_handler = sigint_handler;
    sigemptyset(&siga.sa_mask);
    siga.sa_flags = 0;
    sigaction(SIGINT, &siga, 0);
}


