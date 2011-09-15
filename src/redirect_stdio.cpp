#include "interface/redirect_stdio.hpp"

#include <sys/stat.h>
#include <fcntl.h>


// as default: these are synonyms for stdout, stderr:
int theta::cout_fd = 1;
int theta::cerr_fd = 2;
boost::iostreams::stream<boost::iostreams::file_descriptor_sink> theta::cout(boost::iostreams::file_descriptor_sink(1, false));
boost::iostreams::stream<boost::iostreams::file_descriptor_sink> theta::cerr(boost::iostreams::file_descriptor_sink(2, false));

void theta::redirect_stdio(){
    static bool already_replaced = false;
    if(!already_replaced){
        theta::cout_fd = dup(1);
        theta::cout.close();
        theta::cout.open(boost::iostreams::file_descriptor_sink(theta::cout_fd, false));
        close(1);
        open("/dev/null", O_WRONLY);
        theta::cerr_fd = dup(2);
        theta::cerr.close();
        theta::cerr.open(boost::iostreams::file_descriptor_sink(theta::cerr_fd, false));
        close(2);
        open("/dev/null", O_WRONLY);
        already_replaced = true;
    }
}


