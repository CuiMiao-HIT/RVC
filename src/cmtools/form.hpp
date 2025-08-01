#ifndef FORM_HPP_
#define FORM_HPP_

#include <fstream>
#include <getopt.h>
#include <string.h>

#include "desc.hpp"

#define PATH_LEN 1024

typedef struct{
    int argc;
    char **argv;
    
// Parameters of input
    char chrID[PATH_LEN];
    int Qual;
    char snvinfo_path[PATH_LEN];
    char bed_path[PATH_LEN];
    char snvvcf_path[PATH_LEN];
        
}opt_t;

int opt_parse(int argc, char *argv[], opt_t* p);


#endif /* FORM_HPP_ */