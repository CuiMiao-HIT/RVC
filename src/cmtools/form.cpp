/*
 * form.cpp
 *
 *      Author: miaocui
 */
#include "form.hpp"

int help_usage() {
    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   %s, variant caller\n", PACKAGE_NAME); 
	fprintf(stderr, "Usage:	%s caller [options] -c <chromosome> -i <fin_file> -o <output.vcf>\n\n", PACKAGE_NAME);

    fprintf(stderr, "    -c --chr [chr22]                             The chromosome ID.\n");
    fprintf(stderr, "    -q --qual [default 1]                                Variant quality for filtering.\n");
    fprintf(stderr, "    -i --fin_file <file>                         The input vatiants information format.\n");
    fprintf(stderr, "    -b --bed_file <file>                         The candidates bed format.\n");
    fprintf(stderr, "    -o --fout_vcf <vcf>                          The output vcf format.\n");


    fprintf(stderr, "\n");
    
    return 1; 
}

char *const short_options = "c:i:o:b:q:";
struct option long_options[] = {
    { "chr", 1, NULL, 'c'},
    { "qual", 1, NULL, 'q'},
    { "fin_file", 1, NULL, 'i'},
    { "bed_file", 1, NULL, 'b'},
    { "fout_vcf", 1, NULL, 'o'},
    { 0, 0, 0, 0}
};

int opt_parse(int argc, char *argv[], opt_t *opt) {
	int optc; 
    char *p;
    int option_index=0;
    while((optc = getopt_long(argc, argv, short_options, long_options, &option_index))>=0){ 
        switch(optc){
            case 'c': strcpy(opt->chrID, optarg); break;
            case 'q': opt->Qual = atoi(optarg); break;
            case 'i': strcpy(opt->snvinfo_path, optarg); break;
            case 'b': strcpy(opt->bed_path, optarg); break;
            case 'o': strcpy(opt->snvvcf_path, optarg); break;
            default:
                fprintf(stderr,"Wrong parameters\n");
                return help_usage();
        }
    }
    if(optind +1 != argc) {
        return help_usage();
    }
    opt->argv = argv;
    opt->argc = argc;
    optind++;
    return 0;  
}