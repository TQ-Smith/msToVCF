
// File: Main.c
// Date: 18 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Convert ms style output to VCF format.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include "../lib/ketopt.h"
#include "../lib/zlib.h"
#include "../lib/kstring.h"
#include "../lib/kseq.h"

// Checks that user supplied options are valid.
// Accepts:
//  int ploidy -> The user supplied ploidy.
//  int length -> The user supplied length.
//  double missing -> The user supplied missing genotype probability.
// Returns:
//  int, 0 or 1, for valid user supplied options or invalid options, respectively.
int check_configuration(int ploidy, int length, double missing) {
    if (ploidy < 2) {
        printf("Error! Ploidy must be 2 or greater.\n");
        return 1;
    }
    if (length < 1000) {
        printf("Error! Length must be 1000 or greater to avoid multiple records at the same locus.\n");
        return 1;
    }
    if (missing < 0 || missing >= 1) {
        printf("Error! The probability of a missing genotype must be in [0, 1).\n");
        return 1;
    }
    return 0;
}

// Print the help menu for msToVCF.
// Accepts: void.
// Returns: void.
void print_help() {
    printf("\n");
    printf("msToVCF v1.0 December 2024\n");
    printf("----------------------\n\n");
    printf("Written by T. Quinn Smith\n");
    printf("Principal Investigator: Zachary A. Szpiech\n");
    printf("The Pennsylvania State University\n\n");
    printf("Usage: msToVCF [options]\n");
    printf("Options:\n");
    printf("   -h/--help                  Prints help menu and exits.\n");
    printf("   -v/--version               Prints version number and exits.\n");
    printf("   -p/--ploidy INT            Sets ploidy of samples. Default 2.\n");
    printf("   -l/--length INT            Sets length of segment in number of base pairs. Default 1,000,000.\n");
    printf("   -u/--unphased              If set, the phase is removed from genotypes.\n");
    printf("   -m/--missing DOUBLE        Genotypes are missing with supplied probability. Default 0.\n");
    printf("   -c/--compress              If set, the resulting files are gzipped compressed.\n");
    printf("\n");
}

// Long options used for LODESTAR.
static ko_longopt_t long_options[] = {
    {"help",            ko_no_argument,         'h'},
    {"version",         ko_no_argument,         'v'},
    {"ploidy",          ko_required_argument,   'p'},
    {"length",          ko_required_argument,   'l'},
    {"unphased",        ko_no_argument,         'u'},
    {"missing",         ko_required_argument,   'm'},
    {"compress",        ko_no_argument,         'c'},
    {0, 0, 0}
};

int main(int argc, char *argv[]) {

    // Seed random number generator.
    srand(time(NULL));

    // Single character aliases for long options.
    const char *opt_str = "h:v:p:l:u:m:c";
    ketopt_t options = KETOPT_INIT;
    int c;

    // Print help menu when no options/arguments were given.
    if (argc == 1) {
        print_help();
        return 1;
    }        

    // Pass through options. Check for options requiring an argument that were not given one.
    //  Also, check if any options are supplied that are not defined. If user supplies the help
    //  argument, print help menu and exit, likewise for version.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': printf("Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return 1;
            case '?': printf("Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return 1;
            case 'h': print_help(); return 0;
            case 'v': printf("Version 1.0 December 2024.\n"); return 0;
        }
	}
	options = KETOPT_INIT;

    // Set default option values.
    int ploidy = 2;
    int length = 1000000;
    bool unphased = false;
    double missing = 0;
    bool compress = false;

    // Parse command line arguments.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'p':  break;
            case 'l':  break;
            case 'u':  break;
            case 'm':  break;
            case 'c':  break;
        }
	}
    
    // Check configuration. If invalid argument, exit program.
    if (check_configuration(ploidy, length, missing) != 0) {
        printf("Exiting!\n");
        return 1;
    }

}