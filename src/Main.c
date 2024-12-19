
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
#include "../lib/kvec.h"

// We use kseq to read in from stdin.
#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

// Prints ms replicate to VCF file.
// Accepts:
//  char* fileName -> The name of the input file.
//  int length -> The length of the segment in bp.
//  bool unphased -> If set, the resulting output should be unphased.
//  bool compress -> If set, the resulting files should be compressed.
//  int numReplicate -> The current replicate number.
//  int numSegsites -> The number of segregating sites in the replicate.
//  int numSamples -> The number of samples in the replicate.
//  kvec_t(double) positions -> The list of segregating site positions.
//  kvec_t(kstring_t*) samples -> The list of simulated samples.
// Returns: void.
void toVCF(char* fileName, int length, bool unphased, double missing, bool compress, int numReplicate, int numSegsites, int numSamples, kvec_t(double)* positions, kvec_t(kstring_t*)* samples) {

    // Create the output file name.
    kstring_t outputFileName;
    if (strncmp(fileName + strlen(fileName) - 3, ".ms", 3) == 0) {
        kputsn(fileName, &outputFileName, strlen(fileName) - 3);
    }
    if (strncmp(fileName + strlen(fileName) - 6, ".ms.gz", 6) == 0) {
        kputsn(fileName, &outputFileName, strlen(fileName) - 6);
    }

    // Write to the file
    int prevPosition = 0, pos, leftGeno, rightGeno, temp;
    if (!compress) {
        kputs(sprintf("_rep%d.vcf", numReplicate), &outputFileName);
        FILE* fp = open(outputFileName.s, "w");
        fprintf(fp, "##fileformat=VCFv4.2\n");
        fprintf(fp, "##contig=<ID=chr1,length=%d>\n", length);
        fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (int i = 0; i < numSamples / 2; i++) {
            fprintf(fp, "\ts%d", i);
        } 
        fprintf(fp, "\n");
        for (int i = 0; i < numSegsites; i++) {
            pos = kv_A(*positions, i);
            if ((int) (pos * length) == prevPosition) {
                pos += 1;
            }
            prevPosition = pos;
            fprintf(fp, "chr1\t%d\t.\tA\tT\t.\t.\t.\t.", pos);
            for (int j = 0; j < numSamples / 2; j++) {
                leftGeno = kv_A(*samples, 2 * j) -> s[i];
                rightGeno = kv_A(*samples, 2 * j + 1) -> s[i];
                if (unphased) {
                    if (rand() < 0.5) {
                        temp = leftGeno;
                        leftGeno = rightGeno;
                        rightGeno = leftGeno;
                    }
                    if (missing > 0) {
                        if (rand() < missing) fprintf(fp, "\t."); else fprintf(fp, "\t%d", leftGeno);
                        if (rand() < missing) fprintf(fp, "/."); else fprintf(fp, "/%d", rightGeno);
                    } else {
                        fprintf(fp, "\t%d/%d", leftGeno, rightGeno);
                    }
                } else {
                    if (missing > 0) {
                        if (rand() < missing) fprintf(fp, "\t."); else fprintf(fp, "\t%d", leftGeno);
                        if (rand() < missing) fprintf(fp, "|."); else fprintf(fp, "|%d", rightGeno);
                    } else {
                        fprintf(fp, "\t%d|%d", leftGeno, rightGeno);
                    }
                }
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    } else {
        kputs("_rep", &outputFileName);
        ks_resize(&outputFileName, strlen(outputFileName.s) + ((int) log10(numReplicate) + 1));
        snprintf(outputFileName.s, strlen(outputFileName.s), "%d", numReplicate);
        kputs(".vcf.gz", &outputFileName);
        gzFile fp = gzopen(outputFileName.s, "w");
        gzprintf(fp, "##fileformat=VCFv4.2\n");
        gzprintf(fp, "##contig=<ID=chr1,length=%d>\n", length);
        gzprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (int i = 0; i < numSamples / 2; i++) {
            gzprintf(fp, "\ts%d", i);
        } 
        gzprintf(fp, "\n");
        for (int i = 0; i < numSegsites; i++) {
            pos = kv_A(*positions, i);
            if ((int) (pos * length) == prevPosition) {
                pos += 1;
            }
            prevPosition = pos;
            gzprintf(fp, "chr1\t%d\t.\tA\tT\t.\t.\t.\t.", pos);
            for (int j = 0; j < numSamples / 2; j++) {
                leftGeno = kv_A(*samples, 2 * j) -> s[i];
                rightGeno = kv_A(*samples, 2 * j + 1) -> s[i];
                if (unphased) {
                    if (rand() < 0.5) {
                        temp = leftGeno;
                        leftGeno = rightGeno;
                        rightGeno = leftGeno;
                    }
                    if (missing > 0) {
                        if (rand() < missing) gzprintf(fp, "\t."); else gzprintf(fp, "\t%d", leftGeno);
                        if (rand() < missing) gzprintf(fp, "/."); else gzprintf(fp, "/%d", rightGeno);
                    } else {
                        gzprintf(fp, "\t%d/%d", leftGeno, rightGeno);
                    }
                } else {
                    if (missing > 0) {
                        if (rand() < missing) gzprintf(fp, "\t."); else gzprintf(fp, "\t%d", leftGeno);
                        if (rand() < missing) gzprintf(fp, "|."); else gzprintf(fp, "|%d", rightGeno);
                    } else {
                        gzprintf(fp, "\t%d|%d", leftGeno, rightGeno);
                    }
                }
            }
            gzprintf(fp, "\n");
        }
        gzclose(fp);
    }
    free(outputFileName.s);
}

// Checks that user supplied options are valid.
// Accepts:
//  int length -> The user supplied length.
//  double missing -> The user supplied missing genotype probability.
// Returns:
//  int, 0 or 1, for valid user supplied options or invalid options, respectively.
int check_configuration(int length, double missing) {
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
    printf("Usage: msToVCF <inFile.ms.gz> [options]\n");
    printf("Options:\n");
    printf("   -h/--help                  Prints help menu and exits.\n");
    printf("   -v/--version               Prints version number and exits.\n");
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
    const char *opt_str = "h:v:l:u:m:c";
    ketopt_t options = KETOPT_INIT;
    int c;    

    // File name is the first argument.
    char* fileName = argv[1];
    argv += 2;

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
    int length = 1000000;
    bool unphased = false;
    double missing = 0;
    bool compress = false;

    // Parse command line arguments.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'l': length = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'u': unphased = true; break;
            case 'm': missing = strtod(options.arg, (char**) NULL); break;
            case 'c': compress = true; break;
        }
	}
    
    // Check configuration. If invalid argument, exit program.
    if (check_configuration(length, missing) != 0) {
        printf("Exiting!\n");
        return 1;
    }

    // If the file does not have .ms or .ms.gz extension.
    if (strncmp(fileName + strlen(fileName) - 3, ".ms", 3) != 0 && strncmp(fileName + strlen(fileName) - 6, ".ms.gz", 6)) {
        printf("File does not have .ms or .ms.gz extension. Exiting!\n");
        return 1;
    }

    // Open the input file.
    gzFile file = gzopen(fileName, "r");
    // If file does not exist or is not compressed using gzip, return NULL.
    int errnum;
    gzerror(file, &errnum);
    if (errnum != Z_OK) {
        printf("File does not exist. Exiting!\n");
        return 1;
    }
    kstream_t* stream = ks_init(file);
    kstring_t* buffer = init_kstring(NULL);

    // Initalize memory used to read in a replicate.
    kvec_t(double) positions;
	kv_init(positions);
    kvec_t(kstring_t*) samples;
    kv_init(samples);

    // Eat lines until "segsites:" is encountered.
    do {
        ks_getuntil(stream, '\n', buffer, 0);
    } while (strncmp(ks_str(buffer), "segsites:", 9) != 0);

    int numReplicate = 0;

    // Parse all of the replicates.
    while (true) {

        int segsites = (int) strtol(ks_str(buffer) + 10, (char**) NULL, 10); 

        // Eat lines until "positions:" is encountered.
        do {
            ks_getuntil(stream, '\n', buffer, 0);
        } while (strncmp(ks_str(buffer), "positions:", 10) != 0);

        // Get the positions of the segsites.
        int numSpaces = 0, prevIndex;
        for (int i = 0; i <= ks_len(buffer); i++) {
            if (ks_str(buffer)[i] == ' ') {
                if (numSpaces > 0) {
                    double pos = strtod(ks_str(buffer) + prevIndex, (char**) NULL);
                    if (numSpaces > kv_size(positions)) {
                        kv_push(double, positions, pos); 
                    } else {
                        kv_A(positions, numSpaces - 1) = pos;
                    }
                }
                prevIndex = i;
                numSpaces++;
            }
        }

        // Now, read in all of the samples.
        int numSamples = 0;
        while (ks_getuntil(stream, '\n', buffer, 0) > 0 && strncmp(ks_str(buffer), "segsites:", 9) != 0) {
            if (numSamples >= kv_size(samples)) {
                kstring_t* temp = calloc(1, sizeof(kstring_t));
                kputs(ks_str(buffer), temp);
                kv_push(kstring_t*, samples, temp); 
            } else {
                ks_overwrite(ks_str(buffer), kv_A(samples, numSamples));
            }
            numSamples++;
        }

        // Convert the ms replicate to vcf.
        toVCF(fileName, length, unphased, missing, compress, numReplicate, segsites, numSamples, &positions, &samples);
        
        // If end of file, exit main loop.
        if (ks_eof(stream)) {
            break;
        }

        numReplicate++;

        // More replicates. Read until segsites is encountered.
        do {
            ks_getuntil(stream, '\n', buffer, 0);
        } while (strncmp(ks_str(buffer), "segsites:", 9) != 0);

    }

    // Free memory.
    ks_destroy(stream);
    free(buffer -> s);
    free(buffer);
    kv_destroy(positions);
    for (int i = 0; i < kv_size(samples); i++) {
        free(ks_str(kv_A(samples, i)));
    }
    kv_destroy(samples);
}