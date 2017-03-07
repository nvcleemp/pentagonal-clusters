/*
 * Main developer: Nico Van Cleemput
 * 
 * Copyright (C) 2017 Ghent University.
 * Licensed under the GNU AGPL, read the file LICENSE for details.
 */

/* This program reads pentagonal adjacency graphs of fullerenes from 
 * standard in and searches for realisations of a specific partition of 12.   
 * 
 * 
 * Compile with:
 *     
 *     cc -o appearances_of_clusters -O4 appearances_of_clusters.c
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>


#ifndef MAXN
#define MAXN 12            /* the maximum number of vertices */
#endif
#define MAXE (6*MAXN-12)    /* the maximum number of oriented edges */
#define MAXVAL (MAXN-1)  /* the maximum degree of a vertex */
#define MAXCODELENGTH (MAXN+MAXE+3)

#define FALSE 0
#define TRUE  1

typedef int boolean;

typedef struct e /* The data type used for edges */ {
    int start; /* vertex where the edge starts */
    int end; /* vertex where the edge ends */
    
    struct e *prev; /* previous edge in clockwise direction */
    struct e *next; /* next edge in clockwise direction */
    struct e *inverse; /* the edge that is inverse to this one */
} EDGE;

EDGE *firstedge[MAXN]; /* pointer to arbitrary edge out of vertex i. */
int degree[MAXN];

EDGE edges[MAXE];

int numberOfGraphs = 0;
int numberOfAppearancesOfRequestedPartition = 0;

int nv;
int ne;

int requestedPartition;

/* Below we list how many clusters of a given size are possible in a fullerene.
 * 
 *  1: 0 - 12
 *  2: 0 -  6
 *  3: 0 -  4
 *  4: 0 -  3
 *  5: 0 -  2
 *  6: 0 -  2
 *  7: 0 -  1
 *  8: 0 -  1
 *  9: 0 -  1
 * 10: 0 -  1
 * 11: 0 -  1
 * 12: 0 -  1
 */
int partitions2Number[13][7][5][4][3][3][2][2][2][2][2][2];

int possiblePartitions[77][12];
int possiblePartitionCount;
int currentPossiblePartition[12];

void buildPartitions_recursive(int currentSize){
    int i;
    int sum = 0;
    for(i=0; i<12; i++){
        sum+=(i+1)*currentPossiblePartition[i];
    }
    if(sum==12){
        for(i=0; i<12; i++){
            possiblePartitions[possiblePartitionCount][i] = currentPossiblePartition[i];
        }
        partitions2Number[currentPossiblePartition[0]]
                  [currentPossiblePartition[1]]
                  [currentPossiblePartition[2]]
                  [currentPossiblePartition[3]]
                  [currentPossiblePartition[4]]
                  [currentPossiblePartition[5]]
                  [currentPossiblePartition[6]]
                  [currentPossiblePartition[7]]
                  [currentPossiblePartition[8]]
                  [currentPossiblePartition[9]]
                  [currentPossiblePartition[10]]
                  [currentPossiblePartition[11]] = possiblePartitionCount;
        possiblePartitionCount++;
        return;
    }
    if(sum + currentSize <= 12){
        currentPossiblePartition[currentSize-1]++;
        buildPartitions_recursive(currentSize);
        currentPossiblePartition[currentSize-1]--;
    }
    if(currentSize > 1){
        buildPartitions_recursive(currentSize-1);
    }
}

void buildPartitions(){
    int i;
    for(i=0; i<12; i++){
        currentPossiblePartition[i]=0;
    }
    possiblePartitionCount = 0;
    buildPartitions_recursive(12);
}

void printPartition(FILE *f, int part, boolean newline){
    int i, j;
    
    for(i = 11; i >= 0; i--){
        for(j = 0; j < possiblePartitions[part][i]; j++){
            fprintf(f, "%d ", i+1);
        }
    }
    if(newline)
        fprintf(f, "\n");
}


//////////////////////////////////////////////////////////////////////////////

//=============== Checking for property ===========================

int currentPartition[12];

boolean visited[12];
boolean currentCluster[12];
int currentClusterSize;

//some macros for the stack in the next method
#define PUSH(stack, value) stack[top++] = (value)
#define POP(stack) stack[--top]
#define STACKISEMPTY top==0
#define STACKISNOTEMPTY top>0

void identifyClusters(){
    int i, j, top, currentVertex;
    int stack[12];
    
    for(i = 0; i < 12; i++){
        visited[i] = FALSE;
        currentPartition[i] = 0;
    }
    
    for(i = 0; i < 12; i++){
        if(!visited[i]){
            for(j = 0; j < 12; j++){
                currentCluster[j] = FALSE;
            }
            
            //build cluster containing vertex i
            top = 0;
            PUSH(stack, i);
            visited[i] = currentCluster[i] = TRUE;
            currentClusterSize = 1;
            while(STACKISNOTEMPTY){
                currentVertex = POP(stack);
                if(degree[currentVertex]){
                    EDGE *e, *elast;

                    e = elast = firstedge[currentVertex];
                    do {
                        if (!visited[e->end]) {
                            PUSH(stack, e->end);
                            visited[e->end] = currentCluster[e->end] = TRUE;
                            currentClusterSize++;
                        }
                        e = e->next;
                    } while (e != elast);
                }
            }
            
            currentPartition[currentClusterSize-1]++;
        }
    }
    
    int partitionNumber = partitions2Number[currentPartition[0]]
                  [currentPartition[1]]
                  [currentPartition[2]]
                  [currentPartition[3]]
                  [currentPartition[4]]
                  [currentPartition[5]]
                  [currentPartition[6]]
                  [currentPartition[7]]
                  [currentPartition[8]]
                  [currentPartition[9]]
                  [currentPartition[10]]
                  [currentPartition[11]];

    if(partitionNumber==requestedPartition){
        fprintf(stdout, "%d ", numberOfGraphs);
        numberOfAppearancesOfRequestedPartition++;
    }
    
}

//=============== Reading and decoding planarcode ===========================

EDGE *findEdge(int from, int to) {
    EDGE *e, *elast;

    e = elast = firstedge[from];
    do {
        if (e->end == to) {
            return e;
        }
        e = e->next;
    } while (e != elast);
    fprintf(stderr, "error while looking for edge from %d to %d.\n", from, to);
    exit(0);
}

void decodePlanarCode(unsigned short* code) {
    /* complexity of method to determine inverse isn't that good, but will have to satisfy for now
     */
    int i, j, codePosition;
    int edgeCounter = 0;
    EDGE *inverse;

    nv = code[0];
    codePosition = 1;

    for (i = 0; i < nv; i++) {
        degree[i] = 0;
        if(code[codePosition]){
            firstedge[i] = edges + edgeCounter;
            edges[edgeCounter].start = i;
            edges[edgeCounter].end = code[codePosition] - 1;
            edges[edgeCounter].next = edges + edgeCounter + 1;
            if (code[codePosition] - 1 < i) {
                inverse = findEdge(code[codePosition] - 1, i);
                edges[edgeCounter].inverse = inverse;
                inverse->inverse = edges + edgeCounter;
            } else {
                edges[edgeCounter].inverse = NULL;
            }
            edgeCounter++;
            codePosition++;
            for (j = 1; code[codePosition]; j++, codePosition++) {
                if (j == MAXVAL) {
                    fprintf(stderr, "MAXVAL too small: %d\n", MAXVAL);
                    exit(0);
                }
                edges[edgeCounter].start = i;
                edges[edgeCounter].end = code[codePosition] - 1;
                edges[edgeCounter].prev = edges + edgeCounter - 1;
                edges[edgeCounter].next = edges + edgeCounter + 1;
                if (code[codePosition] - 1 < i) {
                    inverse = findEdge(code[codePosition] - 1, i);
                    edges[edgeCounter].inverse = inverse;
                    inverse->inverse = edges + edgeCounter;
                } else {
                    edges[edgeCounter].inverse = NULL;
                }
                edgeCounter++;
            }
            firstedge[i]->prev = edges + edgeCounter - 1;
            edges[edgeCounter - 1].next = firstedge[i];
            degree[i] = j;
        }
        codePosition++; /* read the closing 0 */
    }

    ne = edgeCounter;
}

/**
 * 
 * @param code
 * @param length
 * @param file
 * @return returns 1 if a code was read and 0 otherwise. Exits in case of error.
 */
int readPlanarCode(unsigned short code[], int *length, FILE *file) {
    static int first = 1;
    unsigned char c;
    char testheader[20];
    int bufferSize, zeroCounter;
    
    int readCount;


    if (first) {
        first = 0;

        if (fread(&testheader, sizeof (unsigned char), 13, file) != 13) {
            fprintf(stderr, "can't read header ((1)file too small)-- exiting\n");
            exit(1);
        }
        testheader[13] = 0;
        if (strcmp(testheader, ">>planar_code") == 0) {

        } else {
            fprintf(stderr, "No planarcode header detected -- exiting!\n");
            exit(1);
        }
        //read reminder of header (either empty or le/be specification)
        if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
            return FALSE;
        }
        while (c!='<'){
            if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
                return FALSE;
            }
        }
        //read one more character
        if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
            return FALSE;
        }
    }

    /* possibly removing interior headers -- only done for planarcode */
    if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
        //nothing left in file
        return (0);
    }

    if (c == '>') {
        // could be a header, or maybe just a 62 (which is also possible for unsigned char
        code[0] = c;
        bufferSize = 1;
        zeroCounter = 0;
        code[1] = (unsigned short) getc(file);
        if (code[1] == 0) zeroCounter++;
        code[2] = (unsigned short) getc(file);
        if (code[2] == 0) zeroCounter++;
        bufferSize = 3;
        // 3 characters were read and stored in buffer
        if ((code[1] == '>') && (code[2] == 'p')) /*we are sure that we're dealing with a header*/ {
            while ((c = getc(file)) != '<');
            /* read 2 more characters: */
            c = getc(file);
            if (c != '<') {
                fprintf(stderr, "Problems with header -- single '<'\n");
                exit(1);
            }
            if (!fread(&c, sizeof (unsigned char), 1, file)) {
                //nothing left in file
                return (0);
            }
            bufferSize = 1;
            zeroCounter = 0;
        }
    } else {
        //no header present
        bufferSize = 1;
        zeroCounter = 0;
    }

    if (c != 0) /* unsigned chars would be sufficient */ {
        code[0] = c;
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant N too small %d > %d \n", code[0], MAXN);
            exit(1);
        }
        while (zeroCounter < code[0]) {
            code[bufferSize] = (unsigned short) getc(file);
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    } else {
        readCount = fread(code, sizeof (unsigned short), 1, file);
        if(!readCount){
            fprintf(stderr, "Unexpected EOF.\n");
            exit(1);
        }
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant N too small %d > %d \n", code[0], MAXN);
            exit(1);
        }
        bufferSize = 1;
        zeroCounter = 0;
        while (zeroCounter < code[0]) {
            readCount = fread(code + bufferSize, sizeof (unsigned short), 1, file);
            if(!readCount){
                fprintf(stderr, "Unexpected EOF.\n");
                exit(1);
            }
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    }

    *length = bufferSize;
    return (1);


}

//====================== USAGE =======================

void help(char *name) {
    fprintf(stderr, "The program %s searches pentagonal adjacency graphs of\nfullerenes for ones having cluster sizes that realise a specific partition\nof 12.\n\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "\nThis program can handle graphs up to %d vertices. Recompile if you need larger\n", MAXN);
    fprintf(stderr, "graphs.\n\n");
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "    -p, --partitions\n");
    fprintf(stderr, "       Print all possible partitions and then exits.\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    /*=========== commandline parsing ===========*/
    
    buildPartitions();

    int c, i;
    char *name = argv[0];
    static struct option long_options[] = {
        {"partitions", no_argument, NULL, 'p'},
        {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "hp", long_options, &option_index)) != -1) {
        switch (c) {
            case 'p':
                for(i = 0; i < possiblePartitionCount; i++){
                    fprintf(stderr, "partition %2d: ", i);
                    printPartition(stderr, i, TRUE);
                }
                return EXIT_SUCCESS;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case '?':
                usage(name);
                return EXIT_FAILURE;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }
    
    requestedPartition = atoi(argv[optind]);
    fprintf(stderr, "Looking for appearances of partition %d: ", requestedPartition);
    printPartition(stderr, requestedPartition, TRUE);

    /*=========== read planar graphs ===========*/

    unsigned short code[MAXCODELENGTH];
    int length;
    while (readPlanarCode(code, &length, stdin)) {
        decodePlanarCode(code);
        numberOfGraphs++;
        identifyClusters();
    }
    if(numberOfAppearancesOfRequestedPartition){
        fprintf(stdout, "\n");
    }
    
    fprintf(stderr, "Read %d graph%s.\n", numberOfGraphs, 
                numberOfGraphs==1 ? "" : "s");
    fprintf(stderr, "Found %d appearance%s of partition %d.\n", 
            numberOfAppearancesOfRequestedPartition, 
            numberOfAppearancesOfRequestedPartition==1 ? "" : "s",
            requestedPartition);
    
    return EXIT_SUCCESS;
}
