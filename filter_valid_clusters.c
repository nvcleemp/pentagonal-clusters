/*
 * Main developer: Nico Van Cleemput
 * 
 * Copyright (C) 2017 Ghent University.
 * Licensed under the GNU AGPL, read the file LICENSE for details.
 */

/* This program reads pentagonal adjacency graphs of fullerenes from standard in
 * and searches for specific clusters of pentagons.   
 * 
 * 
 * Compile with:
 *     
 *     cc -o filter_valid_clusters -O4 filter_valid_clusters.c
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
int numberOfValid = 0;

int nv;
int ne;

/* Below we list how many clusters of a given size are possible in a fullerene.
 * 
 * 1: 0 - 12
 * 2: 0 -  6
 * 3: 0 -  4
 * 4: 0 -  3
 * 5: 0 -  2
 */
int partitions[13][7][5][4][3];

//=============== Checking for property ===========================

int currentPartition[5];

boolean visited[12];
boolean currentCluster[12];
int currentClusterSize;

boolean validateCurrentCluster(){
    if(currentClusterSize < 3) return TRUE;
    if(currentClusterSize > 5) return FALSE;
    
    //build degree frequency table
    int degreeFreqTable[5];
    int i;
    for(i = 0; i < 5; i++) degreeFreqTable[i] = 0;
    for(i = 0; i < 12; i++){
        if(currentCluster[i]){
            degreeFreqTable[degree[i]]++;
        }
    }
    
    if(currentClusterSize == 3) return degreeFreqTable[2] == 3;
    if(currentClusterSize == 4) return degreeFreqTable[2] == 2 && degreeFreqTable[3] == 2;
    if(currentClusterSize == 5) return degreeFreqTable[2] == 2 
            && degreeFreqTable[3] == 2 
            && degreeFreqTable[4] == 1;
    return FALSE;
}

//some macros for the stack in the next method
#define PUSH(stack, value) stack[top++] = (value)
#define POP(stack) stack[--top]
#define STACKISEMPTY top==0
#define STACKISNOTEMPTY top>0

boolean hasValidClusters(){
    int i, j, top, currentVertex;
    int stack[12];
    
    for(i = 0; i < 12; i++){
        visited[i] = FALSE;
    }
    
    for(i = 0; i < 5; i++){
        currentPartition[i] = 0;
    }
    
    for(i = 0; i < nv; i++){
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
            
            //validate cluster
            if(!validateCurrentCluster()) return FALSE;
            
            currentPartition[currentClusterSize-1]++;
        }
    }
    
    partitions[currentPartition[0]][currentPartition[1]][currentPartition[2]][currentPartition[3]][currentPartition[4]]++;
    
    return TRUE;
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
    fprintf(stderr, "The program %s reads pentagonal adjacency graphs of\nfullerenes and searches them for specific clusters.\n\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "\nThis program can handle graphs up to %d vertices. Recompile if you need larger\n", MAXN);
    fprintf(stderr, "graphs.\n\n");
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "    -c, --count\n");
    fprintf(stderr, "       Print the number of times a partition appears.\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    /*=========== commandline parsing ===========*/
    
    boolean printCounts = FALSE;

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
        {"count", no_argument, NULL, 'c'},
        {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "hc", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                break;
            case 'c':
                printCounts = TRUE;
                break;
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
    
    int i, j, k, l, m;
    for(m=2; m>=0; m--){
        for(l=3; l>=0; l--){
            for(k=4; k>=0; k--){
                for(j=6; j>=0; j--){
                    for(i=0; i<13; i++){
                        partitions[i][j][k][l][m] = 0;
                    }
                }
            }
        }
    }

    /*=========== read pentagonal adjacency graphs ===========*/

    unsigned short code[MAXCODELENGTH];
    int length;
    while (readPlanarCode(code, &length, stdin)) {
        decodePlanarCode(code);
        if(nv!=12){
            fprintf(stderr, "This program only supports pentagonal adjacency graphs of fullerenes -- exiting!\n");
            return EXIT_FAILURE;
        }
        if(hasValidClusters()){
            numberOfValid++;
        }
        numberOfGraphs++;
    }
    
    for(m=2; m>=0; m--){
        for(l=3; l>=0; l--){
            for(k=4; k>=0; k--){
                for(j=6; j>=0; j--){
                    for(i=0; i<13; i++){
                        if(partitions[i][j][k][l][m]){
                            if(printCounts){
                                fprintf(stdout, "%d,%d,%d,%d,%d: %d\n", i, j, k, l, m, partitions[i][j][k][l][m]);
                            } else {
                                fprintf(stdout, "%d,%d,%d,%d,%d\n", i, j, k, l, m);
                            }
                        }
                    }
                }
            }
        }
    }
    
    fprintf(stderr, "Read %d graph%s.\n", numberOfGraphs, numberOfGraphs==1 ? "" : "s");
    fprintf(stderr, "Found %d valid cluster%s.\n", numberOfValid, 
                numberOfValid==1 ? "" : "s");
    
    return EXIT_SUCCESS;
}
