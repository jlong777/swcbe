/*
 sw_db.c
 
 farm out one fasta sequence per SPE from a db, and search for a specified
 query, allowing the user to specify the number of top hits to return.
 
 usage: sw_db [options] <scoring matrix file> <fasta file w/1 query> <db fasta file>
 
 example scoring matrix file (comments have a # in the first column, last
 column of letters is optional):
 
# PAM 1

#A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
 2  0 -6  0  0  0 -6  0  0  0  0  0  0  0  0  0  0  0  0 -6  0  0  0  0  0  0  A
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  B
-6  0  2  0  0  0 -6  0  0  0  0  0  0  0  0  0  0  0  0 -6  0  0  0  0  0  0  C
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  D
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  E
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  F
-6  0 -6  0  0  0  2  0  0  0  0  0  0  0  0  0  0  0  0 -6  0  0  0  0  0  0  G
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  H
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  I
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  J
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  K
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  L
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  M
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  N
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  O
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  P
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  Q
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  R
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  S
-6  0 -6  0  0  0 -6  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  T
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  U
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  V
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  W
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  Y
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  Z

 Author: James Long

 Copyright © 2008 University of Alaska Fairbanks. All Rights Reserved.

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this 
    list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, 
    this list of conditions and the following disclaimer in the documentation 
    and/or other materials provided with the distribution.

 3. Neither the name of the University of Alaska Fairbanks nor the names of 
    contributors to the software may be used to endorse or promote products 
    derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cbl.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* CBE includes */
#include <libspe2.h>
#include <pthread.h>
#include "control_block.h"

#define BUFSZ  2048
#define NUM_SPE   6
#define OUTLEN   65 /* length of header to output */
#define Q_LEN  2048 /* qLen <= (32767-g+e)/e, see spu code */
//#define TIMING

#ifdef TIMING
#include <sys/times.h>
#include <time.h>
#endif

typedef struct thread_args
{
  struct spe_context *spe;
  void *argp;
  void *envp;
} thread_args;

void spe_thread(thread_args *arg)
{
  unsigned int runflags=0;
  unsigned int entry=SPE_DEFAULT_ENTRY;

  if(spe_context_run(arg->spe, &entry, runflags, arg->argp, arg->envp, NULL) < 0)
  {
    fprintf(stderr, "sw_db: spe_context_run failure, exiting...\n");
    fflush(NULL);
    exit(1);
  }
  pthread_exit(NULL);
}

void PrintUsage(void)
{
  printf("align a query sequence against a database of fasta entries\n");
  printf("usage: sw_db [options] <scoring matrix file> <query fasta file> [db fasta file]\n");
  printf("              -a output alignments (default = top scores only)\n");
  printf("              -e extend gap penalty (eg), default = 1\n");
  printf("                 1 gap is scored as og + eg\n");
  printf("              -n number of top scores to report (default = 10)\n");
  printf("              -o open gap penalty (og), default = 8\n");
  printf("              -h help\n");
  printf("              <query fasta file> a single query\n");
  printf("              [db fasta file] may be specified or read from stdin\n\n");
  
  return;
}

int main(int argc, char *argv[])
{
  int i, j, k, alignments=0, lineNum, numSPEsWorking=0, numTopScores=10, noMore,
      option, wentToCompute;
  char alglenBuf[128], buf[BUFSZ], env[BUFSZ], header[NUM_SPE][BUFSZ], 
       **hiAlignments, **hiHeaders, path2spuBinary[BUFSZ], *p, *q, 
       record[BUFSZ], *tmpBuf,
       *dbBuf[NUM_SPE] __attribute__ ((aligned(16))), 
       qSeq[Q_LEN]     __attribute__ ((aligned(16)));
  long *algl, *algm, *algs, alglen, algstl, algsts, error=0;
  long dbLen[NUM_SPE], eg=1, eog, *hiAlignLen, hiScore, minHiScore, og=8, qLen,
       smax=0, tmpBufLen, tmpLen, tmpBufSz;
  long  sslookup_long[27][32];
  short sslookup[27][32] __attribute__ ((aligned(128)));
  void *memptr;
  FILE *dbFile, *fp, *matrix, *qFile;

  /* CBE declarations */
  int err, spe;
  void * speEA[NUM_SPE];
  unsigned int mboxData[1];
  pthread_t cellthrd[NUM_SPE];
  spe_context_ptr_t context[NUM_SPE];
  thread_args t_args[NUM_SPE];
  control_block cb[NUM_SPE] __attribute__ ((aligned (128)));
  spe_program_handle_t *spe_code;
  
#ifdef TIMING
  struct tms mytime; /* struct filled by calls to times() */
  float users, userf, syss, sysf, utot=0.0, stot=0.0; /* start and finish times */
  clock_t stime, ftime;
  long totsec=0, totsecCount=0;
  int stotOverflows=0, utotOverflows=0, totsecOverflows=0;
#endif

  /* options parsing */
  while((option = getopt(argc, argv, "e:n:o:ah")) > 0)
  {
    switch(option)
    {
      case 'a':
        alignments = 1;
        break;
	
      case 'e':
        for(i=0; i<strlen(optarg); i++)
          if(!isdigit(optarg[i]))
          {
            fprintf(stderr, "sw_db: the \"extend gap\" argument must be a number, returning...\n");
            return 1;
          }
        eg = atoi(optarg);
        if(eg < 0)
        {
          fprintf(stderr, "sw_db: the \"extend gap\" argument must be >= 0, returning...\n");
          return 1;
        }
        break;
        
      case 'n':
        numTopScores = (int) strtol(optarg, (char **)NULL, 10);
        break;
	
      case 'o':
        for(i=0; i<strlen(optarg); i++)
          if(!isdigit(optarg[i]))
          {
            fprintf(stderr, "sw_db: the \"open gap\" argument must be a number, returning...\n");
            return 1;
          }
        og = atoi(optarg);
        if(og < 0)
        {
          fprintf(stderr, "sw_db: the \"open gap\" argument must be >= 0, returning...\n");
          return 1;
        }
        break;
        
      case 'h':
        PrintUsage();
        break;
        
      default:
        printf("sw_db: unrecognized option -%c, use -h for help\n", option);
        break;
    }
  }
  
  eog = og + eg;
  
  if(argv[optind] == NULL)
  {
    fprintf(stderr, "sw_db: error, scoring matrix file missing, returning...\n\n");
    PrintUsage();
    return 1;
  }
  
  matrix = fopen(argv[optind], "r");
  if(matrix==NULL)
  {
    fprintf(stderr, "sw_db: error opening file %s, returning...\n", argv[optind]);
    return 1;
  }
  
  if(argv[optind+1] == NULL)
  {
    fprintf(stderr, "sw_db: error, query file missing, returning...\n\n");
    PrintUsage();
    return 1;
  }
  
  qFile = fopen(argv[optind+1], "r");
  if(qFile==NULL)
  {
    fprintf(stderr, "sw_db: error opening file %s, returning...\n", argv[optind+1]);
    return 1;
  }
  
  if(argv[optind+2] == NULL)
    dbFile = stdin;
  else
    dbFile = fopen(argv[optind+2], "r");
    
  if(dbFile==NULL)
  {
    fprintf(stderr, "sw_db: error opening file %s, returning...\n", argv[optind+2]);
    return 1;
  }
  
  /* sw_db_spu is in same directory as sw_db. Where is that? */
  fp = NULL;
  strcpy(path2spuBinary, argv[0]);
  if((p = strrchr(path2spuBinary, '/')) != NULL) /* was a path given? */
  {
    *p = '\0';
    strcat(path2spuBinary, "/sw_db_spu");
  }
  else /* binary is in $PATH */
  {
    strcpy(env, getenv("PATH"));
    p = strtok(env, ":");
    if(p!=NULL)
    {
      strcpy(path2spuBinary, p);
      strcat(path2spuBinary, "/sw_db_spu");
    
      while((fp = fopen(path2spuBinary, "rb")) == NULL)
      {
        p = strtok(NULL, ":");
        if(p==NULL) break;

        strcpy(path2spuBinary, p);
        strcat(path2spuBinary, "/sw_db_spu");
      }
    }
    
    if(fp!=NULL)
      fclose(fp);
    else
    {
      fprintf(stderr, "sw_db: error, sw_db_spu not found, exiting...\n");
      return 1;
    }
  }
  
  /* init */
  for(i=0; i<NUM_SPE; i++)
    dbBuf[i] = NULL;
    
  if(alignments)
  {
    hiAlignLen   = (long *)  malloc(numTopScores*sizeof(long *));
    hiAlignments = (char **) malloc(numTopScores*sizeof(char *));
    if(hiAlignLen==NULL || hiAlignments==NULL)
    {
      fprintf(stderr, "sw_db: malloc error, returning...\n");
      return 1;
    }
    
    for(i=0; i<numTopScores; i++)
      hiAlignments[i] = NULL;
  }
  
  hiHeaders = (char **) malloc(numTopScores*sizeof(char *));
  if(hiHeaders==NULL)
  {
    fprintf(stderr, "sw_db: malloc error, returning...\n");
    return 1;
  }
  
  for(i=0; i<numTopScores; i++)
  {
    hiHeaders[i] = (char *) malloc(BUFSZ);
    if(hiHeaders[i]==NULL)
    {
      fprintf(stderr, "sw_db: malloc error, returning...\n");
      return 1;
    }
    
    strcpy(hiHeaders[i], "0 ");
  }

  /* get the scoring matrix */
  for(j=0; j<32; j++)
  {
    sslookup[0][j]      = -1;
    sslookup_long[0][j] = -1;
  }
  
  i = 1; lineNum = 0;
  while(fgets(buf, BUFSZ, matrix) && i<27)
  {
    lineNum++;
    if(buf[0]=='#' || buf[0]=='\n' || strlen(buf)<51) continue;

    sslookup[i][0]      = -1;
    sslookup_long[i][0] = -1;
    j = 1;
    p = buf;
    while(j<27)
    {
      sslookup[i][j] = (short) strtol(p, &q, 10);
      sslookup_long[i][j] = (long) sslookup[i][j];

      if(j)
      {
        if(p==q)
        {
          if(p==buf)
            fprintf(stderr, "sw_db: no digits on line %d of matrix file \"%s\", returning...\n", lineNum, argv[optind]);
          else
            fprintf(stderr, "sw_db: illegal scoring matrix value, returning...\n");
          return 1;
        }
        p = q;
      }
      
      j++;
    }
    
    while(j<32)
    {
      sslookup_long[i][j] = -1;
      sslookup[i][j++]    = -1;
    }
    
    i++;
  }
  fclose(matrix);
  
  /* get the query */
  qLen = 0;
  qSeq[0] = '\0';
  while(fgets(buf, BUFSZ, qFile))
  {
    if(buf[0] == '>')
    {
      printf("query:\n%s", buf);
      
      while(fgets(buf, BUFSZ, qFile))
      {
        if(buf[0]!='>' && buf[0]!='\n' && buf[0]!='\t' && buf[0]!=' ')
        {
          for(i=0; i<strlen(buf); i++)
            if(buf[i] == '\n') {buf[i] = '\0'; break;}
          buf[BUFSZ-1] = '\0';
          
          qLen += strlen(buf);

          printf("%s\n", buf);
          strcat(qSeq, buf);
        }
      }
      if(alignments)
        printf("\ntop %d alignments:\n\n", numTopScores);
      else
        printf("\ntop %d scores:\n\n", numTopScores);
    }
  }
  fclose(qFile);
  
  /* launch spe threads */
  spe_code = spe_image_open(path2spuBinary);
  
  /* create the contexts */
  for(i=0; i<NUM_SPE; i++)
  {
    context[i] = spe_context_create(0, NULL);
    speEA[i]   = spe_ls_area_get(context[i]);
  }
  
  for(i=0; i<NUM_SPE; i++)
  {
    spe_program_load(context[i], spe_code);
  
    cb[i].addr_matrix = (unsigned int) &sslookup[0][0];
    cb[i].addr_qSeq   = (unsigned int) qSeq;
    cb[i].qLen        = (unsigned int) qLen;
    cb[i].eg          = (unsigned int) eg;
    cb[i].eog         = (unsigned int) eog;
    cb[i].speId       = (unsigned int) i;
  
    t_args[i].spe  = context[i];
    t_args[i].argp = &cb[i];
    t_args[i].envp = NULL;
    
    /* spawn a thread for each spe */
    err = pthread_create(&cellthrd[i], NULL, (void *(*)(void *))spe_thread, (void *)(t_args+i));
    if(err)
    {
      fprintf(stderr, "sw_db: pthread creation error, returning...\n");
      return 1;
    }
  }

  /* grab seqs */
  noMore = spe = wentToCompute = 0;
start:
  while(fgets(buf, BUFSZ, dbFile))
  {
    if(buf[0] == '>')
    {
next_header:
      tmpBuf = (char *) malloc(BUFSZ); /* buffer to hold the sequence */
      if(tmpBuf==NULL)
      {
        fprintf(stderr, "sw_db: malloc error, exiting...\n");
        return 1;
      }
      tmpBuf[0] = '\0';
      tmpLen = 0;
      tmpBufSz = BUFSZ;
      strcpy(record, buf);

      while(fgets(buf, BUFSZ, dbFile))
      {
        if(buf[0]!='>' && buf[0]!='\n' && buf[0]!='\t' && buf[0]!=' ')
        {
          for(i=0; i<BUFSZ; i++)
            if(buf[i] == '\n') {buf[i] = '\0'; break;}
        
          tmpBufLen = tmpLen;
          tmpLen += i;
        
          if(tmpLen >= tmpBufSz-128) /* tmpBuf have enough room? */
          {
            tmpBufSz += BUFSZ + tmpBufSz/10;
            tmpBuf = (char *) realloc(tmpBuf, tmpBufSz);
            if(tmpBuf == NULL)
            {
              fprintf(stderr, "sw_db: realloc error, exiting...\n");
              fflush(NULL);
              return (1);
            }
          }
          memcpy(tmpBuf+tmpBufLen, buf, i);
          tmpBuf[tmpBufLen+i] = '\0';
          
        }
        else /* at end of seq, possibly at next header */
        {
compute:
          if(noMore) goto finish;
          
          posix_memalign(&memptr, 16, ((tmpBufSz+15)/16)*16);
          memcpy(memptr, tmpBuf, tmpBufSz);
          free(tmpBuf);
          
#ifdef TIMING
          stime = clock()/CLOCKS_PER_SEC;
          times(&mytime);
          users = (float) mytime.tms_utime;
          syss  = (float) mytime.tms_stime;
#endif
          /*
           loop thru spes; if an spe has data, record, free memory, and hand
           off new data. if the spe had a high score, save the fasta header
          */
          while(1) /* loop until data handed out */
          {
            if(spe_out_mbox_status(context[spe]))
            {
              spe_out_mbox_read(context[spe], mboxData, 1);

              if((int)mboxData[0] >= -1) /* -1 => ready, but no data processed yet */
              {
                if((int)mboxData[0] >= 0) numSPEsWorking--;
//printf("%d\n", (int)mboxData[0]); /* returned smax */
                minHiScore = strtol(hiHeaders[numTopScores-1], (char **)NULL, 10);
                if((long)mboxData[0] > minHiScore)
                {
                  /*
                   save header this spe worked on, and if the -a option
                   was given, save the sequence for potential alignment
                  */
                  for(i=0; i<numTopScores; i++)
                  {
                    hiScore = strtol(hiHeaders[i], (char **)NULL, 10);
                    if((int)mboxData[0] > hiScore) /* bump everybody down */
                    {
                      if(alignments && hiAlignments[numTopScores-1]!=NULL)
                        free(hiAlignments[numTopScores-1]);
                      
                      for(j=numTopScores-1; j>i; j--)
                      {
                        strcpy(hiHeaders[j], hiHeaders[j-1]);
                        if(alignments)
                        {
                          hiAlignments[j] = hiAlignments[j-1];
                          hiAlignLen[j]   = hiAlignLen[j-1];
                        }
                      }
                      
                      sprintf(hiHeaders[i], "%ld\t", (long)mboxData[0]);

                      for(k=0; k<OUTLEN; k++)
                      {
                        if(header[spe][k]=='\n')
                        {
                          header[spe][k] = '\0';
                          break;
                        }
                      }

                      header[spe][OUTLEN] = '\0';
                      strcat(hiHeaders[i], header[spe]+1);
                      
                      if(alignments)
                      {
                        /* save the sequence to compute alignment later */
                        hiAlignments[i] = dbBuf[spe];
                        hiAlignLen[i]   = dbLen[spe];
                        dbBuf[spe] = NULL;
                      }
                      
                      break;
                    }
                  }
                }
                
                /* prepare data for spe */
                memcpy(header[spe], record, BUFSZ);
//printf("record = %s\n", record);
                if(dbBuf[spe] != NULL) free(dbBuf[spe]);
                
                dbBuf[spe] = memptr;
                dbLen[spe] = tmpLen;
                
                /* signal spe that data is ready */
                mboxData[0] = dbLen[spe];
                spe_in_mbox_write(context[spe], mboxData, 1, SPE_MBOX_ANY_NONBLOCKING);
              
                mboxData[0] = (unsigned int)dbBuf[spe];
                spe_in_mbox_write(context[spe], mboxData, 1, SPE_MBOX_ANY_NONBLOCKING);
                
                numSPEsWorking++;

                if(wentToCompute) noMore = 1;
                break;
              }
            }
            
            spe++;
            if(spe==NUM_SPE) spe = 0;
          }
          
#ifdef TIMING
          times(&mytime);
          userf = (float) mytime.tms_utime;
          sysf  = (float) mytime.tms_stime;
          ftime = clock()/CLOCKS_PER_SEC;

          if((userf-users)>=0) utot += userf-users;
          else utotOverflows++;
          if((sysf-syss)>=0) stot += sysf-syss;
          else stotOverflows++;
          if((ftime-stime)>=0) {totsec += ftime-stime; totsecCount++;}
          else totsecOverflows++;
//if(ftime-stime>0) printf("%ld\n", ftime-stime);
#endif
          if(buf[0]=='>')
            goto next_header;
          else if(buf[0]=='\n' || buf[0]=='\t' || buf[0]==' ')
            goto start;
        }
      }
      /* no more input, so process last sequence */
      wentToCompute = 1;
      goto compute;
    }
  }

finish:

#ifdef TIMING
  stime = clock()/CLOCKS_PER_SEC;
  times(&mytime);
  users = (float) mytime.tms_utime;
  syss  = (float) mytime.tms_stime;
#endif

  /* collect remaining results */
  spe = 0;
  while(numSPEsWorking)
  {
    if(spe_out_mbox_status(context[spe]))
    {
      spe_out_mbox_read(context[spe], mboxData, 1);
      if((int)mboxData[0] >= 0) numSPEsWorking--;
//printf("%d\n", (int)mboxData[0]); /* returned smax */
      minHiScore = strtol(hiHeaders[numTopScores-1], (char **)NULL, 10);
      if((long)mboxData[0] > minHiScore)
      {
        /*
         save header this spe worked on, and if the -a option
         was given, save the sequence for potential alignment
        */
        for(i=0; i<numTopScores; i++)
        {
          hiScore = strtol(hiHeaders[i], (char **)NULL, 10);
          if((int)mboxData[0] > hiScore) /* bump everybody down */
          {
            if(alignments && hiAlignments[numTopScores-1]!=NULL)
              free(hiAlignments[numTopScores-1]);
            
            for(j=numTopScores-1; j>i; j--)
            {
              strcpy(hiHeaders[j], hiHeaders[j-1]);
              if(alignments)
              {
                hiAlignments[j] = hiAlignments[j-1];
                hiAlignLen[j]   = hiAlignLen[j-1];
              }
            }
            
            sprintf(hiHeaders[i], "%ld\t", (long)mboxData[0]);

            for(k=0; k<OUTLEN; k++)
            {
              if(header[spe][k]=='\n')
              {
                header[spe][k] = '\0';
                break;
              }
            }

            header[spe][OUTLEN] = '\0';
            strcat(hiHeaders[i], header[spe]+1);
            
            if(alignments)
            {
              /* save the sequence to compute alignment later */
              hiAlignments[i] = dbBuf[spe];
              hiAlignLen[i]   = dbLen[spe];
              dbBuf[spe] = NULL;
            }
            
            break;
          }
        }
      }
    }
    spe++;
    if(spe==NUM_SPE) spe = 0;
  }
  
#ifdef TIMING
  times(&mytime);
  userf = (float) mytime.tms_utime;
  sysf  = (float) mytime.tms_stime;
  ftime = clock()/CLOCKS_PER_SEC;
          
  if((userf-users)>0) utot += userf-users;
  else utotOverflows++;
  if((sysf-syss)>0) stot += sysf-syss;
  else stotOverflows++;
  if((ftime-stime)>0) {totsec += ftime-stime; totsecCount++;}
  else totsecOverflows++;
#endif
  
  if(alignments) /* compute alignments of top scores */
  {
    for(i=0; i<numTopScores; i++)
    {
      if(hiAlignments[i] != NULL)
      {
        cb_swa_fw((long *)qSeq, qLen, (long *)hiAlignments[i], hiAlignLen[i],
                  eg, og, &sslookup_long[0][0], &smax, &algl, &algm, &algs, 
                  &alglen, &algstl, &algsts, &error);
      
        free(hiAlignments[i]);
        
        hiAlignments[i] = (char *) malloc(3*alglen + 32);
        if(hiAlignments[i]==NULL)
        {
          fprintf(stderr, "sw_db: malloc error, returning...\n");
          return 1;
        }
      
        sprintf(alglenBuf, "%ld\n", algsts);
        strcpy (hiAlignments[i], alglenBuf);
        strncat(hiAlignments[i], (char *) algs, alglen); /* db */
        strcat (hiAlignments[i], "\n");
        strncat(hiAlignments[i], (char *) algm, alglen);
        strcat (hiAlignments[i], "\n");
        strncat(hiAlignments[i], (char *) algl, alglen); /* query */
        strcat (hiAlignments[i], "\n");
        sprintf(alglenBuf, "%ld\n", algstl);
        strcat (hiAlignments[i], alglenBuf);
        
        cb_free(algl);
        cb_free(algm);
        cb_free(algs);
      }
    }
  }
  
  /* now output top scores */
  for(i=0; i<numTopScores; i++)
  {
    printf("%s\n", hiHeaders[i]);
    if(alignments)
      printf("%s\n", hiAlignments[i]);
  }

#ifdef TIMING
  printf("\ntotal clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
  printf("total approx secs not in reading db = %d\n", totsec);
  printf("utotOverflows = %d\n", utotOverflows);
  printf("stotOverflows = %d\n", stotOverflows);
  printf("totsecOverflows = %d\n", totsecOverflows);
  printf("totsecCount = %d\n", totsecCount);
#endif
  
  /* cleanup */ 
  mboxData[0] = -1;
  for(i=0; i<NUM_SPE; i++)
  {
    spe_in_mbox_write(context[i], mboxData, 1, SPE_MBOX_ANY_NONBLOCKING);
    pthread_join(cellthrd[i], NULL);
  }
   
  for(i=0; i<NUM_SPE; i++)
    spe_context_destroy(context[i]);
    
  spe_image_close(spe_code);
      
  return 0;
}
