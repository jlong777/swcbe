/*
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

#include <malloc.h>
#include <spu_mfcio.h>
#include <stdio.h>
#include <stdlib.h>
#include "control_block.h"

//#define DEBUG_INFO
#define DIAG_LEN     8 /* HARD CODED number of s-w cells along minor diagonal */
#define MAX_CHUNK 1024 /* <= (LS size) - (spe image & data size), & %128 = 0 */
//#define TIMING
/*
 spu_clock_reads are in TimeBase units, divide by below to get seconds:
 
 Decrement register frequencies
 ==============================
 QS20        14318000
 QS/21/22    26666666
 PS3         79800000 
 
 set this value below in #define FREQ
*/
//#define TRACE
#define VECT_LEN     8 /* how many slots used in a vector, i.e. short=8 */
#define VECTOR

#ifdef TIMING
#define FREQ 79800000
#include <limits.h>
#include "/opt/cell/sdk/prototype/usr/spu/include/spu_timer.h"
#endif

#ifdef TRACE
#include <trace_user.h>
#endif

volatile control_block cb __attribute__ ((aligned (128)));

int main(unsigned long long spuid __attribute__ ((unused)),
         unsigned long long argp  __attribute__ ((unused)),
         unsigned long long envp  __attribute__ ((unused)))
{
  int i, j, jj, k, kk, m, n, chunk, dbLen, err, gotLast, lastI, qLen, savj;
  volatile char q[DIAG_LEN]      __attribute__ ((aligned(128))),
                db[2][MAX_CHUNK] __attribute__ ((aligned(128)));
  int gv[MAX_CHUNK] __attribute__ ((aligned(128)));
  short e, g, *gapH, *max2, smax, *temp[2], *query, last[MAX_CHUNK],
        gh[MAX_CHUNK]  __attribute__ ((aligned(128))), 
        max[DIAG_LEN]  __attribute__ ((aligned(128))), 
        max1[DIAG_LEN],
        max3[DIAG_LEN]   __attribute__ ((aligned(128))), 
        sslookup[27][32] __attribute__ ((aligned(128))),
        swTAB[DIAG_LEN+1][MAX_CHUNK+128] __attribute__ ((aligned(128))),
        zero=0;
  unsigned int dbAdd __attribute__ ((aligned(16)));
  
#ifdef TIMING
  uint64_t clock=0, start=0, totTime=0;
#endif

#ifdef TRACE
  trace_payload_t pl;
  trace_interval_p interval;
#endif

#ifdef VECTOR
  vector signed int vei, *vv_0, *vv_1, vgv_0, vgv_1, 
                    vmax_0, vmax_1, vmax3_0, vmax3_1;
  vector signed short  vgh, *vmax, *vmax1, *vmax2, *vmax3, vsmax, 
                       vdb, vss, vsw, *vh, ve, vg;
  vector unsigned int vu_0, vu_1;
  vector unsigned short vu;
  unsigned char flipPattern[16] __attribute__ ((aligned(128)))
                               ={0x0E,0x0F,  0x0C,0x0D,  0x0A,0x0B,  0x08,0x09,
                                 0x06,0x07,  0x04,0x05,  0x02,0x03,  0x00,0x01};
  unsigned char selectPattern[16] __attribute__ ((aligned(128)))
                               ={0x80,0x80,0x0E,0x0F,  0x80,0x80,0x0C,0x0D,
                                 0x80,0x80,0x0A,0x0B,  0x80,0x80,0x08,0x09};
  unsigned char selectPattern2[16] __attribute__ ((aligned(128)))
                               ={0x80,0x80,0x06,0x07,  0x80,0x80,0x04,0x05,
                                 0x80,0x80,0x02,0x03,  0x80,0x80,0x00,0x01};
  unsigned char selectPattern3[16] __attribute__ ((aligned(128)))
                               ={0x1E,0x1F,  0x1A,0x1B,  0x16,0x17,  0x12,0x13,
                                 0x0E,0x0F,  0x0A,0x0B,  0x06,0x07,  0x02,0x03};
  unsigned char selectPattern4[16] __attribute__ ((aligned(128)))
                               ={0x04,0x05,0x06,0x07,  0x08,0x09,0x0A,0x0B,
                                 0x0C,0x0D,0x0E,0x0F,  0x10,0x11,0x12,0x13};
  vector unsigned char *flip, *select, *select2, *select3, *select4;
#endif

  if(MAX_CHUNK%128)
  {
    fprintf(stderr, "sw_db_spu: MAX_CHUNK = %d is not 0%%128, exiting...\n", MAX_CHUNK);
    fflush(NULL);
    err = 1;
  }
  
  mfc_get(&cb, argp, sizeof(cb), 0, 0, 0);
  mfc_write_tag_mask(0x01);
  mfc_read_tag_status_all();
  
  /* let ppu know we are ready */
  spu_write_out_mbox(-1); /* -1 => ready, but no data processed yet */
  
#ifdef VECTOR
  vmax  = (vector signed short *) max;
  vmax1 = (vector signed short *) max1;
  vmax3 = (vector signed short *) max3;
  flip    = (vector unsigned char *) flipPattern;
  select  = (vector unsigned char *) selectPattern;
  select2 = (vector unsigned char *) selectPattern2;
  select3 = (vector unsigned char *) selectPattern3;
  select4 = (vector unsigned char *) selectPattern4;
#endif

  e    = -(signed short) cb.eg;
  g    = -(signed short) cb.eog;
  qLen = cb.qLen;
  
#ifdef VECTOR
  ve  = spu_splats(e);
  vei = spu_splats((int)e);
  vg  = spu_splats(g);
#endif
  e = -e;
  g = -g;

  max2    = (signed short *) memalign(128, ((qLen+16)/16)*16*sizeof(signed short));
  gapH    = (signed short *) memalign(16,  ((qLen+16)/16)*16*sizeof(signed short));
  temp[0] = (signed short *) memalign(16,  ((qLen+16)/16)*16*sizeof(signed short));
  temp[1] = (signed short *) memalign(16,  ((qLen+16)/16)*16*sizeof(signed short));
  query   = (signed short *) memalign(16,  ((qLen+16)/16)*16*sizeof(signed short));
  if(max2==NULL || gapH==NULL || temp[0]==NULL || temp[1]==NULL || query==NULL)
  {
    fprintf(stderr, "sw_db_spu: malloc error, exiting...\n");
    fflush(NULL);
    err = 1;
  }

  mfc_get((void *)query, cb.addr_qSeq, ((qLen+16)/16)*16, 0, 0, 0);
  mfc_get((void *)sslookup, cb.addr_matrix, 27*32*sizeof(short), 0, 0, 0);

  /*
   the query size is limited, because -g - (qLen-1)*e can
   never be less than -32767, hence qLen <= (32767-g+e)/e
  */
  if((-((int)g) - (qLen-1)*(int)e) < -32767)
  {
    fprintf(stderr, "sw_db_spu: query length must be <= (32767-g+e)/e, exiting...\n");
    fflush(NULL);
    err = 1;
  }

#ifdef TIMING
  spu_clock_start();
#endif

  while(1)
  {
    /* init */
    for(i=1; i<qLen+1; i++)
    {
      max2[i-1] = 0 - g;
      gapH[i-1] = -g - (i-1)*e;
      temp[0][i-1] = 0;
      temp[1][i-1] = 0;
    }
    temp[0][i-1] = 0;
    
    err = smax = 0;
    
#ifdef VECTOR
    vsmax = spu_splats(zero);
#endif

    /* blocking command */
    dbLen = spu_read_in_mbox(); /* db sequence length */
    if(dbLen < 0) 
      goto finish;

    dbAdd = spu_read_in_mbox(); /* a sequence is ready to fetch */

#ifdef TIMING
    start = spu_clock_read();
    while(start > UINT_MAX) start = spu_clock_read(); /* sometimes get a bad value */
#endif
    
    /* compute chunk size */
    chunk = dbLen;
    while(chunk > MAX_CHUNK) chunk /= 2;
    chunk = (chunk/16)*16;
    if(chunk < DIAG_LEN) chunk = 0;

    gotLast = lastI = 0;
    if(chunk)
    {
      savj = (dbLen/chunk)*chunk + 1;
      
      /* get the first chunk of db */
      mfc_get((void *)db[0], dbAdd, chunk, 1, 0, 0);
      mfc_write_tag_mask(0x10);
      mfc_read_tag_status_all();
      
      for(jj=0, n=0; jj<dbLen+1-chunk; jj+=chunk, n++) /* strips */
      {
        for(k=0; k<chunk+1; k++)
          swTAB[DIAG_LEN][k+3] = 0; /* copied to row 0 below, which is 0 initially */
      
        /* get the next chunk (double buffering) */
        if(jj < dbLen+1-2*chunk)
        {
          mfc_get((void *)db[(n+1)%2], dbAdd+jj+chunk, chunk, 2, 0, 0);
        }
        else /* for last strip */
        {
          gotLast = 1;
          mfc_get((void *)db[(n+1)%2], dbAdd+savj-1, ((dbLen+17-savj)/16)*16, 2, 0, 0);
        }

        for(k=0; k<chunk; k++)
          gv[k] = -g - (jj+k)*e;
        
        for(i=1; i<qLen+1-DIAG_LEN; i+=DIAG_LEN)
        {
          j = jj + 1;
#ifdef VECTOR
          vdb = spu_insert((signed int)db[n%2][7], vdb, 0);
          vdb = spu_insert((signed int)db[n%2][6], vdb, 1);
          vdb = spu_insert((signed int)db[n%2][5], vdb, 2);
          vdb = spu_insert((signed int)db[n%2][4], vdb, 3);
          vdb = spu_insert((signed int)db[n%2][3], vdb, 4);
          vdb = spu_insert((signed int)db[n%2][2], vdb, 5);
          vdb = spu_insert((signed int)db[n%2][1], vdb, 6);
          vdb = spu_insert((signed int)db[n%2][0], vdb, 7);
#endif
          for(k=0; k<chunk+1; k++) /* copy bottom row to top */
            swTAB[0][k+3] = swTAB[DIAG_LEN][k+3];
          
          for(k=0; k<DIAG_LEN; k++)
          {
            q[k] = 0x1F & ((char *)query)[i+k-1];
            swTAB[k][0+3] = temp[n%2][i+k-1];
          }

          /* fill upper left corner */
          for(k=0; k<DIAG_LEN-1; k++) /* i offset */
          {
            for(kk=0; kk<DIAG_LEN-1-k; kk++) /* j offset */
            {
              max[k]  = 0x0;
              max1[k] = swTAB[k][kk+3] + sslookup[q[k]][0x1F & db[n%2][kk]];
              max3[k] = swTAB[k][kk+1+3] - g;
      
              if(kk)
                gh[kk] = (max2[i+k-1] > gh[kk-1]-e)   ? max2[i+k-1]: gh[kk-1]-e;
              else
                gh[0]  = (max2[i+k-1] > gapH[i+k-1]-e)? max2[i+k-1]: gapH[i+k-1]-e;
            
              gv[kk] = (max3[k] > gv[kk]-e)? max3[k]: gv[kk]-e;
            
              if(max1[k]>max[k]) max[k] = max1[k];
              if(gh[kk]>max[k])  max[k] = gh[kk];
              if(gv[kk]>max[k])  max[k] = gv[kk];
              if(max[k]>smax)      smax = max[k];
#ifdef DEBUG_INFO
              if(max[k] & 0x8000) printf("sw_db_spu1: cell score exceeds maximum...\n");
#endif
              swTAB[k+1][kk+1+3] = max[k];
              max2[i+k-1] = max[k] - g;
            }
          }

          /* fill the middle, use vector ops and loop unrolling */
#ifndef VECTOR
          for(k=0; k<DIAG_LEN; k++) /* i & j offset */
/* a */     max[k] = 0x0;
            
          for(k=0; k<DIAG_LEN; k++)
/* b */     max1[k] = swTAB[k][DIAG_LEN-k-1+3] + sslookup[0x1F & db[n%2][DIAG_LEN-k-1]][q[k]];
          
          for(k=0; k<DIAG_LEN; k++)
/* c */     max3[k] = swTAB[k][DIAG_LEN-k+3] - g;
            
          for(k=0; k<DIAG_LEN; k++)
            gh[DIAG_LEN-k-1] = (max2[i+k-1] > gh[DIAG_LEN-k-2]-e)  ? max2[i+k-1]: gh[DIAG_LEN-k-2]-e;
/* d */   gh[0] = (max2[i+DIAG_LEN-2] > gapH[i+DIAG_LEN-2]-e) ? max2[i+DIAG_LEN-2]: gapH[i+DIAG_LEN-2]-e;
          
          for(k=0; k<DIAG_LEN; k++)
/* e */     gv[DIAG_LEN-k-1] = (max3[k] > gv[DIAG_LEN-k-1]-e)? max3[k]: gv[DIAG_LEN-k-1]-e;
          
          for(k=0; k<DIAG_LEN; k++)
/* f */     if(max1[k]>max[k])   max[k] = max1[k];
          
          for(k=0; k<DIAG_LEN; k++)
/* g */     if(gh[DIAG_LEN-k-1]>max[k]) max[k] = gh[DIAG_LEN-k-1];
            
          for(k=0; k<DIAG_LEN; k++)
/* h */     if(gv[DIAG_LEN-k-1]>max[k]) max[k] = gv[DIAG_LEN-k-1];
            
          for(k=0; k<DIAG_LEN; k++)
/* i */     if(max[k]>smax) smax = max[k];
#else
          vmax2 = (vector signed short *) (max2+i-1);
          
/* a */   *vmax = spu_splats(zero);
            
/* b */   max1[0] = swTAB[0][10] + sslookup[0x1F & db[n%2][7]][q[0]];
          max1[1] = swTAB[1][9]  + sslookup[0x1F & db[n%2][6]][q[1]];
          max1[2] = swTAB[2][8]  + sslookup[0x1F & db[n%2][5]][q[2]];
          max1[3] = swTAB[3][7]  + sslookup[0x1F & db[n%2][4]][q[3]];
          max1[4] = swTAB[4][6]  + sslookup[0x1F & db[n%2][3]][q[4]];
          max1[5] = swTAB[5][5]  + sslookup[0x1F & db[n%2][2]][q[5]];
          max1[6] = swTAB[6][4]  + sslookup[0x1F & db[n%2][1]][q[6]];
          max1[7] = swTAB[7][3]  + sslookup[0x1F & db[n%2][0]][q[7]];
        
/* c */   max3[0] = swTAB[0][11] - g;
          max3[1] = swTAB[1][10] - g;
          max3[2] = swTAB[2][9]  - g;
          max3[3] = swTAB[3][8]  - g;
          max3[4] = swTAB[4][7]  - g;
          max3[5] = swTAB[5][6]  - g;
          max3[6] = swTAB[6][5]  - g;
          max3[7] = swTAB[7][4]  - g;
          
/* d */   gh[7] = (max2[i-1] > gh[6]-e)    ? max2[i-1]: gh[6]-e;
          gh[6] = (max2[i  ] > gh[5]-e)    ? max2[i  ]: gh[5]-e;
          gh[5] = (max2[i+1] > gh[4]-e)    ? max2[i+1]: gh[4]-e;
          gh[4] = (max2[i+2] > gh[3]-e)    ? max2[i+2]: gh[3]-e;
          gh[3] = (max2[i+3] > gh[2]-e)    ? max2[i+3]: gh[2]-e;
          gh[2] = (max2[i+4] > gh[1]-e)    ? max2[i+4]: gh[1]-e;
          gh[1] = (max2[i+5] > gh[0]-e)    ? max2[i+5]: gh[0]-e;
          gh[0] = (max2[i+6] > gapH[i+6]-e)? max2[i+6]: gapH[i+2]-e;
        
/* e */   gv[7] = (max3[0]   > gv[7]-e)    ? max3[0]  : gv[7]-e;
          gv[6] = (max3[1]   > gv[6]-e)    ? max3[1]  : gv[6]-e;
          gv[5] = (max3[2]   > gv[5]-e)    ? max3[2]  : gv[5]-e;
          gv[4] = (max3[3]   > gv[4]-e)    ? max3[3]  : gv[4]-e;
          gv[3] = (max3[4]   > gv[3]-e)    ? max3[4]  : gv[3]-e;
          gv[2] = (max3[5]   > gv[2]-e)    ? max3[5]  : gv[2]-e;
          gv[1] = (max3[6]   > gv[1]-e)    ? max3[6]  : gv[1]-e;
          gv[0] = (max3[7]   > gv[0]-e)    ? max3[7]  : gv[0]-e;
     
          vu    = spu_cmpgt(*vmax1, *vmax);
/* f */   *vmax = spu_sel(*vmax, *vmax1, vu);
            
          vh    = (vector signed short *) gh;
          vgh   = spu_add(*vh, 0);
          vgh   = spu_shuffle(vgh, vgh, *flip); /* flip elements of vgh */
          vu    = spu_cmpgt(vgh, *vmax);
/* g */   *vmax = spu_sel(*vmax, vgh, vu);
          
          vv_0    = (vector signed int *) gv;
          vv_1    = (vector signed int *)(gv+4);
          vgv_0   = spu_add(*vv_0, 0);
          vgv_1   = spu_add(*vv_1, 0);
          
          vmax_0 = spu_insert((int)max[7], vmax_0, 0);
          vmax_0 = spu_insert((int)max[6], vmax_0, 1);
          vmax_0 = spu_insert((int)max[5], vmax_0, 2);
          vmax_0 = spu_insert((int)max[4], vmax_0, 3);
          vmax_1 = spu_insert((int)max[3], vmax_1, 0);
          vmax_1 = spu_insert((int)max[2], vmax_1, 1);
          vmax_1 = spu_insert((int)max[1], vmax_1, 2);
          vmax_1 = spu_insert((int)max[0], vmax_1, 3);
          
/* h */   vu_0   = spu_cmpgt(vgv_0, vmax_0);
          vmax_0 = spu_sel(vmax_0, vgv_0, vu_0);
          vu_1   = spu_cmpgt(vgv_1, vmax_1);
          vmax_1 = spu_sel(vmax_1, vgv_1, vu_1);
          
/* h */   max[0] = (short) spu_extract(vmax_1, 3);
          max[1] = (short) spu_extract(vmax_1, 2);
          max[2] = (short) spu_extract(vmax_1, 1);
          max[3] = (short) spu_extract(vmax_1, 0);
          max[4] = (short) spu_extract(vmax_0, 3);
          max[5] = (short) spu_extract(vmax_0, 2);
          max[6] = (short) spu_extract(vmax_0, 1);
          max[7] = (short) spu_extract(vmax_0, 0);
          
          vu    = spu_cmpgt(*vmax, vsmax);
/* i */   vsmax = spu_sel(vsmax, *vmax, vu);
#endif

#ifdef DEBUG_INFO
          for(k=0; k<DIAG_LEN; k++)
            if(max[k] & 0x8000) printf("sw_db_spu2: cell score exceeds maximum...\n");
#endif
          
#ifndef VECTOR
          for(k=0; k<DIAG_LEN; k++)
            swTAB[k+1][DIAG_LEN-k+3] = max[k];
            
          for(k=0; k<DIAG_LEN; k++)
            max2[i+k-1] = max[k] - g;
#else
          swTAB[1][11] = max[0];
          swTAB[2][10] = max[1];
          swTAB[3][9]  = max[2];
          swTAB[4][8]  = max[3];
          swTAB[5][7]  = max[4];
          swTAB[6][6]  = max[5];
          swTAB[7][5]  = max[6];
          swTAB[8][4]  = max[7];

          *vmax2 = spu_add(*vmax, vg);

          vsw = spu_insert(swTAB[0][DIAG_LEN+3], vsw, 0);
          vsw = spu_insert(swTAB[1][DIAG_LEN+2], vsw, 1);
          vsw = spu_insert(swTAB[2][DIAG_LEN+1], vsw, 2);
          vsw = spu_insert(swTAB[3][DIAG_LEN  ], vsw, 3);
          vsw = spu_insert(swTAB[4][DIAG_LEN-1], vsw, 4);
          vsw = spu_insert(swTAB[5][DIAG_LEN-2], vsw, 5);
          vsw = spu_insert(swTAB[6][DIAG_LEN-3], vsw, 6);
          vsw = spu_insert(swTAB[7][DIAG_LEN-4], vsw, 7);
#endif
          /* main fill-the-middle loop */
          for(j=jj+DIAG_LEN+1, kk=DIAG_LEN+1; j<jj+chunk+1; j++, kk++)
          {
#ifndef VECTOR
            for(k=0; k<DIAG_LEN; k++) /* i & j offset */
/* a */       max[k] = 0x0;
              
            for(k=0; k<DIAG_LEN; k++)
/* b */       max1[k] = swTAB[k][kk-k-1+3] + sslookup[0x1F & db[n%2][kk-k-1]][q[k]];
            
            for(k=0; k<DIAG_LEN; k++)
/* c */       max3[k] = swTAB[k][kk-k+3] - g;
              
            for(k=0; k<DIAG_LEN; k++)
/* d */       gh[kk-k-1] = (max2[i+k-1] > gh[kk-k-2]-e)? max2[i+k-1]: gh[kk-k-2]-e;
            
            for(k=0; k<DIAG_LEN; k++)
/* e */       gv[kk-k-1] = (max3[k]     > gv[kk-k-1]-e)? max3[k]    : gv[kk-k-1]-e;
              
            for(k=0; k<DIAG_LEN; k++)
/* f */       if(max1[k]>max[k]) max[k] = max1[k];
              
            for(k=0; k<DIAG_LEN; k++)
/* g */       if(gh[kk-k-1]>max[k]) max[k] = gh[kk-k-1];
              
            for(k=0; k<DIAG_LEN; k++)
/* h */       if(gv[kk-k-1]>max[k]) max[k] = gv[kk-k-1];
              
            for(k=0; k<DIAG_LEN; k++)
/* i */       if(max[k]>smax) smax = max[k];

            for(k=0; k<DIAG_LEN; k++)
/* j */       swTAB[k+1][kk-k+3] = max[k];
              
            for(k=0; k<DIAG_LEN; k++)
/* k */       max2[i+k-1] = max[k] - g;

#else
/* a */     *vmax = spu_splats(zero);

            vss = spu_insert(sslookup[0x1F & db[n%2][kk-1]][q[0]], vss, 0);
            vss = spu_insert(sslookup[0x1F & db[n%2][kk-2]][q[1]], vss, 1);
            vss = spu_insert(sslookup[0x1F & db[n%2][kk-3]][q[2]], vss, 2);
            vss = spu_insert(sslookup[0x1F & db[n%2][kk-4]][q[3]], vss, 3);
            vss = spu_insert(sslookup[0x1F & db[n%2][kk-5]][q[4]], vss, 4);
            vss = spu_insert(sslookup[0x1F & db[n%2][kk-6]][q[5]], vss, 5);
            vss = spu_insert(sslookup[0x1F & db[n%2][kk-7]][q[6]], vss, 6);
            vss = spu_insert(sslookup[0x1F & db[n%2][kk-8]][q[7]], vss, 7);
           
            vdb = spu_rlmaskqwbyte(vdb, -2); /* shift right 2 bytes */
            vdb = spu_insert((signed int)db[n%2][kk-1], vdb, 0);
            vdb = spu_and(vdb, 0x001F);
 
/* b */     *vmax1 = spu_add(vss, vsw);
            
            /* vsw for max3, and for next loop */
            vsw = spu_insert(swTAB[0][kk+3], vsw, 0);
            vsw = spu_insert(swTAB[1][kk+2], vsw, 1);
            vsw = spu_insert(swTAB[2][kk+1], vsw, 2);
            vsw = spu_insert(swTAB[3][kk  ], vsw, 3);
            vsw = spu_insert(swTAB[4][kk-1], vsw, 4);
            vsw = spu_insert(swTAB[5][kk-2], vsw, 5);
            vsw = spu_insert(swTAB[6][kk-3], vsw, 6);
            vsw = spu_insert(swTAB[7][kk-4], vsw, 7);
              
/* c */     *vmax3 = spu_add(vsw, vg);
              
            vgh = spu_add(vgh, ve);
            vu  = spu_cmpgt(*vmax2, vgh);
/* d */     vgh = spu_sel(vgh, *vmax2, vu);

#ifdef TRACE
interval = trace_user_interval_entry();
#endif
            vgv_0 = spu_shuffle(vgv_0, vgv_1, *select4);
            vgv_0 = spu_add(vgv_0, vei);
            vgv_1 = spu_slqwbyte(vgv_1, 4); /* shift left 4 bytes */
            vgv_1 = spu_insert(gv[kk-DIAG_LEN+7], vgv_1, 3);
            vgv_1 = spu_add(vgv_1, vei);
            
            vmax3_0 = spu_insert((int)max3[7], vmax3_0, 0);
            vmax3_0 = spu_insert((int)max3[6], vmax3_0, 1);
            vmax3_0 = spu_insert((int)max3[5], vmax3_0, 2);
            vmax3_0 = spu_insert((int)max3[4], vmax3_0, 3);
            vmax3_1 = spu_insert((int)max3[3], vmax3_1, 0);
            vmax3_1 = spu_insert((int)max3[2], vmax3_1, 1);
            vmax3_1 = spu_insert((int)max3[1], vmax3_1, 2);
            vmax3_1 = spu_insert((int)max3[0], vmax3_1, 3);
            
/* e */     vu_0  = spu_cmpgt(vmax3_0, vgv_0);
            vgv_0 = spu_sel(vgv_0, vmax3_0, vu_0);
            vu_1  = spu_cmpgt(vmax3_1, vgv_1);
            vgv_1 = spu_sel(vgv_1, vmax3_1, vu_1);

            gv[kk-8] = spu_extract(vgv_0, 0);

#ifdef TRACE       
trace_user_interval_exit(interval, &pl);
#endif
            vu    = spu_cmpgt(*vmax1, *vmax);
/* f */     *vmax = spu_sel(*vmax, *vmax1, vu);
            
            vu    = spu_cmpgt(vgh, *vmax);
/* g */     *vmax = spu_sel(*vmax, vgh, vu);

            vmax_0 = spu_shuffle((vector signed int)(*vmax), vmax_0, *select);
            vmax_1 = spu_shuffle((vector signed int)(*vmax), vmax_0, *select2);
            
/* h */     vu_0   = spu_cmpgt(vgv_0, vmax_0);
            vmax_0 = spu_sel(vmax_0, vgv_0, vu_0);
            vu_1   = spu_cmpgt(vgv_1, vmax_1);
            vmax_1 = spu_sel(vmax_1, vgv_1, vu_1);

            *vmax = spu_shuffle((vector signed short)vmax_0,
                                  (vector signed short)vmax_1, *select3);
            vu    = spu_cmpgt(*vmax, vsmax);
/* i */     vsmax = spu_sel(vsmax, *vmax, vu);

/* j */     swTAB[1][kk+3] = max[0];
            swTAB[2][kk+2] = max[1];
            swTAB[3][kk+1] = max[2];
            swTAB[4][kk]   = max[3];
            swTAB[5][kk-1] = max[4];
            swTAB[6][kk-2] = max[5];
            swTAB[7][kk-3] = max[6];
            swTAB[8][kk-4] = max[7];
            
/* k */     *vmax2 = spu_add(*vmax, vg);
#endif

#ifdef DEBUG_INFO
            for(k=0; k<DIAG_LEN; k++)
              if(max[k] & 0x8000) printf("sw_db_spu3: cell score exceeds maximum...\n");
#endif
          }
#ifdef VECTOR
          gh[kk-2] = spu_extract(vgh, 0);
          gh[kk-3] = spu_extract(vgh, 1);
          gh[kk-4] = spu_extract(vgh, 2);
          gh[kk-5] = spu_extract(vgh, 3);
          gh[kk-6] = spu_extract(vgh, 4);
          gh[kk-7] = spu_extract(vgh, 5);
          gh[kk-8] = spu_extract(vgh, 6);
          gh[kk-9] = spu_extract(vgh, 7);

          gv[kk-2] = spu_extract(vgv_1, 3);
          gv[kk-3] = spu_extract(vgv_1, 2);
          gv[kk-4] = spu_extract(vgv_1, 1);
          gv[kk-5] = spu_extract(vgv_1, 0);
          gv[kk-6] = spu_extract(vgv_0, 3);
          gv[kk-7] = spu_extract(vgv_0, 2);
          gv[kk-8] = spu_extract(vgv_0, 1);
          //gv[kk-9] = spu_extract(vgv_0, 0); /*this one saved in loop above */
#endif
          gapH[i-1] = gh[kk-2];
										
          /* fill lower right corner */
          j = jj + chunk + 1;
          for(k=1; k<DIAG_LEN; k++) /* i offset */
          {
            for(kk=j-k, m=0; kk<j; kk++, m++) /* j & offset */
            {
              max[k]  = 0x0;
              max1[k] = swTAB[k][chunk-k+m+3] + sslookup[0x1F & db[n%2][chunk-k+m]][q[k]];
              max3[k] = swTAB[k][chunk-k+m+1+3] - g;
      
              gh[chunk-k+m] = (max2[i+k-1] > gh[chunk-k+m-1]-e)? max2[i+k-1]: gh[chunk-k+m-1]-e;
              gv[chunk-k+m] = (max3[k]     > gv[chunk-k+m]-e)  ? max3[k]    : gv[chunk-k+m]-e;

              if(max1[k]>max[k])       max[k] = max1[k];
              if(gh[chunk-k+m]>max[k]) max[k] = gh[chunk-k+m];
              if(gv[chunk-k+m]>max[k]) max[k] = gv[chunk-k+m];
              if(max[k]>smax)            smax = max[k];
#ifdef DEBUG_INFO
              if(max[k] & 0x8000) printf("sw_db_spu4: cell score exceeds maximum...\n");
#endif
              swTAB[k+1][chunk-k+m+1+3] = max[k];
              max2[i+k-1] = max[k] - g;
            }
            gapH[i+k-1] = gh[chunk-k+m-1];
          }
          lastI = i + DIAG_LEN - 1;
        
          for(k=0; k<DIAG_LEN; k++)
            temp[(n+1)%2][i+k] = swTAB[k+1][chunk+3];
        }
      
        /* do bottom rows that parallel section doesn't get */
        for(k=0; k<chunk+1; k++)
          swTAB[1][k+3] = swTAB[DIAG_LEN][k+3]; /* copied to row 0 below */
        
        for(; i<qLen+1; i++)
        {
          for(k=0; k<chunk+1; k++) /* copy bottom row to top */
            swTAB[0][k+3] = swTAB[1][k+3];
          
          q[0] = 0x1F & ((char *)query)[i-1];
          swTAB[0][0+3] = temp[n%2][i-1];
      
          for(j=jj+1, kk=0; j<jj+chunk+1; j++, kk++)
          {
            max[0]  = 0x0;
            max1[0] = swTAB[0][kk+3] + sslookup[0x1F & db[n%2][kk]][q[0]];
            max3[0] = swTAB[0][kk+4] - g;
      
            if(j == jj+1)
              gh[kk] = (max2[i-1] > gapH[i-1]-e)? max2[i-1]: gapH[i-1]-e;
            else
              gh[kk] = (max2[i-1] > gh[kk-1]-e) ? max2[i-1]: gh[kk-1]-e;
          
            gv[kk] = (max3[0] > gv[kk]-e)? max3[0]: gv[kk]-e;
          
            if(max1[0]>max[0]) max[0] = max1[0];
            if(gh[kk]>max[0])  max[0] = gh[kk];
            if(gv[kk]>max[0])  max[0] = gv[kk];
            if(max[0]>smax)      smax = max[0];
#ifdef DEBUG_INFO
            if(max[0] & 0x8000) printf("sw_db_spu5: cell score exceeds maximum...\n");
#endif
            swTAB[1][kk+4] = max[0];
            max2[i-1] = max[0] - g;
          }
          gapH[i-1] = gh[kk-1];
          lastI = i;
          temp[(n+1)%2][i] = swTAB[1][chunk+3];
        }

        /* make sure the double buffer finished */
        mfc_write_tag_mask(0x100);
        mfc_read_tag_status_all();
      }
    
      /* last strip */
      if(lastI)
      {
        for(i=0; i<dbLen+2 - savj; i++)
          last[i] = 0;
 
        if(!gotLast)
        {
          mfc_get((void *)db[n%2], dbAdd+savj-1, ((dbLen+17-savj)/16)*16, 3, 0, 0);
          mfc_write_tag_mask(0x1000);
          mfc_read_tag_status_all();
        }

        for(k=0; k<dbLen+2-savj; k++)
          gv[k] = -g - (savj+k-1)*e;

        for(i=1; i<lastI+1; i++)
        {
          q[0] = 0x1F & ((char *)query)[i-1];
          last[0] = temp[n%2][i-1];
        
          for(j=savj, k=0; j<dbLen+1; j++, k++)
          {
            max[0]  = 0x0;
            max1[0] = last[k] + sslookup[0x1F & db[n%2][k]][q[0]];
            last[k] = max2[i-1] + g;
        
            max3[0] = last[k+1] - g;

            if(k)
              gh[k] = (max2[i-1] > gh[k-1]-e)  ? max2[i-1]: gh[k-1]-e;
            else
              gh[0] = (max2[i-1] > gapH[i-1]-e)? max2[i-1]: gapH[i-1]-e;
        
            gv[k] = (max3[0] > gv[k]-e)? max3[0]: gv[k]-e;

            if(max1[0]>max[0]) max[0] = max1[0];
            if(gh[k]>max[0])   max[0] = gh[k];
            if(gv[k]>max[0])   max[0] = gv[k];
            if(max[0]>smax)      smax = max[0];
#ifdef DEBUG_INFO
            if(max[0] & 0x8000) printf("sw_db_spu6: cell score exceeds maximum...\n");
#endif
            max2[i-1]   = max[0] - g;
          }
          last[dbLen - savj + 1] = max[0];
        }
      }
#ifdef VECTOR
      for(k=0; k<VECT_LEN; k++)
      {
        i = spu_extract(vsmax, k);
        if(i > smax) smax = i;
      }
#endif
    }
    else /* problem was too small for strips */
    {
      /* init */
      mfc_get((void *)db[0], dbAdd, ((dbLen+15)/16)*16, 1, 0, 0);
      mfc_write_tag_mask(0x10);
      mfc_read_tag_status_all();
      
      for(i=0; i<dbLen+1; i++)
      {
        gv[i]   = -g - (i-1)*e;
        last[i] = 0;
      }
      gh[0] = gv[0];
    
      for(i=1;i<qLen+1; i++)
      {
        q[0] = 0x1F & ((char *)query)[i-1];
      
        for(j=1; j<dbLen+1; j++)
        {
          max[0]  = 0x0;
          max1[0] = last[j-1] + sslookup[0x1F & db[0][j-1]][q[0]];
          last[j-1] = max2[i-1] + g;

          max3[0] = last[j] - g;

          gh[j] = (max2[i-1] > gh[j-1]-e)? max2[i-1]: gh[j-1]-e;
          gv[j] = (max3[0]   > gv[j]-e)  ? max3[0]  : gv[j]-e;

          if(max1[0]>max[0]) max[0] = max1[0];
          if(gh[j]>max[0])   max[0] = gh[j];
          if(gv[j]>max[0])   max[0] = gv[j];
          if(max[0]>smax)      smax = max[0];
#ifdef DEBUG_INFO
          if(max[0] & 0x8000) printf("sw_db_spu7: cell score exceeds maximum...\n");
#endif
          max2[i-1] = max[0] - g;
        }
        last[dbLen] = max[0];
      }
    }
    
    spu_write_out_mbox(smax);
    
#ifdef TIMING
    clock = spu_clock_read();
    if(clock < start)
      totTime += UINT_MAX - start + clock;
    else
      totTime += clock - start;
#endif
  }
  
finish:

#ifdef TIMING
  printf("spe %d: total sec = %llu\n", cb.speId, totTime/FREQ);
  printf("spe %d: total sec = %llu (%lf)\n", cb.speId, totTime/FREQ, (double)totTime/(double)FREQ);
  printf("spe %d: clock exit code = %d\n", cb.speId, spu_clock_stop());
#endif

  return 0;
}
