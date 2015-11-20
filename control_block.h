/*
 control block contains the address of things in main memory needed by the spes
*/
typedef struct _control_block
{
  unsigned int addr_matrix;
  unsigned int addr_qSeq;
  unsigned int qLen;
  unsigned int eg;
  unsigned int eog;
  unsigned int speId;
  unsigned char pad[104]; /* pad to a 128-byte cache line */
  
} control_block __attribute__ ((aligned(16)));
