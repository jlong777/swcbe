CC=ppu-gcc -m32
SPECC=spu-gcc
CFLAGS=-O2 -Wall
DEBUG=#-g
INSTALL_DIR=/usr/local/bin

all: sw_db sw_db_spu
sw_db: sw_db.c control_block.h
	$(CC) $(CFLAGS) $(DEBUG) -o sw_db sw_db.c -lspe2 -lpthread -lcbl
sw_db_spu: sw_db_spu.c control_block.h
	$(SPECC) $(CFLAGS) $(DEBUG) -o sw_db_spu sw_db_spu.c
install:
	if [ ! -e ${INSTALL_DIR} ]; then mkdir -p ${INSTALL_DIR}; fi
	mv sw_db     ${INSTALL_DIR}
	mv sw_db_spu ${INSTALL_DIR}
sw_db_spu_timing: sw_db_spu.c control_block.h
	$(SPECC) $(CFLAGS) $(DEBUG) -o sw_db_spu sw_db_spu.c \
        -L/opt/cell/sdk/prototype/usr/spu/lib -lsputimer
timing: sw_db sw_db_spu_timing
clean:
	rm sw_db sw_db_spu
        
