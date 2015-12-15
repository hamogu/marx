#--*-sh-*-------------------------------------------------------------------
#  This Makefile is used to create the programs bdat2tbl and tbl2bdat.  They
#  are used to dump the binary .dat files to and from ASCII tables.
#  This Makefile is written with the directory structure of the source tree,
#  not the installed tree in mind, because these tools are only
#  relevant for MARX developers.
#---------------------------------------------------------------------------
#CC = gcc
#CFLAGS = -g

JDMATH_INC = ../../jdmath/src
JDMATH_LIB = $(JDMATH_INC)/..

LIBS = $(LDFLAGS) -L$(JDMATH_LIB) -ljdmath -lm

all: bdat2tbl tbl2bdat

bdat2tbl: bdat2tbl.c
	$(CC) $(CFLAGS) -I$(JDMATH_INC) bdat2tbl.c -o bdat2tbl $(LIBS)
tbl2bdat: tbl2bdat.c
	$(CC) $(CFLAGS) -I$(JDMATH_INC) tbl2bdat.c -o tbl2bdat $(LIBS)

clean:
	-/bin/rm -f *~ bdat2tbl tbl2bdat
