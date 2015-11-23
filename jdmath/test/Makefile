TEST_PGMS = test_vector test_crout test_bisect
JDMATH_LIB_DIR = ../$(ARCH)objs

INCS = ..
LIBS = -L$(JDMATH_LIB_DIR) -ljdmath -lm

JDMATH_LIB = $(JDMATH_LIB_DIR)/libjdmath.a

runtests: $(TEST_PGMS)
	@echo ""
	@echo "Running tests:"
	@echo ""
	@for X in $(TEST_PGMS); \
	do \
	  ./$$X ; \
	done

test_crout: $(JDMATH_LIB) test_crout.c matrix.c
	$(CC) $(CFLAGS) -o test_crout test_crout.c -I$(INCS) $(LIBS)

test_vector: $(JDMATH_LIB) test_vector.c
	$(CC) $(CFLAGS) -o test_vector test_vector.c -I$(INCS) $(LIBS)

test_bisect: $(JDMATH_LIB) test_bisect.c
	$(CC) $(CFLAGS) -o test_bisect test_bisect.c -I$(INCS) $(LIBS)

$(JDMATH_LIB):
	cd ..; $(MAKE)


clean: 
	-/bin/rm -f *~ *.o $(TEST_PGMS)
