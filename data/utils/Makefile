INCS = -I../../libsrc -I../../../pfile/src -I../../../jdmath/src
marxpipe: marxpipe.c
	$(CC) $(CFLAGS) $(INCS) -o marxpipe marxpipe.c 
clean:
	/bin/rm -f *~ marxpipe
ea_corr:
	@echo "To create the effective area correction files, run:"
	@echo slsh ./mk_hrma_ea.sl
	@echo sh ./run_mkarf.sh
	@echo slsh ./mkcorr.sl
	@echo tbl2bdat corr_1.tbl corr_1.dat
	@echo tbl2bdat corr_3.tbl corr_3.dat
	@echo tbl2bdat corr_4.tbl corr_4.dat
	@echo tbl2bdat corr_6.tbl corr_6.dat
