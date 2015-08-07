HTSDIR=../../htslib
sqlite: sqlite.c
	gcc -o sqlite -I../ sqlite.c -I$(HTSDIR) -Wl,-Bstatic -L../ -lbcftools -L$(HTSDIR) -lhts -Wl,-Bdynamic -lpthread -lz -lm -ldl -lsqlite3

