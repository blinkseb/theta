objects := .bin/lbfgs.o
deps := .bin/lbfgs.d

include ../Makefile.rules

all: ../lib/liblbfgs.so

../lib/liblbfgs.so: $(objects)
	@$(CC) -o $@ -shared $(objects)

clean-subdir:
	@rm -f ../lib/liblbfgs.so
