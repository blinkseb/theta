include ../Makefile.options

sources := $(wildcard *.cpp)
objects := $(patsubst %.cpp,.bin/%.o,$(sources))
deps := $(patsubst %.cpp,.bin/%.d,$(sources))

include ../Makefile.rules

all: ../lib/simplex.so

../lib/simplex.so: $(objects)
	@$(LINK_SO) $(EXTRA_LINK)

#define clean target in any case:
clean-subdir:
	@rm -f ../lib/simplex.so
	
test: .bin/test.o .bin/test_main.o ../lib/simplex.so
	@$(LINK_EXE)
