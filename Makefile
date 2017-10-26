CFLAGS = -g -Wall -Wno-parentheses -Wno-deprecated -O3 -DWORDINDEX_WITH_4_BYTE 

OBJS = Script/main.o Script/maxent.o Script/blmvm.o Script/ibm1.o

all: $(OBJS) extract
extract: $(OBJS)
	g++ $(CFLAGS) $(OBJS) -o $@
clean:
	/bin/rm -r -f $(OBJS) a.out *.o
.cpp.o:
	g++ -c $(CFLAGS) $< -o $@
