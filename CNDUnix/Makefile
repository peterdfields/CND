# makefile for CND CytoNuclearDisequilibria
#
SOURCES= Dmain.c Ddatain.c Dstats.c  \
		 Mmain.c Mdatain.c  Mstats.c  \
		 chisqr.c utils.c exact.c
HEADERS= Dmain.h Mmain.h utils.h chisqr.h
LINTFILES= $(SOURCES:.c=.ln)
OBJECTS= $(SOURCES:.c=.o)
LINTFLAGS= -chu 
LDFLAGS= -lm
#
#  if your compiler doesn't have div and ldiv, then put
#  CPPFLAGS= -DDIVT
# 
CPPFLAGS= 
CFLAGS=
CC=cc
# You can change the BINDIR to whatever is appropriate for your system.
BINDIR= ../bin
#
# You may need to delete these lines if you are not on a Sun
#
.KEEP_STATE:
all debug optim : CNDm CNDd
debug := CFLAGS = -g
optim := CFLAGS = -O
#
# Dependancies
#
utils.o Mmain.o Mstats.o Mdatain.o Dmain.o Dstats.o Ddatain.o : utils.h
Mmain.o Mdatain.o Mstats.o : Mmain.h
Dmain.o Ddatain.o Dstats.o : Dmain.h
#
# Targets
#
CNDm : Mmain.o Mdatain.o Mstats.o utils.o exact.o chisqr.o
	$(CC) -o CNDm Mmain.o Mdatain.o Mstats.o utils.o exact.o chisqr.o $(LDFLAGS)

CNDd : Dmain.o Ddatain.o Dstats.o utils.o exact.o chisqr.o
	$(CC) -o CNDd Dmain.o Ddatain.o Dstats.o utils.o exact.o chisqr.o $(LDFLAGS)
#
install: CNDm CNDd
	cp CNDm CNDd $(BINDIR)
#
lint: $(LINTFILES)
	$(LINT.c) $(LINTFILES)

#
veryclean: clean
	rm -f  CNDm.out CNDm.err CNDd.out CNDd.err CNDm CNDd 
#
clean:
	rm -f $(LINTFILES) $(OBJECTS) *.BAK core
