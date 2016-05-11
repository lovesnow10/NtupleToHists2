PROGNAME      = NtupleToHists2
SRCDIR = src
INCDIR = include
SOURCEFILES   = $(wildcard $(SRCDIR)/*.cxx)
#OBJS          = $(patsubst %.cxx, %.o, $(SOURCEFILES))
OBJS          = $(patsubst $(SRCDIR)/%.cxx, %.o, $(SOURCEFILES))
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
LDFLAGS       = -O
LIBS         += $(ROOTLIBS)
CFLAGS       += $(ROOTCFLAGS)
CFLAGS       += -I$(INCDIR)

#  Not sure why Minuit, TMV aren't being included -- put in by hand
#
LIBS         += -lMinuit
LIBS         += -lTMVA -lTreePlayer
#LIBS          += -L/nfs/scratch0/cowan/TMVA/lib -lTMVA

%.o: $(SRCDIR)/%.cxx 
	clang++ ${CFLAGS} -c  -g -o $@ $<

$(PROGNAME):    $(OBJS)
	clang++ -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

test:
	@echo $(CFLAGS)

clean:	
	-rm -f ${PROGNAME} ${OBJS}
