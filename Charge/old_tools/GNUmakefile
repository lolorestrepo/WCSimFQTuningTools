NOSKLIBRARIES = 1
HEMI_CUDA_DISABLE=1

ifndef NOSKLIBRARIES
ifndef ATMPD_ROOT
ATMPD_ROOT=$(FITQUN_ROOT)
endif
include $(ATMPD_ROOT)/config.gmk
else
include $(FITQUN_ROOT)/config.gmk
endif	
#
#  Objects
#
LOCAL_LIBS = -lfiTQun -lWCSim
ROOTLIBS  := $(shell root-config --glibs) -lTreePlayer -lMinuit

OBJS =

#LOCAL_INC = -I$(ATMPD_ROOT)/src/recon/fitqun
LOCAL_INC = -I$(FITQUN_ROOT)
LDLIBS += -L$(FITQUN_ROOT) -L$(WCSIMDIR)
#
#  Rules for building library 
#

CXXFLAGS += -DHEMI_CUDA_DISABLE

.PHONY:  clean exec

exec: chargePDFdiagnostics

chargePDFdiagnostics: chargePDFdiagnostics.o ${OBJS}
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJS) $(LDLIBS)

clean:
	$(RM) -f *.o *.a *~ chargePDFdiagnostics

