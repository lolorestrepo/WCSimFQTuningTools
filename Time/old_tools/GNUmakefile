export FITQUN_ROOT=/pbs/home/g/gdiazlop/Software/fiTQun/fiTQun
export WCSIMDIR=/pbs/home/g/gdiazlop/Software/WCSim/install/


#ifndef ATMPD_ROOT
#	ATMPD_ROOT=$(FITQUN_ROOT)
#endif

#include $(ATMPD_ROOT)/config.gmk
include $(FITQUN_ROOT)/config.gmk
#
#  Objects
#

LOCAL_LIBS = -lfiTQun -lWCSim
ROOTLIBS  := $(shell root-config --glibs) -lTreePlayer -lMinuit

OBJS =
LOCAL_INC = 
#LOCAL_INC += -I$(ATMPD_ROOT)/src/recon/fitqun
LOCAL_INC += -I$(FITQUN_ROOT) -I$(WCSIMDIR)/include
#LOCAL_INC += -I$(SKOFL_ROOT)/include
LDLIBS += -L$(FITQUN_ROOT) -L$(WCSIMDIR)/lib
SITE_DEFINES += -DHEMI_CUDA_DISABLE


#
#  Rules for building library 
#

.PHONY:  clean exec

exec: makehist #scattpdf

makehist: makehist.o ${OBJS}
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJS) $(LDLIBS)

makehistWCSim: makehistWCSim.o ${OBJS}
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJS) $(LDLIBS)


makehistWCSimIndirect: makehistWCSimIndirect.o ${OBJS}
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJS) $(LDLIBS)

checkTResVsL: checkTResVsL.o ${OBJS}
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJS) $(LDLIBS)

scattpdf: scattpdf.o ${OBJS}
	$(FC) -o $@ scattpdf.o ${OBJS} $(LDLIBS)


clean:

	$(RM) -f *.o *.a *~ makehist scattpdf makehistWCSim makehistWCSimIndirect


