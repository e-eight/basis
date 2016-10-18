################################
# project configuration
################################

libname = basis

# modules -- header-only
modules_h = indexing

# modules -- header-plus-object 
##modules_ho = jt_scheme
modules_ho = operator many_body
modules_ho += lsjt_scheme lsjt_operator jjjt_scheme jjjt_operator
modules_ho += jjjpnorb_scheme jjjpnorb_operator nlj_orbital

# programs
##programs = jt_scheme_test
programs = lsjt_scheme_test lsjt_operator_test jjjt_scheme_test
programs += jjjpnorb_scheme_test jjjpnorb_operator_test
programs += nlj_orbital_test
##programs += write_lsjt_relative
CC := $(CXX)

CXXFLAGS += -std=c++11

# set flag for linking to FORTRAN
# need_fortran = 

################################
# common definitions
################################

COMMON_MAKE_DIR ?= .
include $(COMMON_MAKE_DIR)/common.mk

################################
# options and dependencies
################################

# program linking
CC := $(CXX)

# external libraries
LDFLAGS += -L../am -L../mcutils
LDLIBS +=  -lam -lmcutils -lmcutils -lgsl

