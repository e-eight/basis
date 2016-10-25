$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := indexing
module_units_cpp-h := operator many_body
module_units_cpp-h += lsjt_scheme lsjt_operator jjjt_scheme jjjt_operator
module_units_cpp-h += jjjpn_scheme jjjpn_operator nlj_orbital
# module_units_f := 
module_programs_cpp := lsjt_scheme_test lsjt_operator_test jjjt_scheme_test
module_programs_cpp += jjjpn_scheme_test jjjpn_operator_test
module_programs_cpp += nlj_orbital_test
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
