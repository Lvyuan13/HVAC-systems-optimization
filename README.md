# HVAC-systems-optimization
1.excel files includes the source data. 

coldload_huge.xlsx contains the demand cold load of buildings.chillers.xlsx contains the chillers working data,pump_constant.xlsx contains the constant frequency pump's working data. pump_var.xlsx contains the varible frequency pump's working data.

Don't delete or rename them.

2.HVAC_OPT_chiller.m is the optimization program to optimized the chillers' stop or start and flow disturbtion for chillers.

HVAC_OPT_pump.m is the optimization program to optimize the pumps' stop or start and flow disturbtion of pumps.

3.the two program(HVAC_OPT_chiller.m and HVAC_OPT_pump.m cannot run without OPTI toolbox, which is a open source library for optimization. see https://github.com/jonathancurrie/OPTI  to install it.
