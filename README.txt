This is a thermal-aware floorplanner for FPGA. Before starting, you need a resource file (.res) to specify FPGA resources, and a module file (.module) to specify modules to be mapped on the FPGA. Both files should be written according to Python3 syntax. 'test.res' and 'test.module' are given as examples.

The floorplanner can be executed with Python >= 3.4

Two floorplanning algorithms are implemented: simulated annealing (SA) and modified cluster growth (MCG). The type of floorplanning algorithm need to be specified in the option. To run the floorplanner:
    ./floorplanner.py -t [floorplanner_type] [design_name]
    
design_name is the name of .res and .module file, which is 'test' by default. floorplanner_type can be 'sa' or 'mcg'. There are other options to adjust parameter values in heuristic algorithms, also max thread limit in parallel, please read the code for details.

These file will be generated after floorplanning finished:
    output_files/empty_floorplan.svg    -   empty foorplan according to .res file
    output_files/[design_name]_floorplan_[sa/mcg].svg    -   final floorplan
    output_files/[design_name]_floorplan_[sa/mcg].svg    -   thermal map of the final floorplan
    output_files/json/finall_result   -   module position info of final floorplan in JSON format
    hotspot/[design_name].grid.steady   -   thermal map of the final floorplan raw data
    
The thermal simulator used in this floorplanner if HotSpot 6.0. The './hotspot/hospot' is the binary file of the simulator compiled under x86_64 CentOS 7. The license of HotSpot 6.0 is also included at './hotspot/LICENSE'.
