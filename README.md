# Fault Level in Inverter-Based Power Grids

Souce code for the IBR fault analysis algorithm developed as part of my MEng FYP: "Fault Level in Inverter-Based Power
Grids".

Written and tested using MATLAB R2024b. Must install MATPOWER before use. Make sure both the "Algorithms" and "Case Files" folders are in the path. Full user guide provided in the FYP report.

Algorithms:
* algorithm_v5 - Original iterative algorithm
* algorithm_v6 - Newer version with G74 recommendations implemented
* run_v7 - Final verson in a modular form that can be repeatedly called by another script with different inputs
* G74_algorithm - An implementation of the non-iterative G74 method for comparison

Case Files:
* case39bus - Original IEEE 39-bus system included in MATPOWER library. Modified to also assemble admittance matrix and run powerflow.
* case39bus_2 - Version with some SMs replaced with IBRs
* case39bus_3 - Same as case39bus_2 but can now pass in parameters for "wind" and "solar" to change loading automatically
* case39bus_4 - Copy of case39bus used for temporary testing
* case39bus_5 - Modified case39bus with low IBR penetration
* case39bus_6 - Modified case39bus with high IBR penetration
* case39bus_7 - Modified case39bus with medium IBR penetration
* case57bus - IEEE 57-bus system
* case57bus_2 - Modified case57bus with some SMs removed and some IBRs added