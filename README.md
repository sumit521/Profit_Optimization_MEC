# Profit_Optimization_MEC

In this project we try to use genetic algorithm to solve decoupled optimization problem of deciding offloading strategy and allocation of radio and computational resources. 

[![View Profit_Optimization_MEC on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/95648-profit_optimization_mec)

**New updates:**
Updated the code files for simulating the paper https://ieeexplore.ieee.org/document/10273406. 
Code for "Profit Optimization in MEC using genetic algorithm" can be accessed from the last release of the repository.

Use Main_Script_final.m file  to get the figure for the comparison of algorithms which is the key result of the paper.
Use Main_script_experiment_parameter.m file to generate the figures corresponding to the impact of changes in resource availability on profitability. See the code around "%%Calculating the profits" section and "%% for resource dependency trend" section. Assign the values from vectors ln_vector, B_vector, F_vector to respective variables in different rounds of simulation and also make necessary changes in form of uncommenting/commenting and changing the legend in "%% for resource dependency trend" section of the code.  

Similarly, use Main_script_experiment_parameter.m file  to get the impact of parameters of evolutionary algorithms on the profitability. See the code around "%%Calculating the profits" section and "%% plotting the results" section.

