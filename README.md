# Group 4 Power Sim
The following Power Simulator is a simulation tool for modeling and analyzing power systems. The following tool only works for a 7 bus system. A generator will be placed on bus's 1 and 7, where the slack bus can be set as either of those buses. The simulator is written in Python scripts using Python version 3.7. Each script implements various components of the power system. There are 15 important files in this directory, however the user will not interact with every one. The following lists an overview of the files in the repository.

Files:
* Bus.py: This file implements the Bus class, which represents a bus in the power system.
* DC_Power_Flow_Solver.py: This file implements the DC Power FLow Solver, which solves the power flow equations for a DC power system.
* Group 4 Simulator Documentation.docx: This file is simply the documenation that describes the files in more detail for the user to better understand the program in a Word document format.
* Group 4 Simulator Documentation.pdf: This is the same information that is in the Group 4 Simulator Documentation.docx file, just in pdf format.
* Fast_Decoupled_Solver.py: This file implements the Fast Decoupled Power Flow Solver, which solves the power flow equations for an AC power system using the Fast Decoupled method.
* Fault_Calculation.py: This file implements the Fault Calculation class, which calculates fault currents and voltages in the power system.
* Generator.py: This file implements the Generator class, which represents a generator in the power sytem.
* Grid.py: This file implements the Grid class, which represents the entire power system and is where the componenets of the system will be stored.
* Main.py: This file is the main script that runs the simulation. The user will interact with this file to establish what data for the simulation they want to run.
* Newton_Rhapson_Power_Flow_Solver: This file implements the Newton-Rhapson Power Flow Solver, which solves the power flow equations for an AC power system using the Newton-Rhapson method.
* Sequence_Networks.py: This file implements the Sequence Networks class, which calculates the positive, negative, and zero sequence impedances in the power system.
* Transformer.py: This file implements the Transformer class, which represents a transformer in the power system.
* TransmissionLine.py: This file implements the TransmissionLine class, which represents a transmission line in the power system.
* TransmissionLineBundles.py: This file implements the TransmissionLineBundles class, which holds data on the type of transmission line is being used in the system.
* Zero_Impedance_Network.pdf: This file is a pdf drawing of the zero-impedance sequence network as described in Milestone 5 for project submission.

Usage:
To use the Power Simulator, you will need to install the following Python libraries: numpy, pandas, typing, sys, and math.
Once installed, you can run the Main.py file script. There are preset simulator parameters that can be altered to the user's case, but they are there as an example of how the file can be used as a base case. A user will first create a Grid where everything will be stored. Generators and transformers will be added along with their bus location in the system. Transmission lines will then be added to connect the different buses together as desired. After the user is done adding componenets, the Ybus will be calculated. Bus data will then be set according to their types and desired parameters. Next, the user can select what type of Power Flow solver they would like to implement based by choosing either the Newton Rhapson, Fast Decoupled, or DC Power flow solver. Finally, the user can choose to solve the sequence networks and perform fault calculations if desired. Parameters of interest will be printed to the user's terminal.

Additional Usage Note: The user will first establish all desired buses, transmission lines, generators, and transformers along with their groundings and other important information. Then, the user can run the program where they can select the type of analysis they would like to conduct through input of the terminal, as well as enter necessary inputs regarding what type of fault analysis they would like to conduct.

Conclusion:
The Group 4 Power Simulator is a tool used for modeling and analyzing 7 bus power systems. Install and use the simulator to simulate a power system with your own desired parameters to fit your specific needs.
