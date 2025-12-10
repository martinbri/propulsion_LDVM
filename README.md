# propulsion LDVM


This repo is made to draw the pareto front between trust and efficiency of a flapping wing.
The pitch-plunge kinematics are optimized using the evolutionary algorithm NSGA2. The aerodynamic model is the LDVM model publeshed by Ramesh et al. in 2013. Original code is findable in the folder LDVM_v2_original.5. 
The orginal code is reimplemented in python with no performance loss compared to the original flutter implementation.
The NSGA2 implementation is issued from the pymoo library (Python library.)


The repo is structured as follows. 
- `LDVM_v2_original.5/`: Contains the original LDVM model implementation by Ramesh et al. (2013).
- `python_reimplementation/`: Python reimplementation of the LDVM model.
- `optimization/`: Scripts and configurations for NSGA2-based optimization.
- `results/`: Stores generated Pareto fronts and analysis outputs.
- `docs/`: Documentation and additional resources.