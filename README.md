# propulsion LDVM


This repo is made to draw the pareto front between trust and efficiency of a flapping wing.
The pitch-plunge kinematics are optimized using the evolutionary algorithm NSGA2. The aerodynamic model is the LDVM model publeshed by Ramesh et al. in 2013. Original code is findable in the folder LDVM_v2_original.5. 
The orginal code is reimplemented in python with no performance loss compared to the original flutter implementation.
The NSGA2 implementation is issued from the pymoo library (Python library.)


The repo is structured as follows. 
- `LDVM_v2_original.5/`: Contains the original LDVM model implementation by Ramesh et al. (2013).
- `paper_result/`: gather data obtained by Anderson et al. 1998 for several flapping kinematics. LDVM validation can be made against those data.
- `LDVM`: Contains the code developed for the project.
    - `launsher.sh`: launsh script on slurm supercomputer.
    - `ldvm_back_up.py`: python implementaion of the LDVM method.
    - `validation_ldvm.py`: validation against Theodorsen model and Anderson data 
    - `naca0015_airfoil.dat/sd7003.dat/sd7012.dat`: data points for several profile shapes.
- `LDVM/NSGA2`: Contains the code for kinematics optimization.
    - `pareto_results`: Folder Gathering optim results.
    - `problem_multi_obj.py`: definition of thge optimization problem.
    - `pareto_optim.py` : NSGA2 run. Calculation of the pareto front.
    - `callbacks.py`: defines the post processing to make along with NSGA2 running.
    - `postprocessing.py` : Includes pareto visualization/ flow animation and display optimal designs.
    - `config.yaml`: parameterize the simulation and optimization problem.
    - `launsher.sh`: launsh job on slurm systems.


- `requirements.txt`: all libraries used for computations.