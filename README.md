# Propulsion LDVM

This repo is made to draw the Pareto front between thrust and efficiency of a flapping wing.  
The pitch–plunge kinematics are optimized using the NSGA-II evolutionary algorithm.  
The aerodynamic model is the LDVM model published by Ramesh et al. in 2013. The original code can be found in the folder `LDVM_v2_original.5`.

The original code has been reimplemented in Python with no performance loss compared to the original Flutter implementation.  
The NSGA-II implementation is provided by the *pymoo* Python library.

## Repository Structure

- `LDVM_v2_original.5/`: Contains the original LDVM model implementation by Ramesh et al. (2013).
- `paper_result/`: Gathers data obtained by Anderson et al. (1998) for several flapping kinematics. LDVM validation can be performed against these data.
- `LDVM/`: Contains the code developed for the project.
  - `launcher.sh`: Launch script for SLURM systems.
  - `ldvm_back_up.py`: Python implementation of the LDVM method.
  - `validation_ldvm.py`: Validation against the Theodorsen model and Anderson's data.
  - `naca0015_airfoil.dat / sd7003.dat / sd7012.dat`: Airfoil coordinate files for several profile shapes.
- `LDVM/NSGA2`: Contains the code for kinematics optimization.
  - `pareto_results/`: Folder gathering optimization results.
  - `problem_multi_obj.py`: Definition of the optimization problem.
  - `pareto_optim.py`: NSGA-II execution and Pareto front computation.
  - `callbacks.py`: Defines the post-processing performed during NSGA-II runs.
  - `postprocessing.py`: Contains Pareto visualization, flow animation, and optimal design display.
  - `config.yaml`: Parameters for the simulation and optimization problem.
  - `launcher.sh`: Launch script for SLURM systems.

- `requirements.txt`: Lists all libraries used for the computations.
