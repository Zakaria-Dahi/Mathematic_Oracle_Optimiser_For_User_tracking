# Takeover Time-Driven Adaptive Bi-Phased Evolutionary Algorithm 
**About:** This repository contains the source code of the Takeover Time-driven Adaptive bi-phased Evolutionary Algorithm (TD-EA) for solving the Mobility Management Problem (MMP) in cellular networks [1]. 

**Contact:** zakaria.dahi@uma.es, zakaria.dahi@univ-constantine2.dz

**How To:**
- To run the TD-EA, execute the file **```main.m```**.

**Notes:**
- It is prefered that you use Matlab R2014a to run the code as it is. 
- If you use newer versions of Matlab, you will probably need to update some obselete functions such as ```Mat2Cell```.

**File organisation:**
- **```main.m```:** Allows you to set the algorithm's/experimentation's parameters and launch an execution.
- **```Crossover.m```:** The crossover operator.
- **```Elitist_selection.m```:** The selection operator.
- **```Gbest_Show.m```:** Computes the fitness value of the best individual as well as the network configuration it represents.
- **```HUX_crossover.m```:** Perform the half-uniform crossover operator.
- **```Ini_Pop.m```:** Initialises the population of solutions.
- **```Instance.m```:** Contains all the test network instances that can be used in the experimention.
- **```Mutation.m```:** Implements the mutation operator.
- **```RC_Function.m```:** Computes the fitness function of a given solution.
- **```Vinicity_NRC.m```:** Computes the vicinity of a given cell in the network.
- **```neighbour.m```:** Contains the neighbourdhood configuration of each problem benchmark.
- **```xlwrite.m```:** Implements the Xlwrite function to store the results of experimentations (Ideal for Linux clusters).
- **```Jar```:** Contains the ```jar``` needed to execute ```xlwrite.m```.

> [1] Zakaria Abdelmoiz Dahi, Enrique Alba, Gabriel Luque, A takeover time-driven adaptive evolutionary algorithm for mobile user tracking in pre-5G cellular networks, Applied Soft Computing, Volume 116, 2022, 107992, ISSN 1568-4946, [DOI](https://doi.org/10.1016/j.asoc.2021.107992), [URL](https://www.sciencedirect.com/science/article/pii/S1568494621009145).
