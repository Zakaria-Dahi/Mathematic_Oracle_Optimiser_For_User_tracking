# Quick Description

**Programmer:shipit:**: Zakaria Abdelmoiz DAHI, University of Malaga, Spain. 

**About:** This repositiory implements the strcutred adaptive evolution algorithm devised and studied in [1] that solves the users' mobility tracking in cellular networks.

>  [1] **Z.A. DAHI**, E. Alba, G. Luque, A takeover time-driven adaptive evolutionary algorithm for mobile user tracking in pre-5G cellular networks, Applied Soft Computing, Volume 116, 2022, 107992, ISSN 1568-4946, https://doi.org/10.1016/j.asoc.2021.107992.

## **How :green_book:** 

- Depending on what you want to study or used, you can navigate to the folder `TD-EA` for running the devised algorithm or `SOTA` for running state-of-the-art algorithms.
- Tehn, you just need to run the Matlab script `main.m`.

## **Folders Hiearchy :open_file_folder:**
    
- `TD-EA`: contains the code of the devised approach.
  - `Firs_Set`: uses 25 combination each time.
  - `Second_Set`: uses 50 combination each time..
  - `Third Set`: uses 100 combination each time. 
- `SOTA`: contains the code of the state-of-the-art algorithms.
  - `DE`: implements the differential evolution algorithm.
  - `GPSO`: implements the the geometric particle swarm optimisation algorithm.
  - `OIGA`: implements the oscillatory increasing genetic algorithm.
- `Results`: Once executed the results are stored in excel files with name Network_a_bxc, where a represent the ID of networks of such shape, b and c represent the number of cells in the width and height of the network.
    
## **Demo :movie_camera:**
- Please refer to the original paper [HERE](https://www.sciencedirect.com/science/article/pii/S1568494621009145) for more detailed results and discussions.

## **Files:**
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


## **Notes:**
- It is prefered that you use Matlab R2014a to run the code as it is. 
- If you use newer versions of Matlab, you will probably need to update some obselete functions such as ```Mat2Cell```.