%% This version is the same as the Version 4.2 of Genertional Elitist Genetic Algorithm For Reporting Cell Problem (RCP) 
% Selection :  Binary Tournament selection 
% Size of tournament : Constant (binary tournament)
% Crossover :  Two-Point crossover
% Probability of crossover : constnt always equal 1
% Mutation : Bit Flip Mutation
% Probability of Mutation : deterministic adaptive mutation probability using the sinuosidal increasing amplitude  : 0.3 - (1/L)
% Seed : different seeds of random generator Mersenne Twister 
% Memory :  reduction of memory using the conversion to char type 


%% Author : Zakaria Abd El Moiz DAHI
%% University : Constantine 2, Algeria
%% Vesrio 5.1
clear all
%% --------------- Initialisation of POI Libs To Write Excel Files ----------------------------------------
% Add Java POI Libs to matlab javapath
javaaddpath('Jar/poi-3.8-20120326.jar');
javaaddpath('Jar/poi-ooxml-3.8-20120326.jar');
javaaddpath('Jar/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('Jar/xmlbeans-2.3.0.jar');
javaaddpath('Jar/dom4j-1.6.1.jar');
javaaddpath('Jar/stax-api-1.0.1.jar');
%% -------------------- Starting the execution of the program ---------------------------------------------
for ind=[26,27]
%% -------------   Initialize the parameters of the experiments -------------------------------------------
%%%% -------------- The number of execution ------------------------------
global execution 
       execution = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%%%% -------------- The number of evaluations ----------------------------
global fitness_eval 
       fitness_eval = [175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000 175000];
%% ------ Number of individuals contained in the population --------------
global  indiv
        indiv = 175;
%% Save the dimension of the individual (it corresponds to the size of the network)
global Dimension
       Network_size = size(Instance(ind));
       Dimension = Network_size(1);
%% Number of itterations 
global iter
       iter = fitness_eval(ind);
%% Number of execution 
global Nexe
       Nexe = execution(ind);
%% Instances and neighbourhood of the network
global network 
       network = Instance(ind);
global neighbourhoud
       if ind == 1 || ind == 2 || ind == 3 || ind == 13 || ind == 16 ||  ind == 17
           neighbourhoud = neighbour(4);
       end
       if ind == 4 || ind == 5 || ind == 6 || ind == 14
           neighbourhoud = neighbour(6);
       end
       if ind == 7 || ind == 8 || ind == 9 || ind == 15
           neighbourhoud = neighbour(8);
       end
       if ind == 10 || ind == 11 || ind == 12
           neighbourhoud = neighbour(10);
       end
       if ind == 18 
          neighbourhoud =   neighbour(19);
       end
       if ind == 19
          neighbourhoud =   neighbour(63);
       end
       if ind == 20 
          neighbourhoud =   neighbour(99);
       end
       if ind == 21 
          neighbourhoud =   neighbour(144);
       end
       if ind == 22 
          neighbourhoud =   neighbour(196);
       end
       if ind == 23 
          neighbourhoud =   neighbour(256);
       end
       if ind == 24
          neighbourhoud =   neighbour(324);
       end
       if ind == 25
          neighbourhoud =   neighbour(400);
       end
       if ind == 26 
          neighbourhoud =   neighbour(900);
       end
       if ind == 27
          neighbourhoud =   neighbour(900);
       end
%% ---- Saving the type and the dimension of the network used -------------
global instance 
      if ind == 1
          instance =  'Network_1_4x4';
      end
      if ind == 2
          instance =  'Network_2_4x4';       
      end
      if ind == 3
          instance =  'Network_3_4x4';         
      end
      if ind == 4
          instance =  'Network_1_6x6';         
      end
      if ind == 5
          instance =  'Network_2_6x6';         
      end
      if ind == 6
          instance =  'Network_3_6x6';         
      end
      if ind == 7
          instance =  'Network_1_8x8';         
      end
      if ind == 8
          instance =  'Network_2_8x8';
      end
      if ind == 9
          instance =  'Network_3_8x8';          
      end
      if ind == 10
          instance =  'Network_1_10x10';         
      end
      if ind == 11
          instance =  'Network_2_10x10';          
      end      
      if ind == 12
          instance =  'Network_3_10x10';
      end
      if ind == 13
          instance =  'Network_4_4x4';         
      end
      if ind == 14
          instance =  'Network_4_6x6';          
      end      
      if ind == 15
          instance =  'Network_4_8x8';
      end
      if ind == 16
          instance =  'Network_5_4x4';
      end
      if ind == 17
          instance =  'Network_6_4x4';
      end
      if ind == 18
          instance =  'Network_19_Cells';
      end
      if ind == 19
          instance =  'Network_7x9_Cells';
      end
      if ind == 20
          instance =  'Network_9x11_Cells';
      end
      if ind == 21
          instance =  'Network_12x12_Cells';
      end
      if ind == 22
          instance =  'Network_14x14_Cells';
      end
      if ind == 23
          instance =  'Network_16x16_Cells';
      end
      if ind == 24
          instance =  'Network_18x18_Cells';
      end
      if ind == 25
          instance =  'Network_20x20_Cells';
      end
      if ind == 26
          instance =  'Network_1_30x30_Cells';
      end
      if ind == 27
          instance =  'Network_2_30x30_Cells';
      end
%% Prameters  of the Genetic Algorithm  ----------------------
     %%% to save the crossover and mutation probability
     pcr = 1;  % By default I set it to 1 as recomended by Professor Alba
     pm_orig =  0.3 - ((0.3 - (1/Dimension))/2);  % split the distance into hald between 1/L and 0.3 
     freq = 0.7; %% Frequency of the wave
%% variable to save for each execution  ----------------------
     %%%%   Fitness of the best individual reach after the end of all the ietteration 
     global ALL_EXECUTION
            ALL_EXECUTION = [];
     %%%% Save all the execution Time per run 
     global ALL_TIME
            ALL_TIME = [];
     %%%%% save all the best fitnesses through the generations
     global ALL_FITNESSES
            ALL_FITNESSES = [];
     %%%%% save all the itterations where the fitness value where extracted
     global ALL_ITTERATIONS
            ALL_ITTERATIONS = [];
     %%%% Mean of all the fitness reached  
     global moy
     moy =0;
     %%%% Standard deviation of the fitness reached 
     global ecart 
            ecart = 0;
     %%%% to save all the information relative to the experiment 
     global vect_one 
            vect_one = [];
     global vect_two
            vect_two = [];
     %%%% to save to best individual all along the NEXE execution 
     global FITNESS_GBEST
            FITNESS_GBEST = 1000000000000000000000000000000000000;
     %%%% to save the best individual all along the Nexe exveution 
     global GBEST_ALL
            GBEST_ALL = [];
     %%% In this version we used non-paire value of number of evaluations 
     global it
%% ------------------------------------------------------------------------
for exe=1:Nexe
%% ----- Change the seed of the Mersenne Twister Random Generator ----------
rng shuffle
% ------------------------------------------------------------------------
ti = 0;
tic
%% --- To reset The number Of Fitness Evolution ----------
it = 0;
%% variable used during the calcul ---------------------------------------
     %%% Vector conaiting the fitness of all individuals of the population 
     global fit
            fit = [];
     %%% Vector containing the best individual so far 
     global gbest
            gbest = [];
     %%% contain the index of the best individual so far 
     global indbest_one 
     %%%% The population 
     global population
            population = [];
     %%% To save the itteration and the corresponding fitness value of this itteration (fitness evauation)
     ALL_FITNESSES = [];
     ALL_ITTERATIONS = [];
     %%% To save the probability of mutation used during the execution -------------------
     mutation_evo = [pm_orig];
%% ------------ Generate the intial population ---------------------
population = single(Ini_Pop(indiv,Dimension));
% ------------------------------------------------------------------
%% ------ Compute the initial fitnesses for the initial population -------
for w=1:indiv
    result  = RC_Function(population(w,:),Dimension,network,neighbourhoud);
    fit = [fit , result];
end
% ------------------------------------------------------------------------
%% -- calculate the best individual (it's index , it's fitness and ..) ---
[bestFit_one,indbest_one] = min(fit);
gbest = population(indbest_one,:);
ALL_FITNESSES = [ALL_FITNESSES min(fit)]; 
ALL_ITTERATIONS = [ALL_ITTERATIONS indiv];
%% -------------------------------------------------------------------
%% --------- Incrementing the counter of fitness evaluations by Indiv because we already consumed Indiv fitness evaluations ----------
it = it + indiv;
% -------------------------------------------------------------------------
while it <= (iter - 1)% --------- itteration loop  ------------------------
%% -- Increasing sinuosidal wave for the Mutation Probability evolution ---------------------------------
pm = pm_orig * (sind(2 * pi * freq * (it/indiv)) * ((it/indiv)/(iter/indiv)) + 1);
mutation_evo = [mutation_evo pm];
% -------------------------------------------------------------------------
reste = iter - it; 
if reste > 0
if reste >= indiv 
   individual = indiv;
else
   individual = reste;
end
%% --- Initialise the global fitnesses and population ----------------
   glob_population  = [];
   glob_fitnesses   = [];
% --------------------------------------------------------------------
%% --- Initialise the offsprings population and fitnesses ------------
   offsprings = [];
   fit_offsprings = single(zeros(1,individual));
% --------------------------------------------------------------------
%% ---------  Create N new offsprings --------------------------------------------------------------------------------------------------------
%% ---------  No selection Process Is performed Because All the Population Will Undergo Crossover and Mutation -------------------------------       
  %% -----  Generate N/2 couple of chromosome parents by mean of binary tournament ---------------------------
      % ----- Cr contains the indix of the winners of the binary tournament --------------------
       cr =  zeros(ceil(individual/2),2);  
      % ----- This loop is used to perform the binary tournament --------------------
       for k=1:ceil(individual/2)
           % --- Randomly choose parent 1 and 2 for the first tournament ------
           p_one = randi(indiv);
           p_two = randi(indiv);
           % --- Randomly choose parent 3 and 4 for the second tournament ------
           p_three = randi(indiv);
           p_four  = randi(indiv);
           % ---- First binary tournament -----
           if fit(p_one) < fit(p_two)
              cr(k,1)  =  p_one; 
           else
              cr(k,1) =   p_two;
           end
           % ---- Second binary tournament -----
           if fit(p_three) < fit(p_four)
              cr(k,2)  =  p_three; 
           else
              cr(k,2) =   p_four;
           end
       end
   %% -------  The genetic crossover and genetic mutation  -----------------------------------
        for k=1:ceil(individual/2)
            %% ---- Apply the two-point crossover on the formed couples ----------------
            new_offspring = Crossover(population,Dimension,k,pcr,cr);
            c1 = new_offspring(1,:);
            c2 = new_offspring(2,:);
            % --------------------------------------------------------------------------
            %% ---- Apply Bit-Flip Mutation On the resulted Offsprings -----------------
            new_offsprings = Mutation(c1,c2,Dimension,pm);
            % --------------------------------------------------------------------------
            %% ---- Add the newly produced offsprings to the population of the ofsprings ---
            offsprings = [offsprings ; new_offsprings];
            % ------------------------------------------------------------------------------
        end
        %% ---- this step is to keep only the number of individual --------
        offsprings = single(offsprings(1:individual,:)); 
        % -----------------------------------------------------------------
%% ---------- Evaluation of the Produced Offsprings ------------------------------------------------------
for  w=1:individual% ------------ loop of flower processing  global or local pollination ---------------
              % -------- This part of the code is to avoid having individual with no reporting cells ------
              if offsprings(w,:) == zeros(1,Dimension)
                 %% ---- this line of code create at least one reporting cell in random position -----------
                 offsprings(w,randi([1,Dimension])) =  1; 
                 % -----------------------------------------------------------------------------------------
              end
		         %% ---------------------- Evaluate the new solutions---------------------------------------
                  fit_offsprings(w)  = RC_Function(offsprings(w,:),Dimension,network,neighbourhoud);
                 % -----------------------------------------------------------------------------------------
end
it = it + individual;
end
%% ------------- Apply Generational Elitist Selection : elitis(Lambda + Mu)-------------------
glob_population  = [offsprings; population];
glob_fitnesses   = [fit_offsprings fit];
best_fitness = min(glob_fitnesses);
index = Elitist_selection(glob_fitnesses,glob_population,indiv);
population = single(glob_population(index,:));
fit = glob_fitnesses(index);
% -------------------------------------------------------------------------------------------
%% ------------calculate the new best individual --------------------
[bestFit,indbest] = min(fit);
gbest = population(indbest,:);
fitness_best = fit(indbest);
%% --------------------------------------------------------------------
ALL_FITNESSES = [ALL_FITNESSES fitness_best]; 
ALL_ITTERATIONS = [ALL_ITTERATIONS it];
it
% -------------------------------------------------------------------------
end %--------------end of the loop t=1:iter-------------------------------
if FITNESS_GBEST > fitness_best
FITNESS_GBEST = fitness_best;
GBEST_ALL = gbest;
end
ALL_EXECUTION = [ALL_EXECUTION fitness_best];
exe
oo = toc;
ti = ti + oo;
ALL_TIME = [ALL_TIME ti];
end
%% ----  Writing the result -----------------------------------------------
moy =  mean(ALL_EXECUTION);
ecart = std (ALL_EXECUTION);
datte = mat2cell(date);
timee =  num2str(sum(ALL_TIME));
timme =  mat2cell(timee);
%% ---- Writing on Excel File ---------------------------------------------------------------------------------------------
vect_one = [datte Nexe it indiv Dimension instance timme mean(ALL_TIME) std(ALL_TIME)  min(ALL_EXECUTION) max(ALL_EXECUTION) moy ecart];
headers = {'date', 'Nbr_Execution','Nbr_Fitness_Evaluations','Nbr_Individual','Dimension (cells)','Instance','Execution_Time','Mean','Std','Best','Worst','Mean','Std'};
name = strcat(instance,'.xls');
xlwrite(name,[headers;vect_one],1);
%% --- Recording the results obtained by the best individual  -------------------------------
datee =  num2str(date);
result = Gbest_Show(GBEST_ALL,Dimension,network,neighbourhoud);
vect_two = [datte result(1) result(2) result(3)];
headers = {'date','Fitness','Update Location Cost','Paging Cost'};
xlwrite(name,[headers;vect_two],2);
%% ---- Recording The Reporting Cells ID ------------------------------------------------------
vect_two = [cell2mat(result(4))];
headers = {'Reporting Cell'};
xlwrite(name,vect_two',3);
%% ---- Recording the ID of Non Reporting Cells -----------------------------------------------
vect_two = [cell2mat(result(5))];
headers = {'Non Reporting Cell'};
xlwrite(name,vect_two',4);
xlwrite(name,ALL_EXECUTION,6);
xlwrite(name,ALL_TIME,7);
xlwrite(name,transpose(ALL_FITNESSES),8);
xlwrite(name,transpose(ALL_ITTERATIONS),9);
xlwrite(name,transpose(mutation_evo),10);
end
%% I added this command because Daniel told me that if i don't add it it will ot escape and display the results of the run 
exit;
