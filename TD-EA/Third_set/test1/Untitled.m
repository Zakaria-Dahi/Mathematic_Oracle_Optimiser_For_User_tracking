%% ---- This version is based on variant (rule 1) of task 31. It is our best variant using the mathematical oracle
%% ---- In this code we begin enhancing this version with local search.
%% ---- Basically this variant use the MO as an optimisation process and the CHC as an LS

%% ---- MO ----------------
% Selection:  Binary Tournament selection 
% Size of tournament: Constant (binary tournament)
% Crossover (when performing the mathematical oracle):  Two-point crossover
% Probability of crossover (for both phases): constant always equal 1
% Mutation: Bit Flip Mutation
% Probability of Mutation: deterministic adaptive mutation probability using a modified variant of "declining strategy"
% Info on the affectation rules: It uses two affectation rules :  increasing, decreasing.
% Info on the declining strategy:  It uses a mathematical oracle to drive the value of gamma between two bounds UpperBound = 0.9 to LowerBound = 0.5
% Info on the declining strategy:  It uses a mathematical oracle to drive the value of lambda between two bounds UpperBound = 1.1 to LowerBound =  1.5
% As a Lambda0 and Gamma0 value where we start : it is UB - ((UB - LB)/2)
% Additional Information: It uses a mathematical oracle "Alpha Fitness Best Found"
% Additional Information: It is used to drive the values of Gamma starting from 0.9 till 0.5. and the value of lambda from 1.1 to 1.5
% The value of Alpha is: 1.02 according to previous experiment (see experiment recording)
% N : the application period is set to 30 generatios, according to the authors it is the most important criteria
% As a step of increase decrease we computed: (UB - LB)/(#Evaluations/N) 
% Seed : different seeds of random generator Mersenne Twister 
% Lower and Upper bound for Mathematical Oracle: LB = 60 or 50 and UP = 70 based on the Tasks 25, 26 results
% The local_search for instances under 36 cells started after 2 mathematical oracle periods when
% sigma = 0 and 4 periods for instances above 36 cells
% Memory:  reduction of memory use using the conversion to char type 

%% ---- Local Search 1 : CHC algorithm acting on sqrt(dimension) of the population --------------
% Selection:  Binary Tournament selection 
% Crossover (when performing the local search):  HUX crossover
% Crossover probability: always equal to 1

%% - local search 2 : it is a heuristic that produces the remaining sqrt(Dimension) population -

%% Author : Zakaria Abd El Moiz DAHI
%% University : Constantine 2, Algeria

clear all
%% -------- Open the Pool ---------------------------------------------------------------------------------
% matlabpool open
%% --------------- Initialisation of POI Libs To Write Excel Files ----------------------------------------
% Add Java POI Libs to matlab javapath
javaaddpath('Jar/poi-3.8-20120326.jar');
javaaddpath('Jar/poi-ooxml-3.8-20120326.jar');
javaaddpath('Jar/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('Jar/xmlbeans-2.3.0.jar');
javaaddpath('Jar/dom4j-1.6.1.jar');
javaaddpath('Jar/stax-api-1.0.1.jar');
%% -------------------- Starting the execution of the program ---------------------------------------------
for ind=[7,10,11,12]
%% -------------   Initialize the parameters of the experiments -------------------------------------------
%%%% -------------- The number of execution ------------------------------
global execution 
       execution = [1 1 1 1 1 1 1 1 1 1 1 1];
%%%% -------------- The number of evaluations ----------------------------
global fitness_eval 
       fitness_eval = [1750 1750 1750 1750 1750 1750 1750 1750 1750 1750 1750 1750];
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
%% Prameters  of the Genetic Algorithm  ----------------------
     pcr = 1;  % By default I set it to 1 as recomended by Professor Alba
     pm = 0.5 - ((0.5 - (1/Dimension))/2); % As a starting value we use a value in the middle of adviced intervals in the literature 1/L .. 0.5, I choosed this starting vaue to be fair and give the values the same space to decrease or increase
     lambda_omega = 1.5 - ((1.5 - 1.1)/2); % To subsitute the "lambda" and the "omega" variable by only one variable
     gamma = 0.9 - ((0.9 - 0.5)/2);
%% Parameters of the mathematical Oracle ---------------------
     Alpha = 1.02;
     Period = 30; % This defines the period which we apply the mathematical oracle
     decrease_step_gamma = ((0.9 - 0.5)/2)/ (501/Period);% 501 is the # of itterations the algorithm can perform using this double evaluation mechanism (see # of rows of the experiment excel files)
     increase_step_gamma = ((0.9 - 0.5)/2)/ (501/Period); % The same explanation as above
     decrease_step_lambda = ((1.5 - 1.1)/2)/ (501/Period); % this was done in purpose because based on the previous experiment the value of lambda had no great influence.
     increase_step_lambda = ((1.5 - 1.1)/2)/ (501/Period); % this was done in purpose because based on the previous experiment the value of lambda had no great influence.
     if Dimension <= 36 
	    lower_threshold = 60; % This bound are made based on previous experiment (see experiment 26 recording) 60 is the best for instances of size 4x4 and 6x6 Cells
     else
	    lower_threshold = 50; % This bound are made based on previous experiment (see experiment 26 recording) 60 is the best for instances of size higher than 6x6 cells
	 end
     upper_threshold = 70; % This bound are made based on previous experiment (see experiment recording)
%% Paremetrs of the local search 1 -----------------------------
   global LS_period
   if ind == 11
       LS_period = sqrt(Dimension) + 2;
   else
       LS_period = sqrt(Dimension);
   end
   global catclysmic_rate
   if ind == 11
       catclysmic_rate = 0.6;
   else
       if ind == 10
           catclysmic_rate = 0.75;
       else
           catclysmic_rate = 0.5;
       end
   end
%% Paremetrs of the local search 2 -----------------------------
   global D
          D = round(sqrt(Dimension));
   global combinations
          combinations = str2mat(dec2bin(0:2^D-1))-'0';
   global neighbourhoud_num
          neighbourhoud_num = 1;
   global combinations_amount
          combinations_amount = 25;
          if combinations_amount > (2^(D))
              combinations_amount = 2^(D);
          end
   global combinations_num
          combinations_num = 0;
   %% ----- Extract the First neighbourhoud ------
   seed_1 = 1;
   neighbourhood_1 = (seed_1:seed_1+((sqrt(Dimension)/2) - 1));
   seed_2 = seed_1 + round(sqrt(Dimension));
   neighbourhood_2 = (seed_2:seed_2+((sqrt(Dimension)/2) - 1));
   neighbourdhood_one = [neighbourhood_1 neighbourhood_2];
   %% ----- Extract the Second neighbourhoud -----
   seed_3 = ((round(sqrt(Dimension)/2) - 1) * round(sqrt(Dimension))) + 1;
   neighbourhood_3 = (seed_3:seed_3+((sqrt(Dimension)/2) - 1));
   seed_4 = (round(sqrt(Dimension)/2) * round(sqrt(Dimension))) + 1;
   neighbourhood_4 = (seed_4:seed_4+((sqrt(Dimension)/2) - 1));
   neighbourdhood_two = [neighbourhood_3 neighbourhood_4];
   %% ----- Extract the Third neighbourhoud ------
   seed_5 = ((round(sqrt(Dimension)) - 2) * round(sqrt(Dimension))) + 1;
   neighbourhood_5 =  (seed_5:seed_5+((sqrt(Dimension)/2) - 1));
   seed_6 = round(sqrt(Dimension) - 1) * round(sqrt(Dimension)) + 1;
   neighbourhood_6 =  (seed_6:seed_6+((sqrt(Dimension)/2) - 1));
   neighbourdhood_three = [neighbourhood_5 neighbourhood_6];
   %% ----- Extract the Fourth neighbourhoud -----
   seed_7 = seed_1 + round(sqrt(Dimension)/2);
   neighbourhood_7 = (seed_7:seed_7+((sqrt(Dimension)/2) - 1));
   seed_8 = seed_2 + round(sqrt(Dimension)/2);
   neighbourhood_8 = (seed_8:seed_8+((sqrt(Dimension)/2) - 1));
   neighbourdhood_four = [neighbourhood_7 neighbourhood_8];
   %% ----- Extract the Fifth neighbourhoud -----
   seed_9 = seed_3 + round(sqrt(Dimension)/2);
   neighbourhood_9 = (seed_9:seed_9+((sqrt(Dimension)/2) - 1));
   seed_10 = seed_4 + round(sqrt(Dimension)/2);
   neighbourhood_10 = (seed_10:seed_10+((sqrt(Dimension)/2) - 1));
   neighbourdhood_five = [neighbourhood_9 neighbourhood_10];
   %% ----- Extract the Sixth neighbourhoud -----
   seed_11 = seed_5 + round(sqrt(Dimension)/2);
   neighbourhood_11 = (seed_11:seed_11+((sqrt(Dimension)/2) - 1));
   seed_12 = seed_6 + round(sqrt(Dimension)/2);
   neighbourhood_12 = (seed_12:seed_12+((sqrt(Dimension)/2) - 1));
   neighbourdhood_six = [neighbourhood_11 neighbourhood_12];
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
            FITNESS_GBEST = 1000000000000000000000000000000000000; %% It is made like this to be sure that no intial fitness is higher enough to keep it
     %%%% to save the best individual all along the Nexe exveution 
     global GBEST_ALL
            GBEST_ALL = [];
     %%% In this version we used non-paire value of number of evaluations 
     global it
     %%% To save the recods of the strategies
     global record 
            record = [];
%% ------------------------------------------------------------------------
for exe=1:Nexe
%% ----- Change the seed of the Mersenne Twister Random Generator ---------
rng shuffle
% ------------------------------------------------------------------------6
ti = 0;
tic
%% --- To reset The number Of Fitness Evolution ----------
it = 0;
%% variable used during the calcul ---------------------------------------
     %%% Vector conaiting the fitness of all individuals of the population 
            fit = [];
     %%% Vector containing the best individual so far 
     global gbest
            gbest = [];
     %%%% The population 
            population = [];
     %%% To save the itteration and the corresponding fitness value of this itteration (fitness evauation)
            ALL_FITNESSES = [];
            ALL_ITTERATIONS = [];
     %%% To save the probability of mutation used during the execution -------------------
     mutation_evo = [pm];
     %%% To save the statistical information of the population ----------
     global Stat_Pop
            Stat_Pop = [];
     global current_record
            current_record = [];
    %% To save variables needed to apply the mathematical oracle
    global gamma_values
           gamma_values = [gamma];
    global lambda_values
           lambda_values = [lambda_omega];
    global sigma_values
           sigma_values = [];
    global count_down 
           count_down = Period;
    global Q_t_Period
           Q_t_Period = [];
    global X_i_Period
           X_i_Period = [];
    global global_Q_t 
           global_Q_t = [];
    global growth_curve
           growth_curve = [];
    %% Variable used to rule the locals earch 
    global LS_count_down
           LS_count_down = 0;
    %% Variable used for the HUX crossover 
    global HUX_distance
           HUX_distance = fix(Dimension/4);
    %% Save the local serahc iteration number 
    global LS_loop
           LS_loop = 0;
     gbest = round(rand(1,Dimension));
     fitness_best = RC_Function(gbest,Dimension,network,neighbourhoud);
     LS_count_down = LS_period ;
%% -- Incrementing the counter of fitness evaluations by Indiv because we already consumed Indiv fitness evaluations ----------
it = it + indiv;
%% ******************** Optimisation Process ******************************
while it <= (iter - 4)% --------- itteration loop  ------------------------
 %% ************************ Local Search Process ***********************
%% ---------  Test if we apply the local search or not -----------------
if LS_count_down == LS_period 
%% ----------- Increase the counter of LS-------------------------------
LS_loop = LS_loop + 1;
%% --- Initialise the global fitnesses and population ------------------
glob_population  = [];
glob_fitnesses   = [];
% ----------------------------------------------------------------------
%% --- Initialise the offsprings population and fitnesses --------------
offsprings = [];
fit_offsprings = [];
%% --- Save the previous Gbest -----------------------------------------
gbest_comp = gbest;
fitness_best_comp = fitness_best;
%% ----------------------- Begin Local Search --------------------------
reste = iter - it; 
if reste > 0
if reste >= indiv 
   individual = indiv;
else
   individual = reste;
end
% ----- In case this is the first loop of the Local Search -------------
if LS_loop == 1
    %% ---- Create Pop_1 and Pop_2 --------------------
    allowed_pop1 = individual - 1;
    if neighbourhoud_num <= 6
        test = 2^(round(sqrt(Dimension))) - combinations_num;
        if test > 0
            if test >= combinations_amount
                allowed_pop2 = combinations_amount;
            else
                allowed_pop2 = test;
            end
            if allowed_pop1 > allowed_pop2
                allowed_pop1 = allowed_pop1 - allowed_pop2;
                check_test = 1; 
            else
                check_test = 0;
            end
        else
            neighbourhoud_num = neighbourhoud_num + 1;
            combinations_num = 0;  
            allowed_pop2 = combinations_amount; 
            if allowed_pop1 > allowed_pop2 
                allowed_pop1 = allowed_pop1 - allowed_pop2;
                check_test = 1;
            else
                check_test = 0;
            end
        end
    end
    % ----- Generate the first population of LS 1: pop 1 -----------
    pop_1 = round(rand(allowed_pop1,Dimension));
    population = [gbest;pop_1];
    % ----- Generate the first population of LS 2: pop 2 -----------
    if check_test == 1
        pop_2 = repmat(gbest,allowed_pop2,1); 
        index_1 = combinations_num + 1;
        index_2 = index_1 + (allowed_pop2 - 1);
        if neighbourhoud_num == 1 
            pop_2(:,neighbourdhood_one) = combinations(index_1:index_2,:);
        end
        if neighbourhoud_num == 2 
            pop_2(:,neighbourdhood_two) = combinations(index_1:index_2,:); 
        end
        if neighbourhoud_num == 3 
            pop_2(:,neighbourdhood_three) = combinations(index_1:index_2,:); 
        end
        if neighbourhoud_num == 4 
            pop_2(:,neighbourdhood_four) = combinations(index_1:index_2,:); 
        end
        if neighbourhoud_num == 5 
            pop_2(:,neighbourdhood_five) = combinations(index_1:index_2,:); 
        end
        if neighbourhoud_num == 6 
            pop_2(:,neighbourdhood_six) = combinations(index_1:index_2,:); 
        end
        % --------- Update the ID of the combination reached --------------
         combinations_num = combinations_num + allowed_pop2;
         % ---- concatenate both population of LS1 and LS2: pop 1 and pop 2 ----
         population = [population;pop_2];
    end
    % ------ Evaluate the produced individuals --------------
    fit = [];
    parfor w=1:individual
        result  = RC_Function(population(w,:),Dimension,network,neighbourhoud);
        fit = [fit , result];
    end
    %% ------------ calculate the new best individual --------------------
    if check_test == 0  
        [best_fit,indbest] = min(fit);
        gbest = population(indbest,:);
        fitness_best = fit(indbest);
        if fitness_best_comp ~= fitness_best 
            fitness_best_comp = fitness_best;
            gbest_comp = gbest;
            neighbourhoud_num = 1;
            combinations_num = 0;             
        end
    else
        [best_fit,indbest] = min(fit);
        gbest_prev = population(indbest,:);
        fitness_best_prev = fit(indbest);
        %% ------------ Keep only the first population -------------------
        allowed_pop1 = allowed_pop1 + 1;
        population = population(1:allowed_pop1,:); 
        fit = fit(1:allowed_pop1);  
        [best_fit,indbest] = min(fit); 
        gbest = population(indbest,:); 
        fitness_best = fit(indbest);
        %% ----- Select the best individual from both Pop 1 and Pop 2 ----
        if fitness_best ~= fitness_best_prev
            if fitness_best_prev < fitness_best                  
                population = [gbest_prev; population]; 
                fit = [fitness_best_prev fit]; 
                population = population(1:allowed_pop1,:); 
                fit = fit(1:allowed_pop1); 
                gbest = gbest_prev;
                fitness_best = fitness_best_prev;
            else
                [worst_fit,ind_worst] = max(fit);
                if fitness_best_prev < worst_fit
                    population(ind_worst,:) = gbest_prev; 
                    fit(ind_worst) = fitness_best_prev; 
                end
            end
        end
        if fitness_best_comp ~= fitness_best
            fitness_best_comp = fitness_best;
            gbest_comp = gbest;
            neighbourhoud_num = 1;
            combinations_num = 0;             
        end
    end
    % ----- Increase the NOFE ----------------
    it = it + individual;
end
% ----- In case this is not the first loop of the Local Search -----------
reste = iter - it; 
if reste > 0
if reste >= indiv 
   individual = indiv;
else
   individual = reste;
end
% ----- save the previous population ---------------------
population_save = population;
%% ---- Create Pop_1 and Pop_2 ---------------------------
allowed_pop1 = individual;
if neighbourhoud_num <= 6
    test = 2^(round(sqrt(Dimension))) - combinations_num;
    if test > 0
        if test >= combinations_amount
            allowed_pop2 = combinations_amount;
        else
            allowed_pop2 = test;
        end
        if allowed_pop1 > allowed_pop2
            allowed_pop1 = allowed_pop1 - allowed_pop2;
            check_test = 1; 
        else
            check_test = 0;
        end
    else
        neighbourhoud_num = neighbourhoud_num + 1;
        combinations_num = 0;  
        allowed_pop2 = combinations_amount; 
        if allowed_pop1 > allowed_pop2 
            allowed_pop1 = allowed_pop1 - allowed_pop2;
            check_test = 1;
        else
            check_test = 0;
        end
    end
end
% ----- Perform the binary tournament --------------------
Pre_HUX = [];
selection_list = indiv - combinations_amount;
for k=1:ceil(allowed_pop1/2)
    % --- Randomly choose parent 1 and 2 for the first tournament ------
    p_one = randi(selection_list);
    p_two = randi(selection_list);
    % --- Randomly choose parent 3 and 4 for the second tournament ------
    p_three = randi(selection_list);
    p_four  = randi(selection_list);
    % ---- First binary tournament -----
    if fit(p_one) < fit(p_two)
        Pre_HUX  =  [Pre_HUX;population(p_one,:)]; 
    else
        Pre_HUX  =  [Pre_HUX;population(p_two,:)];
    end
    % ---- Second binary tournament -----
    if fit(p_three) < fit(p_four)
        Pre_HUX  =  [Pre_HUX;population(p_three,:)];
    else
        Pre_HUX  =  [Pre_HUX;population(p_four,:)];
    end
end
%% ---------------------- Apply the HUX crooosver --------------------------
pop_1 = HUX_crossover(Pre_HUX,allowed_pop1,HUX_distance,pcr);
[rows_off,cols_off] = size(pop_1);
if rows_off > allowed_pop1
    rows_off = allowed_pop1;
    pop_1 = pop_1(1:allowed_pop1,:);
end
total_produced = rows_off;
% ----- Generate the first population of LS 2: pop 2 -----------
if check_test == 1 
    pop_2 = repmat(gbest,allowed_pop2,1); 
    index_1 = combinations_num + 1; 
    index_2 = index_1 + (allowed_pop2 - 1); 
    if neighbourhoud_num == 1 
        pop_2(:,neighbourdhood_one) = combinations(index_1:index_2,:);
    end
    if neighbourhoud_num == 2 
        pop_2(:,neighbourdhood_two) = combinations(index_1:index_2,:); 
    end
    if neighbourhoud_num == 3 
        pop_2(:,neighbourdhood_three) = combinations(index_1:index_2,:); 
    end
    if neighbourhoud_num == 4 
        pop_2(:,neighbourdhood_four) = combinations(index_1:index_2,:); 
    end
    if neighbourhoud_num == 5 
        pop_2(:,neighbourdhood_five) = combinations(index_1:index_2,:); 
    end
    if neighbourhoud_num == 6 
        pop_2(:,neighbourdhood_six) = combinations(index_1:index_2,:); 
    end
    % --------- Update the ID of the combination reached --------------
    combinations_num = combinations_num + allowed_pop2;
    %% ----------- Concatenate both produced offsprings ---------------
    offsprings = [pop_1;pop_2];
    total_produced = total_produced + allowed_pop2;
end
%% ----------- Evaluate the newly produced offsprings ----------------------
parfor ll=1:total_produced
    result  = RC_Function(offsprings(ll,:),Dimension,network,neighbourhoud);
    fit_offsprings = [fit_offsprings , result];
end
%% ----------- Extract firstly the best individual for the offsprings ------
[best_fit_prev,indbest] = min(fit_offsprings);
gbest_prev = offsprings(indbest,:);
fitness_best_prev = fit_offsprings(indbest);
%% ----- Apply Generational Elitist Selection : elitis(Lambda + Mu)---------
if check_test == 0  
    glob_population  = [offsprings; population];
    glob_fitnesses   = [fit_offsprings fit];
    best_fitness = min(glob_fitnesses);
    index = Elitist_selection(glob_fitnesses,glob_population,selection_list);
    population = single(glob_population(index,:));
    fit = glob_fitnesses(index);   
    %% ------------ calculate the new best individual ----------------------
    [best_fit,indbest] = min(fit);
    gbest = population(indbest,:);
    fitness_best = fit(indbest);
    %% ----------- Update the neighbourhood --------------------------------
    if fitness_best_comp ~= fitness_best 
        fitness_best_comp = fitness_best;
        gbest_comp = gbest; 
        neighbourhoud_num = 1; 
        combinations_num = 0;         
    end
else
    offsprings = offsprings(1:rows_off,:);
    fit_offsprings = fit_offsprings(1:rows_off);
    glob_population  = [offsprings; population]; 
    glob_fitnesses   = [fit_offsprings fit]; 
    best_fitness = min(glob_fitnesses); 
    index = Elitist_selection(glob_fitnesses,glob_population,selection_list);
    population = single(glob_population(index,:)); 
    fit = glob_fitnesses(index);
    %% ------------ calculate the new best individual ----------------------
    [bestFit,indbest] = min(fit);
    gbest = population(indbest,:);
    fitness_best = fit(indbest);
    %% ----- Select the best individual from both Pop 1 and Pop 2 ----------
    if fitness_best ~= fitness_best_prev
        if fitness_best_prev < fitness_best        
            population = [gbest_prev; population]; 
            fit = [fitness_best_prev fit]; 
            population = population(1:allowed_pop1,:); 
            fit = fit(1:allowed_pop1);    
            % -------- Update the best individual --------
            gbest = gbest_prev;
            fitness_best = fitness_best_prev;
        else
            [worst_fit,ind_worst] = max(fit);
            if fitness_best_prev < worst_fit
                population(ind_worst,:) = gbest_prev; 
                fit(ind_worst) = fitness_best_prev; 
            end
        end           
    end
    %% ----------- Update the neighbourhood --------------------------------
    if fitness_best_comp ~= fitness_best
        fitness_best_comp = fitness_best;
        gbest_comp = gbest; 
        neighbourhoud_num = 1; 
        combinations_num = 0;   
    end 
end
ALL_FITNESSES = [ALL_FITNESSES fitness_best]; 
ALL_ITTERATIONS = [ALL_ITTERATIONS it];
% -------------------- Increase the NOFE -----------------------------------
it = it + total_produced;
%% - Test if the current population and the previous one are identical -----
test_answer = zeros(1,selection_list);
for i=1:selection_list 
    for j=1:selection_list
        answer = isequal(population(i,:),population_save(j,:)); 
        if answer == 1
            test_answer(i) = 1;
            break;
        end
    end
end
if isequal(test_answer,ones(1,selection_list)) == 1
   HUX_distance = HUX_distance - 1;
end
% - if the threshld is equal to zero we engage the cataclysmic restarting --
if HUX_distance == 0
    reste = iter - it; 
    if reste > 0
        if reste >= indiv 
            individual = indiv;
        else
            individual = reste;
        end
        allowed_pop1 = individual - 1;
        if neighbourhoud_num <= 6
            test = 2^(round(sqrt(Dimension))) - combinations_num;
            if test > 0
                if test >= combinations_amount
                    allowed_pop2 = combinations_amount;
                else
                    allowed_pop2 = test;
                end
                if allowed_pop1 > allowed_pop2
                    allowed_pop1 = allowed_pop1 - allowed_pop2;
                    check_test = 1; 
                else
                    check_test = 0;
                end
            else
                neighbourhoud_num = neighbourhoud_num + 1;
                combinations_num = 0;  
                allowed_pop2 = combinations_amount; 
                if allowed_pop1 > allowed_pop2 
                    allowed_pop1 = allowed_pop1 - allowed_pop2;
                    check_test = 1;
                else
                    check_test = 0;
                end
            end
        end
        % ---------- Create the new population by cataclysmic generation ---
        pop_1 =  repmat(gbest,allowed_pop1,1);
        allowed_cata_bits = floor(catclysmic_rate * Dimension);
        for s=1:allowed_pop1
            cata_bits = randi([1,Dimension],1,allowed_cata_bits);
            pop_1(s,cata_bits) = abs(pop_1(s,cata_bits) - 1);
        end
        population = [gbest;pop_1];
        % ----- Generate the first population of LS 2: pop 2----------------
        if check_test == 1 
            pop_2 = repmat(gbest,allowed_pop2,1); 
            index_1 = combinations_num + 1; 
            index_2 = index_1 + (allowed_pop2 - 1); 
            if neighbourhoud_num == 1 
                pop_2(:,neighbourdhood_one) = combinations(index_1:index_2,:);
            end
            if neighbourhoud_num == 2
                pop_2(:,neighbourdhood_two) = combinations(index_1:index_2,:); 
            end
            if neighbourhoud_num == 3
                pop_2(:,neighbourdhood_three) = combinations(index_1:index_2,:);              
            end
            if neighbourhoud_num == 4
                pop_2(:,neighbourdhood_four) = combinations(index_1:index_2,:); 
            end
            if neighbourhoud_num == 5
                pop_2(:,neighbourdhood_five) = combinations(index_1:index_2,:); 
            end
            if neighbourhoud_num == 6
                pop_2(:,neighbourdhood_six) = combinations(index_1:index_2,:); 
            end
            % --------- Update the ID of the combination reached -----------
            combinations_num = combinations_num + allowed_pop2;
            %% ----------- Concatenate both produced offsprings ------------
            population = [population;pop_2];
        end
        % ------ Evaluate the produced individuals -------------------------
        fit = [];
        parfor w=1:individual
            result  = RC_Function(population(w,:),Dimension,network,neighbourhoud);
            fit = [fit , result];
        end
        if check_test == 0     
            %% ------------ calculate the new best individual --------------
            [best_fit,indbest] = min(fit);
            gbest = population(indbest,:);
            fitness_best = fit(indbest);
            %% ----------- Update the neighbourhood ------------------------
            if fitness_best_comp ~= fitness_best 
                fitness_best_comp = fitness_best;
                gbest_comp = gbest; 
                neighbourhoud_num = 1; 
                combinations_num = 0;      
            end
        else
            % ------------calculate the new best individual ----------------
            [best_fit,indbest] = min(fit);
            gbest_prev = population(indbest,:);
            fitness_best_prev = fit(indbest);
            %% ------------ Keep only the first population -----------------
            allowed_pop1 = allowed_pop1 + 1;
            population = population(1:allowed_pop1,:); 
            fit = fit(1:allowed_pop1);  
            [best_fit,indbest] = min(fit); 
            gbest = population(indbest,:); 
            fitness_best = fit(indbest);
            %% ----- Select the best individual from both Pop 1 and Pop 2 --
            if fitness_best ~= fitness_best_prev
                if fitness_best_prev < fitness_best  
                    population = [gbest_prev; population]; 
                    fit = [fitness_best_prev fit]; 
                    population = population(1:allowed_pop1,:); 
                    fit = fit(1:allowed_pop1);  
                    % -------- Update the best individual --------
                    gbest = gbest_prev;
                    fitness_best = fitness_best_prev;
                else
                    [worst_fit,ind_worst] = max(fit);
                    if fitness_best_prev < worst_fit
                        population(ind_worst,:) = gbest_prev;                        
                        fit(ind_worst) = fitness_best_prev; 
                    end
                end
            end
            %% ----------- Update the neighbourhood ------------------------
            if fitness_best_comp ~= fitness_best 
                fitness_best_comp = fitness_best;
                gbest_comp = gbest; 
                neighbourhoud_num = 1; 
                combinations_num = 0;      
            end
        end
        % ------- Reset the threshold of the HUX crossover -----------------
        HUX_distance = fix(catclysmic_rate * (1 - catclysmic_rate) * Dimension);
        % -------------------- Increase the NOFE ---------------------------
        it = it + individual;
    end %-- see if there still enough fitness evauations to perform Catclysmic Restart ---- 
end %-- Condition to perform the cataclysmic Restart ----------------------------------
end % --------------- see if there still enough fitness evauations to perform LS 1 ---
end % --------------- see if there still enough fitness evauations to perform LS  2 --
end % --------------- end of the local search -------------------------------   
size(population)
size(fit)
end
end
end
%% ---- Close the matlabpool ------------------------------------
matlabpool('close');
%% I added this command because Daniel told me that if i don't add it it will ot escape and display the results of the run 
