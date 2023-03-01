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

%% ---- Local Search 1 : CHC algorithm acting on 175 - 25/sqrt(dimension)^2 of the population --------------
% Selection:  Binary Tournament selection 
% Crossover (when performing the local search):  HUX crossover
% Crossover probability: always equal to 1

%% - local search 2 : it is a heuristic that produces the remaining 25/sqrt(dimension)^2 individual -

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
for ind=[13:20]
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
          combinations_amount = 100;
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
[bestFiT,indbest] = min(fit);
gbest = population(indbest,:);
ALL_FITNESSES = [ALL_FITNESSES min(fit)]; 
ALL_ITTERATIONS = [ALL_ITTERATIONS indiv];
%% -- Incrementing the counter of fitness evaluations by Indiv because we already consumed Indiv fitness evaluations ----------
it = it + indiv;
%% ******************** Optimisation Process ******************************
while it <= (iter - 4)% --------- itteration loop  ------------------------
if LS_count_down ~= LS_period %(Normal Mathematical Oracle)
% ------------ Calculating the remaining Allowable Estimation -------------
reste = iter - it; 
if reste > 0
if (floor(reste/2)) >= indiv 
   individual = indiv;
else
   individual = floor(reste/2);
end
%% --- Initialise the global fitnesses and population ----------------
glob_population  = [];
glob_fitnesses   = [];
% --------------------------------------------------------------------
%% --- Initialise the offsprings population and fitnesses ------------
offsprings = [];
fit_offsprings = [];
evolving_pm_direction = zeros(1,(Dimension +1));
% --------------------------------------------------------------------
%% ---------  Create N new offsprings --------------------------------------------------------------------------------------------------------
%% ---------  No selection Process Is performed Because All the Population Will Undergo Crossover and Mutation -------------------------------       
  %% -------  Generate N/2 couple of chromosome parents by mean of binary tournament ---------------------------
      % ----- Cr contains the index of the winners of the binary tournament --------------------
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
   %% -------  The genetic crossover and mutation  -----------------------------------
   evaluated_offsprings = 0;
   for k=1: ceil(individual/2)
            %% ---- Apply the two-point crossover on the formed couples ----------------
            new_offspring = Crossover(population,Dimension,k,pcr,cr);
            c1 = new_offspring(1,:);
            c2 = new_offspring(2,:);
            %% ---- Apply Bit-Flip Mutation / Evaluations On the resulted Offsprings ---
            new_offspring = [];
            new_offsprings = Mutation(c1,c2,Dimension,pm,evolving_pm_direction,lambda_omega,gamma,network,neighbourhoud,evaluated_offsprings,individual);
            % --------------------------------------------------------------------------
            %% ---- Add the newly produced offsprings to the population of the offsprings ---
            if new_offsprings(2,:) == zeros(1,(Dimension + 1))
            %% --- Add the fitness of the only one newly produced offspring --------------------------
            offsprings = [offsprings ; new_offsprings(1,1:Dimension)];
            %% --- Add the fitness of the newly produced offspring --------------------------
            fit_offsprings = [fit_offsprings new_offsprings(1,(Dimension + 1))];
            %% ---- Increment the number of produced offsprings -----------------------------
            evaluated_offsprings = evaluated_offsprings + 1;
            else
            offsprings = [offsprings ; new_offsprings(1:2,1:Dimension)];
            % -------------------------------------------------------------------------------
            %% --- Add the fitness of the newly produced offspring --------------------------
            fit_offsprings = [fit_offsprings transpose(new_offsprings(1:2,(Dimension + 1)))];
            %% ---- Increment the number of produced offsprings -----------------------------
            evaluated_offsprings = evaluated_offsprings + 2;
            end
            %% --- Update the counting of the fittest strategies ----------------------------
            evolving_pm_direction = evolving_pm_direction + new_offsprings(3,:);
			%% ---- Update the statistics of the affectation rules --------------------------
            current_record = [current_record ; evolving_pm_direction(1:2)];
   end
        %% ---- this step is to keep only the number of individual , it is facultatif because we already have made a test --------
           offsprings = single(offsprings(1:individual,:)); 
           fit_offsprings = fit_offsprings(1:individual);
        % -----------------------------------------------------------------
        %% -- Updating The Mutation Probability ---------------------------------
          [value,index] = max(evolving_pm_direction);
          if index == 1
             %% If increasing Pm Value Give better offsprings
             pm = (pm * lambda_omega) ;
             mutation_evo = [mutation_evo pm];
          end
          if index == 2 
             %% If keeping the same Pm value give better results we decrease it
             pm = (pm * gamma);
             mutation_evo = [mutation_evo pm];
          end
%% ---------- Evaluation of the Produced Offsprings ------------------------------------------------------
it = it + (2*individual);
%% ---------- Test if the period to Apply the mathematical oracle came or no --------------------------
if count_down ~= 0
    count_down = count_down - 1;
end
end
%% ------------- Apply Generational Elitist Selection : elitis(Lambda + Mu)-------------------
glob_population  = [offsprings; population];
glob_fitnesses   = [fit_offsprings fit];
best_fitness = min(glob_fitnesses);
index = Elitist_selection(glob_fitnesses,glob_population,indiv);
population = single(glob_population(index,:));
fit = glob_fitnesses(index);
%% ------------calculate the new best individual --------------------
[bestFit,indbest] = min(fit);
gbest = population(indbest,:);
fitness_best = fit(indbest);
ALL_FITNESSES = [ALL_FITNESSES fitness_best]; 
ALL_ITTERATIONS = [ALL_ITTERATIONS it];
%% -- In case it is not time yet to apply the mathematical Oracle ----
% --------------------- Compute P(t) & Q(t) & X(t) -------------------
if count_down ~= 0
    P_t = length(find (fit <=  (Alpha * fitness_best)));
    P_t = (P_t / indiv);
    Q_t = 1 - P_t;
    Q_t_Period = [Q_t_Period Q_t];
    X_i = Q_t * 100;
    X_i_Period = [X_i_Period X_i];
end
%% ---- In case it is time yet to apply the mathematical Oracle --------
if count_down == 0
   sigma_t = sqrt((1/(2 * Period)) * (sum((X_i_Period).^2)));
   %% ----------- Test if it is time to apply the local search or not --
   if  sigma_t == 0
       if LS_count_down < LS_period
          LS_count_down = LS_count_down + 1; 
       end
   end
   % (
   sigma_values = [sigma_values;sigma_t];
   growth_curve_value = (1 - exp((-(round(it/indiv))^2)/(2*(sigma_t^2))));
   growth_curve = [growth_curve;growth_curve_value];
   %% ----- If the convergence is too fast ------------
   if  sigma_t < lower_threshold
       if gamma > 0.5
           gamma = gamma - decrease_step_gamma;
       end
       if lambda_omega > 1.1
           lambda_omega = lambda_omega -  decrease_step_lambda;
       end
   end
   % ----- If the convergence is too slow ------------
   if sigma_t > upper_threshold
       if gamma < 0.9
           gamma = gamma + increase_step_gamma ;
       end
       if lambda_omega < 1.5 
          lambda_omega = lambda_omega +  increase_step_lambda;
       end
   end
   % ----- Save The value of gamma for post experiment analysis ----
   gamma_values = [gamma_values;gamma];
   lambda_values = [lambda_values;lambda_omega];
   global_Q_t = [global_Q_t;00;transpose(Q_t_Period)];
   count_down = Period;
   Q_t_Period = [];
   X_i_Period = [];
end
end % --------------- End of the mathematical oracle phase -------------



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
            if neighbourhoud_num <  6
                neighbourhoud_num = neighbourhoud_num + 1;
                combinations_num = 0;  
                allowed_pop2 = combinations_amount; 
                if allowed_pop1 > allowed_pop2 
                    allowed_pop1 = allowed_pop1 - allowed_pop2;
                    check_test = 1;
                else
                    check_test = 0;
                end
            else     
                
                check_test = 0;
            end
        end
    else
        check_test = 0;
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
    % ------ Evaluate the produced individuals --------------------------------
    fit = [];
    for w=1:individual
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
        if neighbourhoud_num <  6
            neighbourhoud_num = neighbourhoud_num + 1;
            combinations_num = 0;  
            allowed_pop2 = combinations_amount; 
            if allowed_pop1 > allowed_pop2 
                allowed_pop1 = allowed_pop1 - allowed_pop2;
                check_test = 1;
            else
                check_test = 0;
            end
        else
            check_test = 0;
        end
    end
else
    check_test = 0;
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
%% ----------- Concatenate the produced offsprings -------------
if rows_off == 0
    offsprings = [];
else
    offsprings = [offsprings;pop_1];
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
    %% ----------- Concatenate the produced offsprings ---------------
    offsprings = [offsprings;pop_2];
    total_produced = total_produced + allowed_pop2;
end
%% ----------- Evaluate the newly produced offsprings ----------------------
for ll=1:total_produced
    result  = RC_Function(offsprings(ll,:),Dimension,network,neighbourhoud);
    fit_offsprings = [fit_offsprings , result];
end
%% ----------- Extract firstly the best individual for the offsprings ------
[best_fit_prev,indbest] = min(fit_offsprings);
gbest_prev = offsprings(indbest,:);
fitness_best_prev = fit_offsprings(indbest);
%% ----- Apply Generational Elitist Selection : elitis(Lambda + Mu)---------
if check_test == 0  
    if rows_off == 0
        offsprings = [];
        fit_offsprings = [];
    end
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
    if rows_off == 0
        offsprings = [];
        fit_offsprings = [];
    else
        offsprings = offsprings(1:rows_off,:);
        fit_offsprings = fit_offsprings(1:rows_off);
        glob_population  = [offsprings; population];
        glob_fitnesses   = [fit_offsprings fit]; 
        best_fitness = min(glob_fitnesses); 
        index = Elitist_selection(glob_fitnesses,glob_population,selection_list);
        population = single(glob_population(index,:)); 
        fit = glob_fitnesses(index);
    end
    %% ------------ calculate the new best individual ----------------------
    [bestFit,indbest] = min(fit);
    gbest = population(indbest,:);
    fitness_best = fit(indbest);
    %% ----- Select the best individual from both Pop 1 and Pop 2 ----------
    if fitness_best ~= fitness_best_prev
        if fitness_best_prev < fitness_best        
            population = [gbest_prev; population]; 
            fit = [fitness_best_prev fit]; 
            population = population(1:selection_list,:); 
            fit = fit(1:selection_list);    
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
                if neighbourhoud_num <  6
                    neighbourhoud_num = neighbourhoud_num + 1;
                    combinations_num = 0;  
                    allowed_pop2 = combinations_amount; 
                    if allowed_pop1 > allowed_pop2 
                        allowed_pop1 = allowed_pop1 - allowed_pop2;
                        check_test = 1;
                    else
                        check_test = 0;
                    end
                else
                    check_test = 0;
                end
            end
        else
            check_test = 0;
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
        for w=1:individual
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
end %---------------- end of the loop t=1:iter-------------------------------
if FITNESS_GBEST > fitness_best
FITNESS_GBEST = fitness_best;
GBEST_ALL = gbest;
end
ALL_EXECUTION = [ALL_EXECUTION fitness_best];
oo = toc;
ti = ti + oo;
ALL_TIME = [ALL_TIME ti];
%% ---- Record the mean of the updates performed --------------------------
[rows,cols] = size(current_record); 
rule_1 = mean(current_record(1:rows,1));
rule_2 = mean(current_record(1:rows,2));
rules = [rule_1 rule_2];
record = [record ; rules];
end %% ----------- end of the loop of execution ---------------------------
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
xlwrite(name,record,11);
xlwrite(name,gamma_values,12);
xlwrite(name,lambda_values,13);
xlwrite(name,sigma_values,14);
xlwrite(name,global_Q_t,15);
xlwrite(name,growth_curve,16);
xlwrite(name,LS_loop,17);
xlwrite(name,transpose(fit),18);
end
%% I added this command because Daniel told me that if i don't add it it will ot escape and display the results of the run 
exit;
