function [Offspring_resulted] = Mutation(c1,c2,Dimension,pm,evolving_pm_direction,vv,segma,network,neighbourhoud,evaluated_offsprings,individual)
if evaluated_offsprings < individual
if (individual - evaluated_offsprings) == 1
%% --- Initialise the matrix that will contain the 2 new produced best offsprings ---------
Offspring_resulted = [];
%% ---- Initialise the fitness vector containing the fitnesses of the offsprings ----------
fit_offsprings_1 = zeros(1,2);
%% ---- Initialise the population were the produced offsprings are stocked -----------------
offsprings_1 = [];
%% ---- Intalise the resulting offsprings (3 for each offspring) ---
c11 = c1;
c12 = c1;
%% ---- Apply bit-flip to the first variant of offsprings one ------
for n=1:Dimension
    if rand <= (pm * vv)
        c11(n) = 1 -  c11(n);
    end
end
% ------------------------------------------------------------------
%% ---- Apply bit-flip to the second variant of offsprings one -----
for n=1:Dimension
    if rand <= pm
        c12(n) = 1 -  c12(n);
    end
end
% ------------------------------------------------------------------
%% ---- Creating an intermediate population to stock the produced ofssprings ---------------------------
offsprings_1 = [c11;c12];
% ------------------------------------------------------------------------------------------------------
%% --- Evaluating the resulting offsprings (4) = (2 New produced by the mutation)
parfor  off=1:2
		         %% ---------------------- Evaluate the new solutions---------------------------------------
                  fit_offsprings_1(off)  = RC_Function(offsprings_1(off,:),Dimension,network,neighbourhoud);
                 % -----------------------------------------------------------------------------------------
end
%% ---- Extract The best Behaviour Of the New Produced Offsprings ---
[fit_1,index] = min(fit_offsprings_1);
best1 = [offsprings_1(index,:) min(fit_offsprings_1)];
Offspring_resulted = [Offspring_resulted;best1];
evolving_pm_direction(index) = evolving_pm_direction(index) + 1;
%% ---- Add A facultatif Offspring 2 with 0 as soltuion representation and 0 as fitness and do not update the statistica about the Pm strategy value ---- ---
Offspring_resulted = [Offspring_resulted;zeros(1,(Dimension + 1));evolving_pm_direction];
end    

    
    
    
    
    
if (individual - evaluated_offsprings) >= 2
%% --- Initialise the matrix that will contain the 2 new produced best offsprings ---------
Offspring_resulted = [];
%% ---- Initialise the fitness vector containing the fitnesses of the offsprings ----------
fit_offsprings_1 = zeros(1,2);
fit_offsprings_2 = zeros(1,2);
%% ---- Initialise the population were the produced offsprings are stocked -----------------
offsprings_1 = [];
offsprings_2 = [];
%% ---- Intalise the resulting offsprings (2 for each offspring) ---
c11 = c1;
c12 = c1;

c21 = c2;
c22 = c2;
%% ---- Apply bit-flip to the first variant of offsprings one ------
for n=1:Dimension
    if rand <= (pm*vv)
        c11(n) = 1 -  c11(n);
    end
end
% ------------------------------------------------------------------
%% ---- Apply bit-flip to the second variant of offsprings one -----
for n=1:Dimension
    if rand <= pm
        c12(n) = 1 -  c12(n);
    end
end
% ------------------------------------------------------------------
%% ---- Creating an intermediate population to stock the produced ofssprings ---------------------------
offsprings_1 = [c11;c12];
% ------------------------------------------------------------------------------------------------------
%% --- Evaluatiing the resulting offsprings (4) = (3 New produced by the mutation and 1 by the crossover)
parfor  off=1:2
		         %% ---------------------- Evaluate the new solutions---------------------------------------
                  fit_offsprings_1(off)  = RC_Function(offsprings_1(off,:),Dimension,network,neighbourhoud);
                 % -----------------------------------------------------------------------------------------
end
%% ---- Extract The best Behaviour Of the New Produced Offsprings ---
[fit_1,index] = min(fit_offsprings_1);
best1 = [offsprings_1(index,:) min(fit_offsprings_1)];
Offspring_resulted = [Offspring_resulted;best1];
evolving_pm_direction(index) = evolving_pm_direction(index) + 1;
%% --- Apply the same process to the offspring # 2 -------------------------------
%% ---- Apply bit-flip mutation to the first variant of offspring # 2 ------------
for n=1:Dimension
    if rand <= (pm * vv)
        c21(n) = 1 -  c21(n);
    end
end
% ------------------------------------------------------------------
%% ---- Apply bit-flip mutation to the second variant of offspring # 2 ------------
for n=1:Dimension
    if rand <= pm
        c22(n) = 1 -  c22(n);
    end
end
% ------------------------------------------------------------------
%% ---- Creating an intermediate population to stock the produced offsprings ---------------------------
offsprings_2 = [c21;c22];
% ------------------------------------------------------------------------------------------------------
%% --- Evaluatiing the resulting offsprings (4) = (3 New produced by the mutation and 1 by the crossover)
parfor  off=1:2
		         %% ---------------------- Evaluate the new solutions---------------------------------------
                  fit_offsprings_2(off)  = RC_Function(offsprings_2(off,:),Dimension,network,neighbourhoud);
                 % -----------------------------------------------------------------------------------------
end
%% ---- Extract The best Behaviour Of the New Produced Offsprings ---
[fit_2,index] = min(fit_offsprings_2);
best2 = [offsprings_2(index,:) min(fit_offsprings_2)];
evolving_pm_direction(index) = evolving_pm_direction(index) + 1;
Offspring_resulted = [Offspring_resulted;best2];
Offspring_resulted = [Offspring_resulted;evolving_pm_direction];
end
end
end

