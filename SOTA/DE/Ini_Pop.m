function [population] = Ini_Pop(indiv,Dimension,Prob)

%% -------- Generate the intial population According to a Standard Uniform Disribution U(0,1) ------------------
population = rand(indiv,Dimension);
% --------- Generate the initial Binary Population (0,1) According to the Probability "Prob" -------------------
population = (population <= Prob) ;
% --------- Eliminate the non feasible solutions ---------------------------------------------------------------
for w=1:indiv
    if population(w,:) == zeros(1,Dimension)
        % ------ Assign a reporting cell at random -----
        population(w,randi(Dimension)) = 1;
    end
end
% --------------------------------------------------------------------------------------------------------------
end

