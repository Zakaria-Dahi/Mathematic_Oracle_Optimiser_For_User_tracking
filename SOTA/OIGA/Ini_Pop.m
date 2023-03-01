function [population] = Ini_Pop(indiv,Dimension)

%% ------------ Generate the intial population According to a Standard Uniform Disribution U(0,1)---------------------
population = round(rand(indiv,Dimension));
% --------------------------------------------------------------------------------------------------------------
end

