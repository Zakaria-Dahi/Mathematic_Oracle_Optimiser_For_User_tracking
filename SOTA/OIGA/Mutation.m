function [Offspring_resulted] = Mutation(c1,c2,Dimension,pm)
%% ---- Apply bit-flip to the first offsprings ---------------------
for n=1:Dimension
    if rand <= pm
        c1(n) = 1 -  c1(n);
    end
end
% ------------------------------------------------------------------
%% ---- Apply bit-flip mutation to the second offspring ------------
for n=1:Dimension
    if rand <= pm
        c2(n) = 1 -  c2(n);
    end
end
% ------------------------------------------------------------------
Offspring_resulted = [c1;c2];
end

