function [Offspring_resulted] = Mutation(c1,c2,Dimension,pm,evolving_pm_direction,vv,segma,network,neighbourhoud,evaluated_offsprings,individual)
if evaluated_offsprings < individual
if (individual - evaluated_offsprings) == 1
Offspring_resulted = [];
fit_offsprings_1 = zeros(1,2);
offsprings_1 = [];
c11 = c1;
c12 = c1;
for n=1:Dimension
    if rand <= (pm * vv)
        c11(n) = 1 -  c11(n);
    end
end
for n=1:Dimension
    if rand <= pm
        c12(n) = 1 -  c12(n);
    end
end
offsprings_1 = [c11;c12];
parfor  off=1:2
                  fit_offsprings_1(off)  = RC_Function(offsprings_1(off,:),Dimension,network,neighbourhoud);
end
[fit_1,index] = min(fit_offsprings_1);
best1 = [offsprings_1(index,:) min(fit_offsprings_1)];
Offspring_resulted = [Offspring_resulted;best1];
evolving_pm_direction(index) = evolving_pm_direction(index) + 1;
Offspring_resulted = [Offspring_resulted;zeros(1,(Dimension + 1));evolving_pm_direction];
end    

    
    
    
    
    
if (individual - evaluated_offsprings) >= 2
Offspring_resulted = [];
fit_offsprings_1 = zeros(1,2);
fit_offsprings_2 = zeros(1,2);
offsprings_1 = [];
offsprings_2 = [];
c11 = c1;
c12 = c1;

c21 = c2;
c22 = c2;
for n=1:Dimension
    if rand <= (pm*vv)
        c11(n) = 1 -  c11(n);
    end
end
for n=1:Dimension
    if rand <= pm
        c12(n) = 1 -  c12(n);
    end
end
offsprings_1 = [c11;c12];
parfor  off=1:2
                  fit_offsprings_1(off)  = RC_Function(offsprings_1(off,:),Dimension,network,neighbourhoud);
end
[fit_1,index] = min(fit_offsprings_1);
best1 = [offsprings_1(index,:) min(fit_offsprings_1)];
Offspring_resulted = [Offspring_resulted;best1];
evolving_pm_direction(index) = evolving_pm_direction(index) + 1;
for n=1:Dimension
    if rand <= (pm * vv)
        c21(n) = 1 -  c21(n);
    end
end
for n=1:Dimension
    if rand <= pm
        c22(n) = 1 -  c22(n);
    end
end
offsprings_2 = [c21;c22];
parfor  off=1:2
                  fit_offsprings_2(off)  = RC_Function(offsprings_2(off,:),Dimension,network,neighbourhoud);
end
[fit_2,index] = min(fit_offsprings_2);
best2 = [offsprings_2(index,:) min(fit_offsprings_2)];
evolving_pm_direction(index) = evolving_pm_direction(index) + 1;
Offspring_resulted = [Offspring_resulted;best2];
Offspring_resulted = [Offspring_resulted;evolving_pm_direction];
end
end
end

