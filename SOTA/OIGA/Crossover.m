function [Offspring_resulted] = Crossover(population,Dimension,k,pcr,cr)
%% ---- Apply the two-point crossover on the formed couples ---------------
c1 = [];
c2 = [];
if rand <= pcr
    cross_point_one = randi(Dimension);
    cross_point_two = randi(Dimension);  
    % -- This condition is to avoid that the first point is bigger than the indice of the second one because it   will generate a wrong behaviour or a bug ---
    if cross_point_one > cross_point_two
        a = cross_point_one;
        cross_point_one = cross_point_two;
        cross_point_two = a;
    end
    c1 = [population(cr(k,1),1:cross_point_one) population(cr(k,2),(cross_point_one+1):cross_point_two) population(cr(k,1),(cross_point_two + 1):Dimension)];
    c2 = [population(cr(k,2),1:cross_point_one) population(cr(k,1),(cross_point_one+1):cross_point_two) population(cr(k,2),(cross_point_two + 1):Dimension)];
else
    c1 = population(cr(k,1),:);
    c2 = population(cr(k,2),:);
end
Offspring_resulted = [c1;c2];
end

