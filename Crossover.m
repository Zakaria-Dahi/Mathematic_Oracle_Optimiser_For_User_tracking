function [Offspring_resulted] = Crossover(population,Dimension,k,pcr,cr)
c1 = [];
c2 = [];
if rand <= pcr
    cross_point = sort(randi([1,Dimension],1,2));
    cross_point_one = cross_point(1);
    cross_point_two = cross_point(2);  
    c1 = [population(cr(k,1),1:cross_point_one) population(cr(k,2),(cross_point_one+1):cross_point_two) population(cr(k,1),(cross_point_two + 1):Dimension)];
    c2 = [population(cr(k,2),1:cross_point_one) population(cr(k,1),(cross_point_one+1):cross_point_two) population(cr(k,2),(cross_point_two + 1):Dimension)];
else
    c1 = population(cr(k,1),:);
    c2 = population(cr(k,2),:);
end
Offspring_resulted = [c1;c2];
end

