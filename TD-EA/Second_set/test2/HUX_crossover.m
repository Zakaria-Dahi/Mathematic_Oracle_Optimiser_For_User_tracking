function [offsprings] = HUX_crossover(population,indiv,distance,cr)
%% ----- Recombinaison des individus HUX ---------------------------
offsprings = [];
offspring_one = [];
offspring_two = [];
if mod(indiv,2) == 0
   for i=1:2:indiv
   if rand <= cr
    soustract = population(i,:) - population((i+1),:);
    differ = find (soustract ~= 0);
    hamming_dist = length(differ);
    threshold = (fix(hamming_dist/2));
    if threshold > distance
       bit_flip = randi([1,hamming_dist],[1,threshold]);
       uniqueX = unique(bit_flip);
       countOfX = hist(bit_flip,uniqueX);
       indexToRepeatedValue = (countOfX~=1);
       longeur = length(find(indexToRepeatedValue == 1));       
       while longeur ~= 0
           bit_flip = randi([1,hamming_dist],[1,threshold]);
           uniqueX = unique(bit_flip);
           countOfX = hist(bit_flip,uniqueX);
           indexToRepeatedValue = (countOfX~=1);
           longeur = length(find(indexToRepeatedValue == 1));
       end
       offspring_one = population(i,:);
       offspring_two = population(i+1,:);
       offspring_one(differ(bit_flip)) = abs(offspring_one(differ(bit_flip)) - 1);
       offspring_two(differ(bit_flip)) = abs(offspring_two(differ(bit_flip)) - 1);
       offsprings = [offsprings;offspring_one;offspring_two];
    end
   end
   end
else
   for i=1:2:(indiv - 1)
   if rand <= cr
    soustract = population(i,:) - population((i+1),:);
    differ = find (soustract ~= 0);
    hamming_dist = length(differ);
    threshold = (fix(hamming_dist/2));
    if threshold > distance
       bit_flip = randi([1,hamming_dist],[1,threshold]);
       uniqueX = unique(bit_flip);
       countOfX = hist(bit_flip,uniqueX);
       indexToRepeatedValue = (countOfX~=1);
       longeur = length(find(indexToRepeatedValue == 1));       
       while longeur ~= 0
           bit_flip = randi([1,hamming_dist],[1,threshold]);
           uniqueX = unique(bit_flip);
           countOfX = hist(bit_flip,uniqueX);
           indexToRepeatedValue = (countOfX~=1);
           longeur = length(find(indexToRepeatedValue == 1));
       end
       offspring_one = population(i,:);
       offspring_two = population(i+1,:);
       offspring_one(differ(bit_flip)) = abs(offspring_one(differ(bit_flip)) - 1);
       offspring_two(differ(bit_flip)) = abs(offspring_two(differ(bit_flip)) - 1);
       offsprings = [offsprings;offspring_one;offspring_two];
    end
   end
   end
   if rand <= cr
    soustract = population(indiv,:) - population(1,:);
    differ = find (soustract ~= 0);
    hamming_dist = length(differ);
    threshold = (fix(hamming_dist/2));
    if threshold > distance
       bit_flip = randi([1,hamming_dist],[1,threshold]);
       uniqueX = unique(bit_flip);
       countOfX = hist(bit_flip,uniqueX);
       indexToRepeatedValue = (countOfX~=1);
       longeur = length(find(indexToRepeatedValue == 1));       
       while longeur ~= 0
           bit_flip = randi([1,hamming_dist],[1,threshold]);
           uniqueX = unique(bit_flip);
           countOfX = hist(bit_flip,uniqueX);
           indexToRepeatedValue = (countOfX~=1);
           longeur = length(find(indexToRepeatedValue == 1));
       end
       offspring_one = population(indiv,:);
       offspring_one(differ(bit_flip)) = abs(offspring_one(differ(bit_flip)) - 1);
       offsprings = [offsprings;offspring_one];
    end
   end
end
end

