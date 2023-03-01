population_offsrings =  rand(1,10)
zeros = find(population_offsrings < 0.5)
ones  = find(population_offsrings >= 0.5)
population_offsrings(ones) = 1;
population_offsrings(zeros) = 0;
population_offsrings