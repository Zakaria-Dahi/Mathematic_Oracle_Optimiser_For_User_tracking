function [index] = Elitist_selection(temp_fit,temp_pop,indiv)
[n,m] = size(temp_pop);
sorted_temp_fit = sort(temp_fit,'ascend');
uniqueX = unique(sorted_temp_fit);
countOfX = hist(sorted_temp_fit,uniqueX);
indexToRepeatedValue = (countOfX~=1);
longeur = length(find(indexToRepeatedValue == 1));
Member = find(countOfX > 1);
Doubles = uniqueX(Member);
Counting = zeros(1,length(Doubles));
index = [];
individual = 0;
f = 1;
while individual < indiv
   check = 0;
   for k=1:length(Doubles)
       if Doubles(1,k) == sorted_temp_fit(f)
          Counting(1,k) =   Counting(1,k) + 1; 
          numero = k;
          check = 1;
       end
   end
   if check == 0
      for m=1:n
         if sorted_temp_fit(f) == temp_fit(m)
            index = [index m];
            individual =  individual + 1;
         end
      end
   else
      if Counting(1,numero) == 1
        for m=1:n
          if sorted_temp_fit(f) == temp_fit(m)
            if Counting(1,numero) == 1
               index = [index m];           
               individual =  individual + 1;
               Counting(1,numero) =   Counting(1,numero) + 1; 
            end
          end
        end
      end
   end
   f = f + 1;
end
end

