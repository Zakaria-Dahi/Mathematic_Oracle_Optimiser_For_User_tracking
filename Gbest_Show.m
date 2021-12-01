function [Final_Result] = Gbest_Show(individual,Dimension,network,neighbourhoud)
Beta = 10;
Lu = 0;
Reporting_cells = [];
Vinicity_vector = [];
Non_reporting_cells_vector = [];
for cell=1:Dimension
    if individual(cell) == 1
        Lu = Lu + network(cell,1);
        Reporting_cells = [Reporting_cells cell];
    end
end
Size_RC = length(Reporting_cells);
for RC=1:Size_RC 
    Non_reporting_cells = [];
    NRC_ID = [];
    vinicity = 1;
    for i=2:7
      if neighbourhoud(Reporting_cells(RC),i) ~= 1000
          Response =  isempty(find((neighbourhoud(Reporting_cells(RC),i) ==  Reporting_cells) == 1));
          if Response == 1
              Non_reporting_cells = [Non_reporting_cells  neighbourhoud(Reporting_cells(RC),i)];
              NRC_ID = [NRC_ID neighbourhoud(Reporting_cells(RC),i)];
              vinicity = vinicity + 1;
          end
      end
    end
    if isempty(Non_reporting_cells) == 0
     Size_NRC = length(Non_reporting_cells);
        for j=1:Size_NRC
             result = 0;
             result = Vinicity_NRC(Non_reporting_cells(j),Reporting_cells,neighbourhoud,NRC_ID);
             vinicity = vinicity +  result(1); 
             NRC_ID = result(2:length(result));
        end
    end  
Vinicity_vector = [Vinicity_vector vinicity];
Non_reporting_cells_vector = [Non_reporting_cells_vector mat2cell(NRC_ID)];
end
Vinicity_Reporting_Cells = [Reporting_cells ;Vinicity_vector];
Non_reporting_cell_ID = [];
Non_reporting_cell_Vinicity = [];
for cell=1:Dimension
    if individual(cell) == 0
           Non_reporting_cell_ID = [Non_reporting_cell_ID cell];
    end
end
for k=1:length(Non_reporting_cell_ID)
    vinicity_record = [];
    for z=1:length(Non_reporting_cells_vector)
        Neighbour_RC = cell2mat(Non_reporting_cells_vector(z));
        Response =  (isempty(find((Non_reporting_cell_ID(k) == Neighbour_RC) == 1)));
        if Response == 0
           vinicity_record = [vinicity_record Vinicity_Reporting_Cells(2,z)];
        end
    end
    Non_reporting_cell_Vinicity = [Non_reporting_cell_Vinicity max(vinicity_record)]; 
end
Vinicity_Non_Reporting_Cells = [Non_reporting_cell_ID;Non_reporting_cell_Vinicity];
Paging_Vinicity_Data = [Vinicity_Non_Reporting_Cells Vinicity_Reporting_Cells];
Pc = 0;
for id_one=1:Dimension
    for id_two=1:Dimension
        if Paging_Vinicity_Data(1,id_two) == id_one
           Pc = (network(id_one,2) * Paging_Vinicity_Data(2,id_two)) + Pc;
        end
    end
end
result = (Beta * Lu) + Pc;
Final_Result = [result Lu Pc mat2cell(Reporting_cells) mat2cell(Non_reporting_cell_ID)];
end

