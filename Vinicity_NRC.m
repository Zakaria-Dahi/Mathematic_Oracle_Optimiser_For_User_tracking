function [Final_Result] = Vinicity_NRC(Actual_NRC,Actual_RC,Actual_neighbourhoud,History_NRC_ID)
Actual_vinicity = 0;
Actual_NRC_ID   = [];
Non_reporting_cells_Neighbourhoud = [];
for i=2:7 
    check = 0;
    for w=1:length(History_NRC_ID)
        if Actual_neighbourhoud(Actual_NRC,i) == History_NRC_ID(w)
            check = 1;
        end
    end
    if check == 0 
        if (Actual_neighbourhoud(Actual_NRC,i) ~= 1000) && (Actual_neighbourhoud(Actual_NRC,i) ~= Actual_NRC) 
            Response =  isempty(find((Actual_neighbourhoud(Actual_NRC,i) ==  Actual_RC) == 1));
            if Response == 1
                        Non_reporting_cells_Neighbourhoud = [Non_reporting_cells_Neighbourhoud  Actual_neighbourhoud(Actual_NRC,i)];
                        Actual_NRC_ID = [Actual_NRC_ID Actual_neighbourhoud(Actual_NRC,i)];
                        Actual_vinicity = Actual_vinicity + 1;
                        History_NRC_ID = [History_NRC_ID Actual_neighbourhoud(Actual_NRC,i)];
            end
        end
    end
end
if isempty(Non_reporting_cells_Neighbourhoud) == 0
    Size_NRC = length(Non_reporting_cells_Neighbourhoud);
    for j=1:Size_NRC
        result = 0;
        result = Vinicity_NRC(Non_reporting_cells_Neighbourhoud(j),Actual_RC,Actual_neighbourhoud,History_NRC_ID);
        Actual_vinicity = Actual_vinicity +  result(1); 
        History_NRC_ID = result(2:length(result));
    end
end
Final_Result = [Actual_vinicity History_NRC_ID];
end

