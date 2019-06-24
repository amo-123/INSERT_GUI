function [Phys_order,ASIC_order]=INSERT_reorder(modality)

%% Preclinical & Little INSERT
if strcmp(modality,'Preclinical') || strcmp(modality,'Little_INSERT')
    % from ADC to ASIC order
    ASIC_order = [20  2 19  1 22  4 21  3 24  6 23  5 26  8 25  7 28 10 27  9 30 12 29 11 32 14 31 13 34 16 33 15 36 18 35 17];
    if strcmp(modality,'Preclinical')
        % from ASIC order to Physical order;
        Phys_order = [35 30 36 34 29 28 23 22 24 18 16 17 10 11 4 6 12 5 33 25 31 32 26 27 20 21 19 13 15 14 9 8 3 1 7 2];
    elseif strcmp(modality,'Little_INSERT')
        % from ASIC order to Physical order;
        Phys_order = [13 1 2 3 6 5 4 16 22 19 20 24 34 29 27 26 25 30 31 32 33 36 35 28 12 7 8 9 11 10 18 14 15 17 23 21];
    end
end

%% Clinical
if strcmp(modality,'Clinical')
    % from ADC to ASIC order
    ASIC_order = [55 37 19 01 56 38 20 02 57 39 21 03 58 40 22 04 59 41 23 05 60 42 24 06 61 43 25 07 62 44 26 08 63 45 27 09 64 46 28 10 65 47 29 11 66 48 30 12 67 49 31 13 68 50 32 14 69 51 33 15 70 52 34 16 71 53 35 17 72 54 36 18];
    % from ASIC order to Physical order;    
    Preclinical_order = [35 30 36 34 29 28 23 22 24 18 16 17 10 11 4 6 12 5 33 25 31 32 26 27 20 21 19 13 15 14 9 8 3 1 7 2];
    %Preclinical_order = [34 36 33 15 18 16 35 32 31 13 14 17 28 30 29 11 12 10 27 25 26 8 7 9 20 23 24 6 5 2 21 22 19 4 1 3];
    
    Phys_order = [Preclinical_order Preclinical_order+36];
end

end










 

    
 