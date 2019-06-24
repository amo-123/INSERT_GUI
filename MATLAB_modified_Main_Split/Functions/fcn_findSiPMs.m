function [active_sipms] = fcn_findSiPMs( x_scint, y_scint, Par, n_sipms)
col_vect = {'r', 'y', 'g', 'b'};
value_max = Par.cryst_lung_x;
%% Funzione per selezionare gli n_sipms SiPMs collocati nell'intorno della posizione di scintillazione

x_scint = x_scint*ones(1,Par.num_detect_x*Par.num_detect_y);
y_scint = y_scint*ones(1,Par.num_detect_x*Par.num_detect_y);

x_dist = x_scint - Par.X_c;
y_dist = y_scint - Par.Y_c;
dist = sqrt(x_dist.^2 + y_dist.^2);
sorted_dist = sort(dist, 'ascend');
active_sipms = cell(numel(n_sipms),1);

for i = 1: numel(n_sipms)
    threshold = sorted_dist(n_sipms(i));
    active_sipms{i,1} = find(dist <= threshold); 
    sorted_dist = sorted_dist(n_sipms(i)+1:end);
    dist(active_sipms{i,1}) = value_max;
    
end

PlotGeometry(Par);
hold on
sipm_dim = Par.detect_dim*ones(1,Par.num_detect_x*Par.num_detect_y);
sipm_corners_x = Par.X_c - sipm_dim/2;
sipm_corners_y = Par.Y_c - sipm_dim/2;

for j = numel(n_sipms):-1:1
    for i = 1:numel(active_sipms{j})
        rectangle('Position', [sipm_corners_x(active_sipms{j,1}(i)), sipm_corners_y(active_sipms{j,1}(i)), Par.detect_dim, Par.detect_dim]', 'FaceColor', char(col_vect(j)))
    end 
end


plot(x_scint, y_scint, '*', 'MarkerSize',15, 'MarkerEdgeColor','b')

end