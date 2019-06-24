function PlotGeometry(Par, varargin)

sipm_dim = Par.detect_dim*ones(1,Par.num_detect_x*Par.num_detect_y);
sipm_corners_x = Par.X_c - sipm_dim/2;
sipm_corners_y = Par.Y_c - sipm_dim/2;

figure, hold on
axis([0 Par.cryst_lung_x 0 Par.cryst_lung_y])
set(gca,'Ydir','reverse')
for i = 1 :Par.num_detect_x*Par.num_detect_y
    rectangle('Position', [sipm_corners_x(i), sipm_corners_y(i), Par.detect_dim, Par.detect_dim]')
end 

plot_spots = isempty(varargin);

if (plot_spots == 0) 
    spots = varargin{1};
    for i = 1:numel(spots)
        plot(spots{i}.hole(1), spots{i}.hole(2), '*')
        text(spots{i}.hole(1), spots{i}.hole(2) + 0.5, num2str(i))
    end
end
 return