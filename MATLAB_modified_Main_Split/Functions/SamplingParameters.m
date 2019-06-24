function [ Par ] = SamplingParameters( Par )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Par.x_knot=Par.sampling/2:Par.sampling:Par.cryst_lung_x;%definition of coordinates of a finite number of nodes where the LRF has to be computed by the estimator
Par.y_knot=Par.sampling/2:Par.sampling:Par.cryst_lung_y;
% Par.LRF_x_knot=Par.pixel_min/2:Par.pixel_min:Par.cryst_lung_x;
% Par.LRF_y_knot=Par.pixel_min/2:Par.pixel_min:Par.cryst_lung_y;
Par.x_coord=reshape(repmat(Par.x_knot,length(Par.y_knot),1),1,length(Par.x_knot)*length(Par.y_knot));
Par.y_coord=repmat(Par.y_knot,1,length(Par.x_knot));
end

