%% Load Data
%clear
%close all
%load 1KSphereUniformSampled
load cowPointCloud
%% Generate the pair-wise distance matrix
num_pt = size(pt,1);
warning off
% output off
Dist=double(zeros(num_pt,num_pt));
rate=1;
for index_1=1:num_pt
    Dist_Truth(index_1,:) = sqrt(sum((bsxfun(@minus,pt(index_1,:),pt)).^2,2));
end;
%Dist_Truth is a symmetric matrix with pair-wise distance
%% Generate weight for random-missing distance (Weight==1 means available distance)
rate = 0.2; %Assume only 10% random sampled Distance are available
Weight=rand(num_pt,num_pt);
Weight(Weight>1-rate)=1;
Weight(Weight<1)=0;
Weight(Weight>0)=1;
for i=1:num_pt
    Weight(i,i)=1;
    for j=i+1:num_pt
        Weight(i,j)=Weight(j,i);
    end;
end;
Dist = Dist_Truth.*Weight;
% The distance we know is denoted as Dist, The diagonal component must be
% 0, other value equals to 0 means unknown distances