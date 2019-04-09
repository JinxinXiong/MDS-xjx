function [Dist] = DistMatrix(pt)
    num_pt = size(pt,1);
    warning off
    % output off
    Dist=double(zeros(num_pt,num_pt));
    for index_1=1:num_pt
        Dist(index_1,:) = sqrt(sum((bsxfun(@minus,pt(index_1,:),pt)).^2,2));
    end;
end