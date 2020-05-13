function [ind] = arrlocarr(A,B)
%arrlocarr finds index of common elements array A and B. Index is for B
    v = intersect(A,B);
    ind = find(ismember(B,v));
end

