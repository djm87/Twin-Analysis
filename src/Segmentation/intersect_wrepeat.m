function [ind] = intersect_wrepeat(A,B)
%arrlocarr finds index of common elements array A and B
    v = intersect(A,B);
    ind = find(ismember(B,v));
end

