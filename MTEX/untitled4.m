
v = vector3d.rand(100);
v.x = 0;
v = rotation.rand * v + vector3d.rand;

vsave = v;

%%

s = v(1);

% shift to origin
v = v - s;

% compute normal vector
[n,~] = eig3(v*v);


r = rotation.map(n(1),zvector);

v = r * v;

plot3(v.x,v.y,v.z,'.')

%%

P = [v.x(:),v.y(:)];
T = delaunayn(P);
n = size(T,1);
W = zeros(n,1);
C=0;
for m = 1:n
    sp = P(T(m,:),:);
    [null,W(m)]=convhulln(sp);
    C = C + W(m) * mean(sp);
end

C = vector3d(C(1),C(2),0)./sum(W);



%%

hold on
plot3(C.x,C.y,C.z,'MarkerSize',10)
hold off






