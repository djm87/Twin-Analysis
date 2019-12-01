function G = InitializeGraph(ebsd,grains,twin,Mistol,meanMistol,...
    meanMistolRelaxed,doplot)
%This function initializes graph objects for grains. 
%Initialization entails computation of all local grain properties such as
%orientation, centroid, aspectRatio, etc.. and is stored in edge pair
%(intergranular) or nodes (grains). 

%Compute grain neighbors

[counts,pairs] = neighbors(grains);

%Initialize graph
s=pairs(:,1);
t=pairs(:,2);
G=graph(s,t);

% Add node labels and other grain level information 
G.Nodes.Name=cellstr(int2str(grains.id));
G.Nodes.Id=str2num(cell2mat(G.Nodes.Name));
G.Nodes.Area=grains.area; 
G.Nodes.Perimeter=grains.perimeter; 
G.Nodes.AspectRatio=grains.aspectRatio;
G.Nodes.Paris=grains.paris;
G.Nodes.centroids=grains.centroid;
G.Nodes.meanOrientation=grains.meanOrientation;
G.Nodes.Properties.UserData.mineral=grains.mineral; %For single phase material

% %DS COMMENT
% % Compute grain boundary for each grain fragment (i.e. node) 
% G.Nodes.Gb=cell(length(grains),1);
% 
% for i=1:length(grains)
%     G.Nodes.Gb{i}=grains(i).boundary;
% end

% Add intergranular information
G.Edges.pairs=pairs; %this contains the grain ids corresponding to grains.id

%Compute misorientation between all pairs
mori=inv(grains(G.Edges.pairs(:,1)).meanOrientation).*...
    grains(G.Edges.pairs(:,2)).meanOrientation; 

%Test if mean misorientation is a twin type so we can cluster grains
[combine,type] = TestTwinRelationship(mori,meanMistol,twin);

G.Edges.type=type; %Twin relation type (# from twin list definitions)
G.Edges.combine=combine; %whether pairs should be grouped into grains

% %DS COMMENT
% % Make a list of boundaries shared between grains connected by an edge.
% % This is used in relative boundary fraction calculations
% G.Edges.Gb=cell(length(G.Edges.pairs),1);
% G.Edges.ebsdId=cell(length(G.Edges.pairs),1);
% 
% grain1=grains(G.Edges.pairs(:,1));
% grain2=grains(G.Edges.pairs(:,2));
% 
% 
% for i=1:length(grain1)
%     id1=grain1(i).boundary.ebsdId;
%     id2=grain2(i).boundary.ebsdId;
%     [boundaryEbsdId,loc]=intersect(id1,id2,'rows');
%     G.Edges.Gb{i}=grain1.boundary(loc);
%     G.Edges.ebsdId{i}=boundaryEbsdId;
% end
% %DS COMMENT

%Overlayer graph on grains
if doplot==true
    figure;
    h=plot(grains,grains.meanOrientation,'Micronbar','off');
    set(gca,'Units','normalized');
    hold on 
    p=plot(G,'XData',G.Nodes.centroids(:,1),'YData',G.Nodes.centroids(:,2));
    hold off
    p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
end

end

