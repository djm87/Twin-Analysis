function [ebsd,G_Complete] = fragmentReconstruction(G_Complete,ebsd,grains)
%fragmentReconstruction combines fragments that appear are one twin but are
%broken up due to resolution issue. Helps Fix count stats.. 
    angleDiffTol=5*degree;
    [omega,a_mag,b_min] = fitEllipse(grains);

    G_Complete.Nodes.MergeTwin = zeros(length(G_Complete.Nodes.Id),1);
    cnt=0;
    if true
        %loop through each group/grain
        for i=1:max(G_Complete.Edges.Group) 
            egroupId = find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
            ngroupId = find((i==G_Complete.Nodes.Group)==true);
            nFamily = G_Complete.Nodes.FamilyID(ngroupId);
            nType = G_Complete.Nodes.Type(ngroupId); 
            nId = G_Complete.Nodes.Id(ngroupId);
            eType = G_Complete.Edges.type(egroupId);
            eVote = G_Complete.Edges.Vote(egroupId,:);
            ePairs = G_Complete.Edges.pairs(egroupId,:);
            eFamily = G_Complete.Edges.FamilyID(egroupId,:);
            eGlobalId = G_Complete.Edges.GlobalID(egroupId);

            MergeTwin=zeros(length(nId),1);
            for j=1:length(nId)
                actFam=nFamily(j);
                if nType(j)>0 %If the fragment is not the parent
                    for k=1:length(nId)
                       if k~=j && nFamily(k)==actFam %Only the same family if the same orientation
                          %Test if grains should be merged 
                          centroidj=grains(nId(j)).centroid;
                          centroidk=grains(nId(k)).centroid;
                          radiusSum=a_mag(nId(k))+a_mag(nId(j));
                          angleDiffEllipse=abs(omega(nId(k))-omega(nId(j)));

                          %Construct a line and find min distance from other
                          %centroid
                          m = tan( omega(nId(j)));
                          b = centroidj(2) - (m * centroidj(1) );
                          B=-m;
                          C=-b;
                          A=1;

                          distFromLine=abs(A*centroidk(2)+B*centroidk(1)+C)/sqrt(A^2+B^2);

                          %Convert line to vector 
                          centroidDiff=norm(centroidk-centroidj);


                          if  angleDiffEllipse<angleDiffTol && distFromLine<1.25*b_min(nId(j)) && centroidDiff<5*radiusSum
                                if MergeTwin(j)~=0 && MergeTwin(k)==0
                                   MergeTwin(k)= MergeTwin(j);
                                elseif MergeTwin(k)~=0 && MergeTwin(j)==0
                                   MergeTwin(j)=MergeTwin(k);
                                elseif MergeTwin(k)==0 && MergeTwin(j)==0
                                   cnt=cnt+1;
                                   MergeTwin(j)=cnt;
                                   MergeTwin(k)=cnt;
                                end
                          end
                       end %
                    end %loop k over fragments

                end %Type

            end
            G_Complete.Nodes.MergeTwin(nId)=MergeTwin;

        end
    end
end

