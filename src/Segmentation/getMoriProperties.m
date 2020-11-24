function [moriOut,nMori,moriUnknown,moriLAGB] = getMoriProperties(mori)  

    %Main loop
    cnt=1;
    for i=1:length(mori)
        
        %If the misorientation relationship is a twin define the type
        if mori{i}.istwin
            %Define the transformation types
            tType{1}=rotation.byMatrix([-1  0  0;0 -1  0;0  0 1],mori{i}.CS);
            tType{2}=rotation.byMatrix([1  0  0;0 -1  0;0  0 -1],mori{i}.CS);   
        end
        %Loop over the mapping specified (j>1 is double mis, tripple mis,
        %etc...)
        for j=1:length(mori{i}.u1)
            %Define the Rtw such that Rtw transforms crystal to twin or
            %phase to phase if not twin
            mori{i}.Rtw(j)=orientation.map(mori{i}.u1,mori{i}.v1,...
                mori{i}.u2,mori{i}.v2); 

            if mori{i}.istwin
                RMTtmp=mori{i}.Rtw(j)'*tType{mori{i}.actType(j)}*mori{i}.Rtw(j);
            else
                RMTtmp=mori{i}.Rtw(j);
            end
            
            % Create the full misorientation
            if j==1
                mori{i}.RMT=RMTtmp;
            else
                mori{i}.RMT=mori{i}.RMT*mori{i}.CS*RMTtmp;
            end

            %Specify a generation multiplier for identifying
            %generation during family tree calculations
            mori{i}.genMultiplier=j;
        end

        % Compute misorientation axis and angle
        mori{i}.axis=round(mori{i}.RMT.axis);
        mori{i}.angle=angle(mori{i}.RMT)/degree;            

        %For each misorientation compute variants and relevant quantities
        for k=1:length(mori{i}.variantsToUse)
            uVar=mori{i}.variantsToUse(k);
            moriOut{cnt}=mori{i};
            moriOut{cnt}.RMT=mori{i}.RMT(uVar);
            moriOut{cnt}.angle=mori{i}.angle(uVar);
            moriOut{cnt}.axis=mori{i}.axis(uVar);
            moriOut{cnt}.sS=slipSystem(mori{i}.u2(end),mori{i}.u1(end)).symmetrise('antipodal');
            moriOut{cnt}.axisVariants=moriOut{cnt}.axis.symmetrise('unique');
            moriOut{cnt}.Rtw=orientation.map(moriOut{cnt}.sS.n,mori{i}.CS.cAxis,moriOut{cnt}.sS.b,mori{i}.CS.aAxis); 
            moriOut{cnt}.name=mori{i}.name{k};
            moriOut{cnt}.variantsToUse=1;
            cnt=cnt+1;
        end
    end
    moriOut{cnt}.name='unknown internal';

    %store nTwin and twinUnknown index
    nMori=cnt-1;
    moriUnknown=cnt;
    moriLAGB=cnt+1;
    
end

