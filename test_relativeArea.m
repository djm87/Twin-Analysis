    for j=1:length(pairs)
        n1=pairs(j,1)
        n2=pairs(j,2)
        n1F=FamilyID(n1==nodeID)
        n2F=FamilyID(n2==nodeID)
        RFA1(j)=(FA(n1F)-FA(n2F))/(FA(n1F)+FA(n2F))
        RFA2(j)=(FA(n2F)-FA(n1F))/(FA(n1F)+FA(n2F))
    end