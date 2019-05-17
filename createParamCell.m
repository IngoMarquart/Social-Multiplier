%%
% % createParamCell
% This function simply creates a list of parameters to run
% @param: Lists - the lists of configurations
% @param: paramsDefault - starting struct of parameters that are constant
% @return: pCell - Cell of params structures over which to loop
%%
function pCell=createParamCell(PCscale,Wscale,thetascale,mList,nList,eList,consList,paramsDefault,symmetric)

%% Create a cell array of type vectors
gammaVec={};
iC=1;
for watchP = Wscale
    for scale = PCscale
        normC=scale.*(1-watchP);
        normS=(1-scale).*(1-watchP);
        gammaVec{iC}=[normC,watchP,normS];
        if (normC == normS) || symmetric~=1
            iC=iC+1;
        else
            gammaVec{iC+1}=[normS,watchP,normC];
            iC=iC+2;
        end
        
    end
end
clear iC normC normS scale watchP

%% Create a cell array of theta vectors
thetaVec={};
iz=1;
for scale = thetascale
    if scale==2
        thetaVec{iz}=[2,scale,1,1];
        iz=iz+1;
    elseif symmetric==1
        thetaVec{iz}=[2,scale,1,1];
        thetaVec{iz+1}=[scale,2,1,1];
        iz=iz+2;
    elseif symmetric==-1
        thetaVec{iz}=[scale,2,1,1];
        iz=iz+1;
    else
        thetaVec{iz}=[2,scale,1,1];
        iz=iz+1;
    end
end
clear scale iz

pCell={};
i=1;
for m=mList
    for n=nList
        for e=eList
            for cons=consList
                for theta=thetaVec
                    for gamma=gammaVec
                        params=paramsDefault;
                        params.n=n;
                        params.m=m;
                        params.e=e;
                        params.gamma=gamma{:};
                        params.thetaD=theta{:};
                        params.cons=cons;
                        pCell{i}=params;
                        i=i+1;
                    end
                end
            end
        end
    end
end

end
