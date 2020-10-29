function pCell=createParamCell(PCscale,Wscale,thetascale,mList,nList,eList,consList,paramsDefault,symmetric,conUtil,conParam,learningRates)
% createParamCell - This function simply creates a list of parameters to run
%%
% @param: Lists - the lists of configurations
% @param: paramsDefault - starting struct of parameters that are constant
% @return: pCell - Cell of params structures over which to loop
%%
%% Create a cell array of type vectors
gammaVec={};
iC=1;
ceoActStartT=0;
for watchP = Wscale
    for scale = PCscale
        normC=scale.*(1-watchP);
        normS=(1-scale).*(1-watchP);
        if normC<0 || normC>1 || normS <0 ||normS>1 || watchP<0 || watchP>1
            error("Type probabilities out of bounds. Check config!")
        end
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
                    for ceoStartT=ceoActStartT
                        for learningRate=learningRates
                            for gamma=gammaVec
                                params=paramsDefault;
                                params.learningRate=learningRate;
                                params.n=n;
                                params.m=m;
                                params.e=e;
                                params.gamma=gamma{:};
                                params.thetaD=theta{:};
                                params.cons=cons;
                                params.ceoStartT=ceoStartT;
                                params.conParam=conParam;
                                % Check which concavity to use for social benefit
                                if conUtil == 1 % Concave Benefit
                                    params.conUtil=1;
                                    pCell{i}=params;
                                    i=i+1;
                                elseif conUtil == 0 % Linear Benefit
                                    params.conUtil=0;
                                    pCell{i}=params;
                                    i=i+1;
                                else % Both and Convex benefit
                                    params.conUtil=1;
                                    pCell{i}=params;
                                    i=i+1;
                                    params.conUtil=0;
                                    params.conParam=0;
                                    pCell{i}=params;
                                    i=i+1;
                                    
                                    params.conUtil=1;
                                    params.conParam=1./conParam;
                                    pCell{i}=params;
                                    i=i+1;
                                    
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
disp(['Finished creating parameter cell!'])
end
