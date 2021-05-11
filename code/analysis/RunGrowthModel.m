function [B,b] = RunGrowthModel(Input,D,PD,Params,density_val,InitialSEED)

% This function will generate a network while accounting for changes in the
% distances of nodes over time.
%
% Written by Stuart Oldham, 2018

%% Check inputs
if iscell(D)
    N = length(D{1});
    Dists = D;
else 
    N = length(D);
    Dists{1} = D;
end

TopoTypes = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

TopoType = TopoTypes{Input.ModelNum};

if iscell(D)
    if length(D) ~= length(density_val)  
        error('D and density_val must be of the same length')
    end
end

Exponents = [Params(1) Params(2) Params(5)];
Alphas = [1 Params(3) Params(4)];

modelvar{1} = Input.DistFunc;
modelvar{2} = Input.TopoFunc;

if isfield(Input,'PDFunc')
    modelvar{3} = Input.PDFunc;
end

normsum = Input.normsum;

if isempty(InitialSEED)
    InitialSEED = zeros(N);
end

if ~iscell(PD) && ~isempty(PD)
    PD = num2cell(PD,[1 2]);
end

%% Run generative model

% Create network

B = InitialSEED;

for I = 1:length(Dists)
        
    % Find the desired number of edges to add. If it is a value less than 1
    % (i.e., signifying a proportion of edges to add), this is converted to
    % an exact number of edges.
    
        if density_val(I) <= 1
            desiredEdges = ceil(density_val(I)*N*(N-1)/2);
        else
            desiredEdges = density_val(I);
        end 
        
            % Combine the Dists and PDM matrices into a cell
            
                PD2 = [Dists{I} PD];
                
                % If PD is not a cell, make it one
                
                if ~iscell(PD2) && ~isempty(PD2)
                    PD2 = num2cell(PD2,[1 2]);
                end
                          
                % Use the desired model formulation
                
                switch Input.AddMult
                    case 'Add'
            
                        if normsum == 0
                            [B,btemp] = GenMdlAdditive_normmax(B,PD2,desiredEdges,TopoType,modelvar,Exponents,Alphas,Input.epsilon);
                        else
                            [B,btemp] = GenMdlAdditive_normsum(B,PD2,desiredEdges,TopoType,modelvar,Exponents,Alphas,Input.epsilon);          
                        end

                    case 'Mult'
                [B,btemp] = GenMdlMultiplicative(B,PD2,desiredEdges,TopoType,{{modelvar{1},modelvar{3}},modelvar{2}},Exponents,epsilon);             
         
                end
                
        if I == 1 
			b = btemp;
		else
        		if density_val(I-1) <= 1
            			PreviousDesiredEdges = ceil(density_val(I-1)*N*(N-1)/2);
        		else
            			PreviousDesiredEdges = density_val(I-1);
        		end 

			b = [b; btemp(PreviousDesiredEdges+1:desiredEdges)];
		end
        
end

end