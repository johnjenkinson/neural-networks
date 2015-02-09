% call: MPnetwork_JAJ_AM.m
% John Jenkinson, UTSA ECE 2014
% Azima Mottaghi, ITSA 2014
%
% Modeling a McCulloch-Pitts network
% with 100 neurons arranged linearly
% in space with connections only
% running from right to left.
% number of neurons: nneurons 
% number of states: states 
% activity of the first neuron:
% initial_condition
%
function(active,noisyactive,weight,noisyweight)...
=mpnetwork(nneurons,states,initial_condition)

weight=zeros(states);
for k=2:states
    weight(k,k-1)=1; % weight matrix j->i
end
noisyweight=weight;
active=zeros(states); % activation matrix
active(1,1)=initial_condition;
Vth=0; % offset voltage
for k=2:states
    for n=1:s-1
        if(weight(k,:)*active(:,n)-Vth>0)
            active(k,n+1)=1;
        end
    end
end

% modeling weight matrix corrupted 
% by noise
noisyactive=zeros(states); % noisy activity matrix
noisyactive(1,1)=initial_condition;
for b=2:states
    if(rand(1)>0.85)
        noisyactive(b,1)=1;
    end
    for y=1:s-1
        if(noisyweight(b,:)*noisyactive(:,y)-Vth>0)
            noisyactive(b,y+1)=1;
        end
    end
end
