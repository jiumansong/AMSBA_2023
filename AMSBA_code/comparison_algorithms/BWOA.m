% The Whale Optimization Algorithm
function [Leader_pos,Convergence_curve]=BWOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems
K = SearchAgents_no;

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=[];
m=2500;
t=0;% Loop counter
FEs=0;
% Main loop
while FEs<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            FEs=FEs+1;
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end
        
    end
    
    a=2-FEs*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+FEs*((-1)/Max_iter);
    
    % Update the Position of search agents
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
                
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                
            end
            
        end
        
        
        %          r= cauchyrnd(0, 1, 1);t=randperm(1); cauchy=r(t(1),:);
        %         Moth_pos_m_cauchy=Positions(i,:)*(1+1*cauchy);
        %                 Moth_fitness_m_cauchy=fobj(Moth_pos_m_cauchy);
        %         Moth_fitness_s=fobj(Positions(i,:));
        %
        %         Moth_fitness_comb=[Moth_fitness_m_cauchy,Moth_fitness_s];
        %         [~,m]=min(Moth_fitness_comb);
        %         if m==1
        %             Positions(i,:)=Moth_pos_m_cauchy;
        %
        %         else
        %             Positions(i,:)=Positions(i,:);
        %         end
        
        Moth_pos_m_levy1=Positions(i,:)*(1+Levy1(1));
              
        Flag4ub=Moth_pos_m_levy1>ub;
        Flag4lb=Moth_pos_m_levy1<lb;
        Moth_pos_m_levy1=(Moth_pos_m_levy1.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            FEs=FEs+2;
        Moth_fitness_m_levy1=fobj(Moth_pos_m_levy1);
        Moth_fitness_s=fobj(Positions(i,:));
        
        Moth_fitness_comb=[Moth_fitness_m_levy1,Moth_fitness_s];
        [~,m]=min(Moth_fitness_comb);
        if m==1
            Positions(i,:)=Moth_pos_m_levy1;
        else
            Positions(i,:)=Positions(i,:);
        end
        
        
        
    end
            setCan = 1-power(abs((FEs-1)/FEs),m);

    % setCan = (Max_iter-t+1)/Max_iter;
    x = rand();
    while(~(x~= 0.25 && x~= 0.5 && x~=0.75))
        x = rand();
    end
    
    ch(1) = x;
    
    for ant=1:(K)
        ch(ant+1)=4*ch(ant)*(1-ch(ant));
        CH(ant,:) = lb+ch(ant)*(ub-lb);    %ub大
        V = (1-setCan)*Leader_pos+setCan*CH(ant);
        % ObjValV=feval(objfun,V);         %计算函数值
        Flag4ub=V>ub;
        Flag4lb=V<lb;
        V=(V.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            FEs=FEs+1;
        [FitnessV]=fobj(V);%计算适应度值
        %                 [FitnessV,~,~]=fobj(Positions(i,:)',funcNum);
        
        %    if (ObjValV<GlobalMin)
        if (FitnessV<Leader_score)
            Leader_score = FitnessV;
            Leader_pos = V;
            break;
        end
    end
    
    t=t+1;
    Convergence_curve(t)=Leader_score;
    %     [t Leader_score]
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function o=Levy1(d)
% d=1;
beta=3/2;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);

% Eq. (3.9)
o=step;

end




