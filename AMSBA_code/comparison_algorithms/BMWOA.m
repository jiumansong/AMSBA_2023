%An enhanced associative learning-based exploratory whale optimizer for global optimization
%  An enhanced WOA algorithm with 
% 1- BHC local search, 
% 2-Associative learning mechanism with immediate memory (ALIM)

%  Developed in MATLAB R2013b
%
%  Author and programmer: Ali Asghar Heidari
%
%         e-Mail: as_heidari@ut.ac.ir

% ------------------------------------------------------
 function [Leader_pos,Convergence_curve]=BMWOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
% disp('WOA5 (EWOA=WOA+BHC+ALIM) is now estimating')
% tic
% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems


%Initialize the X of search agents
X=initialization(SearchAgents_no,dim,ub,lb);

% Convergence_curve=zeros(1,Max_iter);
Convergence_curve=[];
FEs=0;
t=0;% Loop counter

% Main loop
while FEs<Max_iter  
  %  t
    for i=1:size(X,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
          FEs=FEs+1;
        % Calculate objective function for each search agent
        fitness =fobj(X(i,:));
        
        
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=X(i,:);
        end
        
    end
    
    a=2-FEs*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+FEs*((-1)/Max_iter);
    
    % Update the Position of search agents
    for i=1:size(X,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
                  b=1;               %  parameters in Eq. (2.5)
%         b=0.3063489;             %golden spiral
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(X,2)
            
            if p<0.5
                
                if abs(A)>=1
                    p4=rand();
                    if p4<0.5
                        
                        rand_leader_index = floor(SearchAgents_no*rand()+1);
                        X_rand = X(rand_leader_index, :);
                        D_X_rand=abs(C*X_rand(j)-X(i,j)); % Eq. (2.7)
                        X(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    elseif p4>=0.5
                        %% -----Associative learning mechanism with immediate memory (ALIM)
                        
                        gx=Leader_pos;
                        g=FEs;
                        j2=SearchAgents_no;
                        j3=j2;
                        G=Max_iter;
                        rand_leader_index = floor(SearchAgents_no*rand()+1);
                        X_rand = X(rand_leader_index, :);
                        X(j3,:)=X(j3,:)+0.001*unifrnd(lb-X(j3,:),ub-X(j3,:))+(2*g/G)*rand(size(gx)).*(gx(1,:)-X(j3,:))+(1-g/G)*rand(size(gx)).*(X_rand(1,:)-X(j3,:));%step away like PSO
                        %                         X(j3,:)=2*rand()*X_rand(j)+0.001*unifrnd(lb-X(j3,:),ub-X(j3,:))+(2*g/G)*rand(size(gx)).*(gx(1,:)-X(j3,:))+(1-g/G)*rand(size(gx)).*(X_rand(1,:)-X(j3,:));%step away like PSO
                        
                    end
                    
                elseif abs(A)<1
                    p2=rand();
                    if p2<0.5
                        
                        D_Leader=abs(C*Leader_pos(j)-X(i,j)); % Eq. (2.1)
                        X(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                        
                        
                    elseif p2>=0.5
                        %
                        
                        %% ------------- BHC LOCAL SEARCH-------------
                        
                        bw=0.001;   % N-operator probability
                        Beta=0.1;  % B-operator probability
                        New_X = X;
                        
                        rndindex = round((dim-1)*rand()+1);
%                           New_X(1,rndindex)=X(1,rndindex)+(2*rand()-1)*bw;
                          New_X(1,rndindex)=X(1,rndindex)+(2*rand()-1)*bw;
%                         New_X(rndindex, :) = X(rndindex, :)+bw*(2*rand()-1);
                        
                        for k=1:dim
                            if rand()<Beta
                                % New_X(1,k)=lb + (ub-lb)*rand();
                                New_X(k, :)=lb+(ub-lb)*rand();
                                
                            end
                        end;
                          FEs=FEs+2;
                        newfitness =fobj(New_X(i,:));
                        fitness=fobj(X(i,:));
                        if newfitness < fitness
                            X = New_X;
                            fitness= newfitness;
                        end % Finish updating the harmony memory
                             
             
                    end
                end
                
            elseif p>=0.5
                
                distance2Leader=abs(Leader_pos(j)-X(i,j));
                % Eq. (2.5)
                
                 %% ------------- standard spiral motions-------------
                 
                     X(i,j)= distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                
            end
            
        end
    end
    t=t+1;
    Convergence_curve(t)=Leader_score;
    [t Leader_score];
 end
% toc;


