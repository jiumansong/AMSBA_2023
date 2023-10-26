% Main programs starts here
%2020.12.by SongJiuman
function [best,Convergence_curve]=AMSBA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
n=SearchAgents_no;  
A=rand(1,n)+ones(1,n);      % Loudness  (constant or decreasing)
r=rand(1,n);      % Pulse rate (constant or increasing)
% This frequency range determines the scalings
% You should change these values if necessary
Qmin=0;         % Frequency minimum
Qmax=2;         % Frequency maximum

% Dimension of the search variables
d=dim;           % Number of dimensions 
% Lower limit/bounds/ a vector
Lb=lb.*ones(1,d);
% Upper limit/bounds/ a vector
Ub=ub.*ones(1,d);
% Initializing arrays
Q=zeros(n,1);   % Frequency
v=zeros(n,d);   % Velocities
% Initialize the population/solutions
for i=1:n
  Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);   
  Fitness(i)=fobj(Sol(i,:));
end

[fmin,I]=min(Fitness);
best=Sol(I,:);     

Convergence_curve=[];
FEs=0;
t=1;
%分为三个阶段，第一个阶段用于探索、第二个阶段用于混合过度，第三个阶段用于开发
perc_iterat=0.2;
perc_iterbt=0.8;
delta_k=MaxFEs*perc_iterat;   
beta_k=MaxFEs*perc_iterbt;


% Main loop
while  FEs < MaxFEs
% for t=1:N_gen, 
% Loop over all bats/solutions
       for i=1:n,
          Q=Qmin*ones(1,d)+(Qmin-Qmax)*rand(1,d);
          w=abs(fmin/Fitness(i));
          v(i,:)=(w+0.1)*v(i,:)+rand()*(Sol(i,:)-best).*Q;
          S(i,:)=Sol(i,:)+v(i,:);

          Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
          if(FEs<delta_k)              
             if rand>r(i)
                    g=randperm(n);
                    g(g==i)=[];       
                    S(i,:)=Sol(i,:)+rand()*(Sol(g(1),:)-Sol(g(2),:))+rand()*(Sol(g(3),:)-Sol(g(4),:));                                       
              end
          elseif (FEs>=delta_k && FEs<beta_k) 
              if rand>r(i)                 
                         S(i,:)=best+A(i)/5*cauchy1(0,0.5,d)*rand();                     
              end
          else  
          
                   Fitness(i)=fobj(S(i,:));
                   [ SSFitness,SSPosiiton ] = MRM(S(i,:),Sol,n,dim,fobj );
                   if SSFitness<Fitness(i)
                      S(i,:)=SSPosiiton;
                      Fitness(i)=SSFitness;
                   end
                               
          end   
     % Evaluate new solutions
           S(i,:)=simplebounds(S(i,:),Lb,Ub);
            if FEs<MaxFEs
                FEs=FEs+1;
                Fnew=fobj(S(i,:));
               % Update if the solution improves, or not too loud
               if (Fnew<=Fitness(i))
                    Sol(i,:)=S(i,:);
                    Fitness(i)=Fnew;
               end
               if (Fnew<fmin) && (rand<A(i))
                   A(i) = 0.9*A(i);
                   r(i) = r(i)*(1-exp((-0.9)*t));
               end    
            else
                break;
            end

          % Update the current best solution
          if Fnew<=fmin
                best=S(i,:);
                fmin=Fnew;
          end
        end
        Convergence_curve(t)=fmin;
        t=t+1;
end

end
% Application of simple limits/bounds
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound vector
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bound vector 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
end