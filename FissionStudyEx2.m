% "Fission as a source of variation..." Example2
% along with HybridRate4.m with x0 = .2, b0 = .1, sigmab = 3
epsilon = 10^-10; e0 = .00005; s = .000025; y0 = .8;
f0 = .0005; Msigma = .01; mutprob = .1;
mu = .000; d = .001; r = 0;
% b=0; c=0;
FissBins = 20;
MaxG = 500; grid = 100; Xmin = 0; Xmax = 1; dx = (Xmax-Xmin)/grid;
Groups = 100; last = 100; next = 101;
for i = 1:grid
    M123(i) = i;
end
mutdist0 = normcdf(dx/2,0,Msigma) - normcdf(-dx/2,0,Msigma);
mutrange = ceil(4*Msigma/dx); % mutation distribution
for i = 1:mutrange
    mutdist(i) = normcdf((i+.5)*dx,0,Msigma) - normcdf((i-.5)*dx,0,Msigma);
end

state = zeros(MaxG,grid); alive = zeros(1,MaxG); 
for i = 1:Groups                        % initialize groups
    state(i,.45*grid:.55*grid) = 200; alive(i) = 1; pop(i) = sum(state(i,:))*dx;
end
T = 20000; dt = .1; t = 0; iter = 0;
while t < T
    t = t + dt; iter = iter + 1;
    if floor(iter/250) == iter/250
        k = iter/250; time(k) = t;
        intraGroupAve(k,:) = AveCoopLevel;  % E(Y_k|G_k)
        intraGroupVar(k,:) = VarCoopLevel;  % Var(Y_k|G_k) 
        groups(k) = Groups;
        GroupAve(k) = betweenGroupAve;    % E(Y_k)
        GroupVar(k) = betweenGroupVar;    % Var(E(Y_k|G_k))
        GroupVar2(k) = betweenGroupVar2;   %group-centric average
        Pop(k) = sum(pop)/Groups;   % Ave group size
        EVarCoopLevel(k) = sum(VarCoopLevel)/Groups; % E(Var(Y_k|G_k))
        [t/T, r, Groups, betweenGroupAve, betweenGroupVar2, EVarCoopLevel(k)] %E(Y), sd(E(Y|G)), E(Var(Y|G))
    end
    pop = zeros(1,MaxG); 
    AveCoopLevel = zeros(1,MaxG); VarCoopLevel = zeros(1,MaxG);
    for i = 1:last
        if alive(i) == 1
            pop(i) = sum(state(i,:))*dx; % Ave, var for each group
            AveCoopLevel(i) = sum(state(i,:).*M123)/(pop(i)*grid^2);
            VarCoopLevel(i) = sum(state(i,:) .* M123.^2)/(pop(i)*grid^3);
            VarCoopLevel(i) = VarCoopLevel(i) - AveCoopLevel(i)^2;
        end
    end
    totpop = sum(pop);
    betweenGroupAve = sum(AveCoopLevel.*pop)/totpop; % E(E(Y_k|G_k))
    betweenGroupVar = sum((AveCoopLevel.^2) .* pop)/totpop; 
    betweenGroupVar = betweenGroupVar - betweenGroupAve^2; % Var(E(Y_k|G_k))
    betweenGroupAve2 = sum(AveCoopLevel)/Groups;   %group-centric statistics
    betweenGroupVar2 = sum(AveCoopLevel.^2)/Groups;
    betweenGroupVar2 = betweenGroupVar2 - betweenGroupAve2^2;
    avepop = sum(pop)/Groups;
    
    for i = 1:last             % update group i state
        if alive(i) == 1 
            Z = HybridRate4(state(i,:),dx,grid,mutrange,mutprob,mutdist0,mutdist,d);
            state(i,:) = state(i,:) + Z(1,:)*dt;
            for j = 1:grid
                if state(i,j) < epsilon
                    state(i,j) = 0;
                end
            end  
            
            if rand < f0*pop(i)*dt        % does i fission?
                alive(next) = 1;
                for j = 1:FissBins
                    cat = sum(state(i,(j-1)*grid/FissBins + 1:j*grid/FissBins));
                    if cat*dx < .01  % if # in category is small, assume fission is uniform
                        u = .5;
                        for jj = (j-1)*grid/FissBins+1:j*grid/FissBins
                            state(next,jj) = state(i,jj)*u; 
                            state(i,jj) = state(i,jj)*(1-u);
                        end
                    else  % daughter group#1 is N(#/2, #/4 * r^2)
                        done = 0;
                        while done == 0
                            Norm = .5*cat*dx + .5*sqrt(cat*dx)*r*randn;
                            if Norm < cat*dx & Norm > 0
                                done = 1;
                            end
                        end
                        frac = Norm/(cat*dx);
                        for jj = (j-1)*grid/FissBins+1:j*grid/FissBins
                            state(next,jj) = state(i,jj)*frac;
                            state(i,jj) = state(i,jj)*(1-frac);
                        end
                    end
                end
                pop(i) = sum(state(i,:))*dx; pop(next) = sum(state(next,:))*dx;
                AveCoopLevel(i) = sum(state(i,:).*M123)/(pop(i)*grid^2);
                AveCoopLevel(next) = sum(state(next,:).*M123)/(pop(next)*grid^2);
                Groups = Groups + 1;    
                if next > last % change last
                    last = next;
                end
                found = 0; z = next;
                while found == 0  % change next
                    z = z + 1;
                    if alive(z) == 0
                        next = z; found = 1;
                    end
                end
            end
            
            y = sum(state(i,:).*M123)*dx/sum(state(i,:));  % does i die?
            if rand < (e0 + s*(y-y0)^2) * Groups 
                alive(i) = 0; Groups = Groups - 1;
                if i == last % change last
                    found = 0; z = i-1;
                    while found == 0
                        if alive(z) == 1
                            last = z; found = 1;
                        end
                        z = z - 1;
                    end
                end
                if i < next  % change next
                    next = i;
                end
            end
        end
    end
end

xx = size(intraGroupAve); % clean up for display
I = xx(1); J = xx(2); delta = 1/grid;
for i =1:I 
    for j = 1:J 
        if intraGroupAve(i,j) ==0 
            intraGroupAve(i,j) = NaN;
        end
    end
end 


