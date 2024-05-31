function Z = HybridRate4(state,dx,grid,mutrange,mutprob,mutdist0,mutdist,d)
% takes current state vector, and other model parameters, 
% and returns the vector of derivatives. No within-group games! 

x0 = .5; b0 = .1; sigmab = 2; % Ex1: x0=.5, b0=.1, sigmab=.3
b1 = .001;
pop = sum(state)*dx;    % total population
for i = 1:grid  % find total birth rate for each type
    x = i/grid;
    birthrate(i) = b0 * exp(-(x-x0)^2/(2*sigmab^2));
end
birthrate = birthrate .* state * dx;
deathrate = state * d * pop * dx; %total deathrate for each type
for i = 1:grid
    new(i) = 0;      
    for j = max(1,i-mutrange):min(grid,i+mutrange) % total birthrate due to j + mutation
        if j == i
            new(i) = new(i) + birthrate(i)*(1-mutprob + mutprob*mutdist0);
        else
            new(i) = new(i) + birthrate(j)*mutprob*mutdist(abs(i-j));
        end
    end
    Z(1,i) = (new(i) - deathrate(i))/dx;
    Z(2,i) = (new(i) + deathrate(i))/dx; % variance for sde
end