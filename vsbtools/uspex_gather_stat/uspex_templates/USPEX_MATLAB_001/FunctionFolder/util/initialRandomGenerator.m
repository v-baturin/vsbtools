function initialRandomGenerator()

% used to initial the random number 


rand('state',sum(100*clock));
c = clock;
rand_init = round(c(5)+c(6));
for i = 1:rand_init
    dummy = rand;
    dummy = randperm(5);
end
