 

function [outputvalue,bnvalues] = cuckoosc4(n,y,xr,sd)
if nargin<1
    n = 25;
end
q = size(xr);
%tmax = 2*sd*sqrt(2*log(q(2)));
number_of_solution = 4;
Lb = [ 0,0,.1,.1];
Ub = [1000,1,20,4];

for i=1:n
    nest(i,:) = Lb + (Ub - Lb).*rand(size(Lb));
end

pa = .25;           %discovery rate

fitness = 10^10.*ones(n,1);
[fmin,bestnest,nest,fitness] = get_best_nest(nest,nest,fitness,xr,y);

N_iterTotal = 40; %total number of iterations

N_iter = 0;

for iter=1:N_iterTotal
    
    new_nest = get_cuckoos(nest,bestnest,Lb,Ub);
    [fnew,best,nest,fitness] = get_best_nest(nest,new_nest,fitness,xr,y);
    
    N_iter = N_iter + n;
    
    new_nest = empty_nests(nest,Lb,Ub,pa);      %discovery and Randomization
    [fnew,best,nest,fitness] = get_best_nest(nest,new_nest,fitness,xr,y);
    
    N_iter = N_iter + n;
    
    if fnew < fmin
        fmin = fnew;
        bestnest = best;
    end
    
    bnvalues(iter,:) = bestnest;
end
 u = bestnest
    mxr = [];
            q = size(xr)
            thr = u(1);
            k = u(2);
            d = u(3);
            n = u(4);
             for m=1:q(2)
                if xr(1,m)>thr
                    mxr(1,m) = xr(1,m) -(((.5)*(thr^d)*(k))/(xr(1,m)^(d-1))) + (k-1)*thr;
                elseif abs(xr(1,m)) <= thr
                    mxr(1,m) = ((0.5)*(k*((abs(xr(1,m)))^n))*sign(xr(1,m)))/(thr^(n-1));
                elseif xr(1,m) < -thr
                    mxr(1,m) =xr(1,m) + (((.5)*((-thr)^d)*k)/(xr(1,m)^(d-1))) - (k-1)*thr;
                end
             end
            outputvalue = mxr;

%subfunctions

%cuckoos using levy flight
function nest = get_cuckoos(nest,best,Lb,Ub)
n = size(nest,1);
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n
    s=nest(j,:);
    u =randn(size(s))*sigma;
    v =randn(size(s));
    step = u./abs(v).^(1/beta);
    
    stepsize = .01*step.*(s - best);
                                                %stepsize = 1
    s = s + randn(size(s));

    nest(j,:) = SimpleBounds(s,Lb,Ub);
end

%current best nest

function [fmin,best,nest,fitness] = get_best_nest(nest,newnest,fitness,xr,y)

for j = 1:size(nest,1)
    fnew = fobj(newnest(j,:),xr,y);
    if fnew < fitness(j)
        fitness(j) = fnew;
        nest(j,:) = newnest(j,:);
    end
end
%current bests
[fmin,K] = min(fitness);
best = nest(K,:);

%replacing some nests by constructing new solutions/nests
function new_nest = empty_nests(nest,Lb,Ub,pa)
n = size(nest,1);

k = rand(size(nest)) > pa;
stepsize = rand*(nest(randperm(n),:) - nest(randperm(n),:));
 new_nest = nest + stepsize.*k ;
 
 for j=1:size(new_nest)
     s = new_nest(j,:);
     new_nest(j,:) = SimpleBounds(s,Lb,Ub);
 end
 
 %function for simple bounds
    function s  = SimpleBounds(s,Lb,Ub)
        ns_temp = s;
        I = s<Lb;
        ns_temp(I) = Lb(I);
        
        J = s>Ub;
        ns_temp(J) = Ub(J);
        
        s = ns_temp;
        
 %objective Function
 
        function z = fobj(u,xr,y)
            mxr = [];
            q = size(xr);
            thr = u(1);
            k = u(2);
            d = u(3);
            n = u(4);
            for m=1:q(2)
                if xr(1,m)>thr
                    mxr(1,m) = xr(1,m) - (((.5)*(thr^d)*(k))/(xr(1,m)^(d-1))) + (k-1)*thr;
                elseif abs(xr(1,m)) <= thr
                    mxr(1,m) = ((0.5)*(k*((abs(xr(1,m)))^n))*sign(xr(1,m)))/(thr^(n-1));
                elseif xr(1,m) < -thr
                    mxr(1,m) =xr(1,m) + (((.5)*((-thr)^d)*k)/(xr(1,m)^(d-1))) - (k-1)*thr;
                end
            end
         z = (.5)*sum((y-mxr).^2);
         