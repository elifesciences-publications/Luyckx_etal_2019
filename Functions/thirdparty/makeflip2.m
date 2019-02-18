function [rr1, tr1] = makeflip2(ninputs,noutputs,nstim)

%  Input
base = randn(ninputs,1)';
shuf = reshape(Shuffle(1:ninputs),ninputs/2,2);
for n = 1:nstim
    flip = ones(1,length(base));
    flip(shuf(n,:)) = -1;
    base = base.*flip;
    rr1(:,n) = base;
end

% Output signal
base = -0.5+(rand(noutputs,1))';
shuf = Shuffle(1:noutputs);
for n = 1:nstim
    flip = ones(1,length(base));
    flip(shuf(n)) = -1;
    base = base.*flip;
    tr1(n,:) = base;
end

tr1 = tr1+0.5;