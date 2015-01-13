
function [sa_N,f] = create_population(N,sig)
if nargin < 2
    sig = 1;
end

global rmax sigma s ex_cells
sigma = sig;
rmax = 50;
s = -5:0.1:5;
sa = -5:5;
sa_N = linspace(min(sa),max(sa),N);   % uniformly spaced preferred stimulus values

for nn = 1:N
for ii = 1:length(s)

    f(ii,nn) = rmax * exp(-.5*(((s(ii)-sa_N(nn))/sigma)^2));
    
end
end
% flat = sum(f,2);

% figure;
% hold on
% plot(s,f./rmax,'k','LineWidth',0.25);
% % plot(s,flat./(2*max(flat)),'k','LineWidth',2)
% xlim([min(s) max(s)])

% ex_cells = round(linspace(40*N/100,60*N/100,N/14));
ex_cells_zero = find(abs(sa_N)==nthroot(min(sa_N.^20),20));
% ex_cells = [ex_cells_zero-22:11:ex_cells_zero+22];
ex_cells = [ex_cells_zero-24:5:ex_cells_zero+25];

end