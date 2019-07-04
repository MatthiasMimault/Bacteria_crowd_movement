function plotTime(Tspend, n, Nt)
%% Time estimation
% Tspend is toc, n is the current step and Nt the total number of steps
Tend = Tspend*(Nt-n)/(n+1);
date = datestr(datetime('now','TimeZone','local','Format',...
    'd-MMM-y HH:mm:ss Z')+seconds(Tend));
fprintf(['Estimated end : ' date ' (%g s remaining)\n'],Tend);
end

