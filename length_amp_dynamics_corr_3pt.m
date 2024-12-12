
function [acorr,meansigsqr,plateau_value, plateau_data,sigmasquared ] = length_amp_dynamics_corr_3pt(ncomp,maxL,winstart,winend,Yp,draw)


clear acorr;
clear acorr2;
clear plateau_data;
clear plateau_value;

ncomp = 400; % - # of points to compare @ 50 ms = 20 seconds. 
for(times=1:1)
close all;
clear acorr,
clear plateau_value;
clear sigmasquared;
    
    winstart
    winend
    
for(Len=maxL:-1:1)
clf;
winstart;
winend;

    dataThisWindow = Yp(Len,[winstart:winend]); 
    
    ns= length(dataThisWindow);

    for(tstep=1:ncomp)
        
        idx = [tstep:ns];
        A2=dataThisWindow(idx); % - get the values
        A3=dataThisWindow(idx-tstep+1); % - get the values
        idx_good = find(~(isnan(A2)|isnan(A3))); % - remove the NaN values
        A2f = A2(idx_good); % remove NaN
        A3f = A3(idx_good); % - remove NaN
        acorr(Len,tstep)= mean((A2f-A3f).*(A2f-A3f)); % - get the acorr value
       
    end
   meansigsqr(Len) = mean(~isnan(dataThisWindow)); % - gets the mean of the fluctuations

% _ extra stuff for plateau value checking.....
 plateau_data(Len,[1:1+ncomp/2])  =  acorr(Len,[0.5*ncomp:ncomp]); % - the plateau at the  last1/2 of data

 plateau_value(Len) =  mean(plateau_data(Len,[1:1+ncomp/2])); % - the value of the plateau a thsi point

min_plat_data=min(plateau_data(Len,:));
max_plat_data=max(plateau_data(Len,:));

if(  (min_plat_data > plateau_value(Len)-.20*plateau_value(Len)) && (max_plat_data < plateau_value(Len) +.20*plateau_value(Len)))
    %- if that is true this is a good legit plateau add it to the sigma sqr

    sigmasquared(1,Len) = plateau_value(Len);
else
    
sigmasquared(1,Len) = NaN; 
end


if(draw==1)
subplot(4,1,1)
   plot(dataThisWindow,'k'); hold on
     xlabel('time');
     ylabel('perpendicular displacement this window');


subplot(4,1,2)
   plot(acorr(Len,1:end-2),'k'); hold on
   plot( (1:ncomp),ones(1,ncomp)*plateau_value(Len),'r');
    plot( (1:ncomp),1.15.*ones(1,ncomp)*plateau_value(Len),'g-');
     plot( (1:ncomp),0.85.*ones(1,ncomp)*plateau_value(Len),'g-');
     xlabel('position on tube');
     ylabel('autocorrelation value');

   subplot(4,1,3)
   plot(1:maxL,sigmasquared)
   xlabel('position on tube');
   ylabel('sigma squared');

   subplot(4,1,4)
   plot(1:maxL,sigmasquared.^(1/3),'b');
   xlabel('position on tube');
   ylabel('sigma squared^1/3');
   % here we will calculate the slope of sigma^2/3 vs L should have an
   % offset that is the measurement error, fit it back until 

   pause(.1)
end       

end


end % - this loop just keep it running over night so I don't lose the ML license

sigmasquared=(sigmasquared')/2; %  uncorrleatd points

% sigmaSquaredScaled = sigmasquared.*108e-9.*108e-9; % puts the correct units on sigma squared
% Xscaled = [1:size(sigmaSquaredScaled,2)].*108e-9;