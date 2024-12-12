% - A Script to Analyze simulated data. % MODIFIED 10/13/23 to now account
% for taking the error at j-2,j+2 averages, becuase we are concerned about
% correlations in the noise due to teh PSF being ~2pixels. Test this and
% see if it makes a difference or not...
function [acorr,meansigsqr,plateau_value, plateau_data,sigmasquared, numSample ,meandY3point,noise_estimate_3point,noise_corr3point,noise_corr_abs3point,...
    sigma_3_twoThirds3pt,f, gof, slope, yintercept, Xintercept_seed_positions,Lp,Lplower,Lpupper,dYsqr_3point...
    ] = analyze_backbones(data)
import java.awt.Robot;
mouse=Robot;


%_ a function take in a backbone and get out the values of Lp


YpAll=data.Yp;; % - gets the Yp's

YpAllOrig=YpAll;


pixel_of_existence = 120; % based on pixel size thi should correspond to ~12um in real space
pixel_of_existenceYp=YpAll(pixel_of_existence,:);
num_nonNaN=6; % - number of consecuive measurements to consider this tube as existing at this length
% - now find first column where we get 4 non -NaN values to begin

for(h=1:size(YpAll,2))
    if(  ~isnan(pixel_of_existenceYp(h)) & ~isnan(pixel_of_existenceYp(h+1)) & ~isnan(pixel_of_existenceYp(h+2)) & ~isnan(pixel_of_existenceYp(h+3)))
        startingPoint=h;
        break; % - found starting point, get out of the loop
    end
end


startingPoint% - the frame where the MT is at least 12 um long



% - now to get the ending point

sizes=size(data.Yp);
maxTubeLen=sizes(1);
%- in cases where it doesnt get tobe 24 um long...

if ( maxTubeLen > 240 )
    endLength =240;
else
    endLength = sizes(1) %maxTubeLe if not up to 24 take what we got and take note of it!
end

% - end @ 24 um so look at 12-24 um only
endingYp=YpAll(endLength,:);


num_nonNaN=6; % - number of consecuive measurements to consider this tube as existing at this length
% - now find first column where we get 4 non -NaN values
for(h=startingPoint:size(YpAll,2))
    if(  ~isnan(endingYp(h)) & ~isnan(endingYp(h+1)) & ~isnan(endingYp(h+2)) & ~isnan(endingYp(h+3)) & ~isnan(endingYp(h+4)))
        EndingPoint=h;
        break; % - found ending point, get out of the loop
    end
end


EndingPoint % - the frame where the tube is 24um ( 240px) long.



Yp=YpAll(1:pixel_of_existence,startingPoint:EndingPoint); % - trim the data to what we want %-NOW Yp is from 12 to 24 um the frames
maxL=size(Yp,1); % MT lax length
winstart =1 ;  % whcih frame to start at
winend = size(Yp,2); % which frame to end at

ncomp = 400; % - how many frames over which to computer the acorr values. 400 is more than enough for our data.


% - calculate the  sigsqr values from the plateau of the acorr(MSD)
[acorr,meansigsqr,plateau_value, plateau_data,sigmasquared ] = length_amp_dynamics_corr_3pt(ncomp,maxL,1,winend,Yp,1);

for(k=size(acorr,1):-1:1)
    %- so k is iterating over positions

    k;
    max_val(k)=max(acorr(k,:));
    if(max_val(k)<=0 |isnan(max_val(k)))
        thalf(k)=0;
    else
        thalf(k) = max ( find( acorr(k,:)< max_val(k)/2) );
    end


end

%
%
numSample=(thalf(120)./size(Yp(:,winstart:winend),2))^-1;

%  - a way to measure the measurement noise, idea is that neighbors are
%  close by and all points have about the same noise levels then we look at
% [yj - (yj-2 + yj+2)/2]^2 

%%%%%%%%%%%%%%%%%%%%%%%%%


for(t=1:size(Yp,2)) % - loop over the length now..
    % - got the data now..
    ymeasure= Yp(:,t);
    pixel_distnace_for_error=2; %- here define the use of 2 pixel gap for error estimation
    for(K=1:pixel_of_existence) % - now find the differences


        if (K<pixel_distnace_for_error+1) % was K==1
            'one';
            dYsqr_3point(K,t) = [ ymeasure(K)- ymeasure(K+pixel_distnace_for_error)].^2;% was ymeasure(K+1)
        elseif (K>=pixel_of_existence-pixel_distnace_for_error   ) % - was K==pixel_of_existence
            'poe';
            dYsqr_3point(K,t) = [ ymeasure(K)- ymeasure(K-pixel_distnace_for_error)].^2;% was ymeasure(K-2)

        else
            'other';
            K;
            dYsqr_3point(K,t) =[ ymeasure(K)- ( ymeasure(K-pixel_distnace_for_error)+ymeasure(K+pixel_distnace_for_error) )/2].^2; % was ymeasure(K-2)+ymeasue(K+2)% CAN CHANGE THE POINTS WE USE HERE!!!
        end

    end
end

meandY3point = nanmean(dYsqr_3point,2);
noise_estimate_3point = meandY3point*(2/3);

noise_corr3point= (sigmasquared - noise_estimate_3point);
noise_corr_abs3point= abs(sigmasquared - noise_estimate_3point);


sigma_3_twoThirds3pt = (noise_corr3point).^(1/3);

% - now do the fit

Tfit = 1;% the fitting threshold on sig 2/3 v L plots to keep it linear generlaly between 0.8-1




X=[1:pixel_of_existence];
Xfi=[max(find(sigma_3_twoThirds3pt<Tfit)):numel(sigma_3_twoThirds3pt)-1]'; % - data for fitting
Yfi=sigma_3_twoThirds3pt(Xfi);
%- check in case data is bad
if( sum(isnan(Yfi)~=0)) % - if there are NaN values fit once we ge past them all
    'there are NaN!!!'


    % - so here try and find the and strip the out
    idxNotNan= find(isnan(Yfi)==0) ; % - here are the indicied of all the non -nans values so

    Yf=Yfi(idxNotNan);
    Xf=Xfi(idxNotNan);


else
    Xf=Xfi;
    Yf=Yfi;

end
Xf(1);
Xf(end);

clf;
scatter(X',sigma_3_twoThirds3pt,'x');
hold on


if (license('checkout','Curve_Fitting_Toolbox')==0) % - just to make sure we have a licence and don't go to sleep waitin for it
    A=0;
    while(A==0)
        A = license('checkout','Curve_Fitting_Toolbox')
        mouse.mouseMove(round(1000*rand),round(1000*rand));
    end
end

[f,gof]=fit(Xf,Yf,'poly1');
plot(f,Xf,Yf,'g');

slope=f.p1;
yintercept=f.p2;
Xintercept_seed_positions = -yintercept/slope
Lp = (108e-9) * ( 1/slope.^3/3)
ci = confint(f);
Lpupper = (108e-9) * ( 1/ci(2).^3/3)
Lplower = (108e-9) * ( 1/ci(1).^3/3)



end

