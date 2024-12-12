

%   Inputs   -      fname -> file name of movie
%                   fit_params.cr_start  -> [column, row] of point located before
%                                       start of microtubule.
%                             .cr_end    -> [column, row of point located beyond
%                                        end of microtubule.
%                             .contour_width  -> Width of contour in pixels
%                                                Tube must stay within contour
%                             .psf_width ->  Gaussian PSF width in pixels
%                             .pixel_size -> size of pixel in meters.
%                             .contour_extra_length -> Extra number of pixels
%                                           beyond fixed and free ends of
%                                           tube
%                             .Nframe_avg_growth -> Average growth over +-
%                             frames, e.g. if 1, then average over n-1, n
%                             and n+1
%
%
%   Output -
%               result.Yp(:,Nimg) -> column vector (per image) of
%                                   perpendicular displancement of backbone
%                                   as function of distance along backbone.
%                                   Undefined values are NaN
%               result.Ip(:,Nimg) -> column vector (per image) of
%                                   integrated intensity versus
%               result.l_fixed_end -> Nimg row vector of positions of fixed
%                                       end
%               result.l_free_end  -> Nimg row vector of positions of free
%                                       end

function result = fit_backbone_v3_4_f_real_data_for_fig_ch8mt1_06262019(fname, fit_params)
    
    % File Name of Movie
    if (nargin < 1)
        % fname = '909242018_ch3_9um_187181_5_fr1to100.tif';
        fname = 'example_data.tif';
    end
    
    % Start and end positions of tube
    if (nargin < 2)



    fit_params.cr_start=[58,324] % - for this file : 06292019_ch8_mt3
    fit_params.cr_end = [43,26]
    rs=[420,222]; % best guess for seed, not used at this time

    end
    
    % Width of contour (+/- in pixels)
    if ~isfield(fit_params,'contour_width')
        fit_params.contour_width = 10; % - might need to try different values here...
    end
    
    % Width of gaussian for PSF.
    if ~isfield(fit_params,'psf_width')
        fit_params.psf_width = 2.2; %- should be 1.9 make it 2 just to fit gooder
    end
    
    % Pixel size
    if ~isfield(fit_params,'pixel_size')
        fit_params.pixel_size = 108 * 1e-9; % Pixel size in m
    end
    
    % Pixel size
    if ~isfield(fit_params,'contour_extra_length')
        fit_params.contour_extra_length =10; % Number of extra pixels beyond ends of tube.
    end
    
    % Number of frames to average growth over
    if ~isfield(fit_params,'Nframe_avg_growth')
        fit_params.Nframe_avg_growth = 10; % Average growth over +- Nframe_avg_growth frames
    end


    [baseName, path]=uigetfile({'*.tif';'*.tiff'});
fullFileName = fullfile(path, baseName);

fname = fullFileName; % to make it easier than typing fullFileName each time...

  Movie=TIFFStack(fname);
    
    Imax = 1150;   % - looks good for new Darkfield data but might need adjusting
    Imin = 650;
    img=double(Movie); % - I've already loaded the movie elsewhere some image processing requires it be in doubles ormat so this is just a precaution
    % Initialize Coordiantes
    [Ny, Nx, Nimg] = size(img);
    Xv = ones(Ny,1) * [1:Nx];
    Yv = [1:Ny]' * ones(1,Nx);
    

    
    rb=fit_params.cr_start % - set up above
    rf=fit_params.cr_end
   % rs=[193,393] % for 12 um sim data with no noise location of the seed. don't need this now so just fill it with these numbers and do notihgn with it.

    Width = 40; % Width of contour (+/- in pixels to look for backbone - must be bigger than largest observable defleciton
    Length = round(norm(rf - rb))
    vn = rf - rb; vn = vn / norm(vn); % Normal vector
    vp = [-1 * vn(2), vn(1)] ; % Perpendicular vector
    nv = (Xv-rb(1))*vn(1) +   (Yv - rb(2))*vn(2); % Distance along tubule
    pv = (Xv-rb(1))*vp(1) +   (Yv - rb(2))*vp(2);
    
    % Now map pixels to positions along contour.
    for j=1:Length
        idx{j} = find( (j-0.5 <= nv) &( nv < j+0.5)&(abs(pv)<=Width) );
        %index of all pixels closest to position j along the tube.
        xp{j} = pv(idx{j}); % Distance of each of these pixels from the contour.
    end
    
    % Initialize Fit
    xc  = 0; % Initial guess for center position
    psf_width = fit_params.psf_width; % Rough guess for gaussian witdth
    Yp = zeros(Length,Nimg); % Displacement of each position
    Ip = zeros(Length,Nimg); % Intensity at each position.
    signalSTD = 20;
    
    
    meanTubeInt=1000; % measured from the movie, 
    % Initialize figure to show fit
    figure(1); clf;
    k = 1;
    img_handle = imagesc(img(:,:,k),[Imin Imax ]); axis image; hold on; colormap('gray'); % Plot image
    plot([rb(1), rf(1)], [rb(2), rf(2)], 'g'); % Add on contour director
    fit_handle = plot( 0, 0,'r'); % Initialize fitted contour plot.
    t_handle = title('0');
    
    for k=1:Nimg % Fit each image  %!!!!!!!!!!!!! CHANGE STARt POINT HERE!!!!!
        
        a = img(:,:,k); % Grab image data.
        xc = 0;
        
        for j = 1:Length; % Fit each point along contour from starting point to ending point
            %xp{j} is holding distance from contour of these pixels
            %a(idx{j} is holding the index of the pixels
            
            % - try  sliding the window  so we are not recenterg on the
            % entire length
            frameToFrameDisplacement = 15 ; % how far it can move between frames then
            %  'xp{j}'
            %     xp{j}
            %     'a(idx{j})'
            %       a(idx{j})
            pixelsWithContour = find(abs( xp{j} - xc ) < frameToFrameDisplacement);

             
             



                [xcn, resid, yfit, A] = recenterg(xp{j}(pixelsWithContour), a(idx{j}(pixelsWithContour)) , xc, psf_width,meanTubeInt); % Find center and average intensity
                
                Yp(j,k) = xcn; % Record center position
                Ip(j,k) = A;  % Record average
%             else
%                 xc=NaN;
%                 Yp(j,k) = NaN % Record center position
%                 Ip(j,k) = NaN
%             end
            %figure(2); clf;
            %plot(xp{j}, a(idx{j}),'bo'); hold on
            %plot(xp{j}, yfit,'k');
            %pause
            
            
            if isnan(xcn); else; xc = xcn ; end  % If center not defined, restore guess to 0.
            
        end
        
        % Now display result
        [xfit,yfit] = contour_to_xy(rb, vn, vp, Yp(:,k));
        set(img_handle, 'CData',a);
        set(fit_handle,'XData',xfit, 'YData',yfit);
        set(t_handle, 'String',sprintf('%d',k));
        drawnow;
    end
    
    %     % Initialize Fit
    %     xc  = 0; % Initial guess for center position
    %     psf_width = fit_params.psf_width; % Rough guess for gaussian witdth
    %     Yp = zeros(Length,Nimg); % Displacement of each position
    %     Ip = zeros(Length,Nimg); % Intensity at each position.
    %
    %     % Initialize figure to show fit
    %     figure(1); clf;
    %     k = 1;
    %     img_handle = imagesc(img(:,:,k),[2400 3561  ]); axis image; hold on; colormap('gray'); % Plot image
    %     plot([rb(1), rf(1)], [rb(2), rf(2)], 'g'); % Add on contour director
    %     fit_handle = plot( 0, 0,'r'); % Initialize fitted contour plot.
    %     t_handle = title('0');
    %
    %     for k=1:Nimg % Fit each image
    %
    %         a = img(:,:,k); % Grab image data.
    %
    %         for j = 1:Length; % Fit each point along contour
    %
    %             [xc, resid, yfit, A] = recenterg(xp{j}, a(idx{j}) , xc, psf_width); % Find center and average intensity
    %             Yp(j,k) = xc; % Record center position
    %             Ip(j,k) = A;  % Record average
    %
    %             %figure(2); clf;
    %             %plot(xp{j}, a(idx{j}),'bo'); hold on
    %             %plot(xp{j}, yfit,'k');
    %             %pause
    %
    %
    %             if isnan(xc); xc = 0; end  % If center not defined, restore guess to 0.
    %
    %         end
    %
    %         % Now display result
    %             [xfit,yfit] = contour_to_xy(rb, vn, vp, Yp(:,k));
    %             set(img_handle, 'CData',a);
    %             set(fit_handle,'XData',xfit, 'YData',yfit);
    %             set(t_handle, 'String',sprintf('%d',k));
    %             drawnow;
    %    end
    
    
    
    
    
    % Plot total fit
    figure(2); clf
    
    % Plot intensity along tube - i.e. kymograph.
    subplot(2,1,1);
    imagesc(Ip); colormap('gray'); colorbar; hold on;
    
    % Plot displacement along tube
    subplot(2,1,2);
    imagesc(Yp); colorbar; hold on
    
    
    subplot(2,1,1);
    % Find ends of tubule
     x_fixed = find_edge(Ip, 1,fit_params.contour_extra_length , fit_params.Nframe_avg_growth  )
     x_free = find_edge(Ip, -1,fit_params.contour_extra_length , fit_params.Nframe_avg_growth )
%     
%     
%     %- change this so that x_fixed is the rb point I have selected..
     %x_fixed(1:end)=rb(1);
%     
%     
%     
    % Mark on ends
    plot([1:length(x_fixed)], x_fixed,'r');
    plot([1:length(x_free)], x_free,'r');
    
    % Label plot
    xlabel('Image Number');
    ylabel('x (pixels)');
    title('Intensity')
    
    subplot(2,1,2);
    
    % Mark on ends
    plot([1:length(x_fixed)], x_fixed,'r');
    plot([1:length(x_free)], x_free,'r');
    
    % Label Plot
    xlabel('Image Number');
    ylabel('x (pixels)');
    title('Displacment (pixels)');
    
    % Now extract detailed backbone description.
    Nl = size(Yp,1)
    Ns = size(Yp,2)
    %pause;
    
    % First trim tube to length
    Ypnn = NaN * ones(size(Yp));
    for j=1:Ns
        j % - getting errors here if there is scum on the seed modifif 7/23/19 by jeff
        if (ceil(x_fixed(j))<0 && floor(x_free(j))<0)
         idx = [ceil(x_fixed(j-1)):floor(x_free(j-1))] 
         '1'
        elseif (floor(x_free(j))<0)
            idx = [ceil(x_fixed(j)):floor(x_free(j-1))]
            '2'
        elseif (ceil(x_fixed(j))<0 )
         idx = [ceil(x_fixed(j-1)):floor(x_free(j))] 
         '3'
        else
             idx = [ceil(x_fixed(j)):floor(x_free(j))];
             '4' ;
        end

        Ypnn(idx,j) = Yp(idx,j);
        %pause
    end
    
    % Now calculate means, rms, extremes
    Ymean = zeros(Nl, 1) * NaN;
    Ymax  = zeros(Nl, 1) * NaN;
    Ymin  = zeros(Nl, 1) * NaN;
    Yrms = zeros(Nl, 1) * NaN;
    Yp_ss = Yp;
    for j=1:Nl
        idx = find( ~isnan(Ypnn(j,:))) % Find all defined values at this length
        %-added by jeff for testing purposes..to see how many NAN we got
        numMeans(j)=numel(idx);
        % pause
        Ymean(j) = mean(Yp(j,idx));
        if (length(idx)>0)
            Ymax(j) = max(Yp(j,idx));
            Ymin(j) = min(Yp(j,idx));
        end
        Yp_ss(j,:) = Yp(j,:) - Ymean(j);
        Yrms(j) = std(Yp_ss(j,idx));
    end
    
    % Generate scaled and mean-subtracted displacements
    x_pix = [1:Nl]* fit_params.pixel_size;
    % And extract individual backbones
    for j=1:Ns
        idx = find(~isnan(Ypnn(:,j))); % Remove undetermined values
        tube.x{j} = x_pix(idx) ;
        tube.y{j} = Yp_ss(idx,j) * fit_params.pixel_size;
        tube.L(j) = x_free(j) * fit_params.pixel_size;
    end
    
    % Plot contours on image
    figure(1); clf;
    imagesc(img(:,:,Nimg)); axis image; hold on; colormap('gray'); % Plot image
    plot([rb(1), rf(1)], [rb(2), rf(2)], 'g'); % Plot path
    
    % Mark on Contour
    [xfit,yfit] = contour_to_xy(rb, vn, vp, Ymean);  plot(xfit,yfit,'r');
    [xfit,yfit] = contour_to_xy(rb, vn, vp, Ymean+Yrms);  plot(xfit,yfit,'r:');
    [xfit,yfit] = contour_to_xy(rb, vn, vp, Ymean-Yrms);  plot(xfit,yfit,'r:');
    [xfit,yfit] = contour_to_xy(rb, vn, vp, Ymax);  plot(xfit,yfit,'b');
    [xfit,yfit] = contour_to_xy(rb, vn, vp, Ymin);  plot(xfit,yfit,'b');
    
    
    % Plot contours on own
    figure(3); clf;
    subplot(2,1,1);
    plot(Ymean,'b'); hold on
    plot(Ymean+Yrms,'b:')';
    plot(Ymean-Yrms,'b:')';
    plot(Ymax,'r:')';
    plot(Ymin,'r:')';
    xlabel('Distance along contour (pixels)');
    ylabel('Displacment (pixels)');
    
    subplot(2,1,2);
    plot(Yrms,'bo');
    xlabel('Distance along contour (pixels)');
    ylabel('RMS Displacement (pixels)');
    
    % Save results
    result.Yp = Yp;
    result.Ip = Ip;
    result.x_fixed = x_fixed;
    result.x_free = x_free;
    result.tube = tube;
    result.x_pix = x_pix;
    result.sigma = Yrms * fit_params.pixel_size;
    %result.Xclick=Xclick;
    %result.Yclick=Yclick;
    result.rb=rb;
    result.rf=rf;
    result.rs=rs;
    result.numMeans=numMeans;
    

    
    
function [xc, resid, yfit, A] = recenterg(x, y, xc, width,meanTubeInt)
    % Fit y = c + m* x + A * exp(- (x-xc)^2 / (2 * width^2)) )
    %           c -> average background intensity
    %           m -> slope of background intensity
    %           A -> Amplitude of peak
    %           xc -> peak center
    %           width -> width of peak (essential psf)
    %
    % Input -> x,y -> vector of positions
    %           xc -> initial guess for peak position
    %           width ->  width of psf
    %           meanTubeInt --> mean intensity of the MT measured in teh
    %           movie, used for rejecting scum
    %    
    % Outputs ->
    %           xc -> refined peak center position
    %           resid -> rms of fitted intensities
    %           yfit -> fitted values of y
    %           A  -> fitted peak amplitude
    %
    %********MODIFIED 3/27/19by Jeff to do an actual gaussian fit after the
    %recentering procedure is finished to see if there is a difference
    
    
    
    % Pick initial fit values so loop excutes at least once
    
    Niter = 0;
    Niter_max = 50; % Limit to 50 iterations
    xc_cut = 0.01; % Limit to fit precision of 0.01 pixels
    xc_min = min(x) + width; % Limits on center position
    xc_max = max(x) - width;
    %     figure(5)
    %     plot(x,y)
    %     pause(0.1)
    while (Niter<=Niter_max)
        
        Niter = Niter + 1;
     
            % Perform linear regression
            % of y = c + m* x + A * f + (A * xc_err) * df/dxc
            %        c = cbest(1); m = cbest(2) ; A = cbest(3); A * xc_err =     cbest(4)
            
            f = exp( - (x - xc).^2 / (2*width^2))/ width ;% Peak function
%             figure(5)
%             plot(x,f,'r')
%             xc
%             figure(6)
%             plot(x,y) for testing purposes 

            M = [ones(size(x)), x, f, f.*(x-xc)/width^2];

            cbest = (M'*M) \ M' * y;
            rmsd = std(y - M*cbest); % RMS error
            xc_err = cbest(4)/cbest(3);
            % - HERE CHECK FOR AND REJECT SCUM
            outOfRange = sum( y  > meanTubeInt*1.2); % the number of brigher pixels
            if( outOfRange<4 )
                if (cbest(3)> 2*rmsd)&(abs(xc_err)<width)
                    % Peak signal to noise > 2 sigma
                    % Correction in peak position < width of peak
                    xc = xc + xc_err; % Update center position
                    yfit = M*cbest;
                    resid = sqrt(mean((y - yfit).^2));
                    A = cbest(3);
                
                % End fit if convert to within precision of xc_cut
                if (abs(xc_err)<xc_cut)
                    Niter = Niter_max + 1;
                end
            else
                %xc = xc + 0.3 * xc_err;
                Niter = Niter_max +1 ; % Abort fit
                xc = NaN;
                resid = NaN;
                yfit = NaN;
                A = cbest(3);
            end
            
        else
            Niter = Niter_max +1 ; % Abort fit
            xc = NaN;
            resid = NaN;
            yfit = NaN;
            A = cbest(3);
        end
        %****** START HERE GAUSSIAN FITTING PROCEDURE
        
        %
        %             fo=fitoptions('Method','NonlinearLeastSquares',...
        %                       'Lower',[0 0 0 0],...
        %                       'Upper',[Inf Inf Inf Inf],...
        %                       'StartPoint', [ Ai Bi Ci Di]); % - setup the fitting options
        %
        %                  ft=fittype( 'a+(b-a)*exp(-(x-c)*(x-c)/(2*d*d))','options',fo); % make the fit type
        %
        %                  [CURVE,GOF] = fit(X,Y,ft); % - do the fit!
        %
        %                 widths(k)=CURVE.d;
        %                 centers(k)=CURVE.c;
        %                 noise(k) = CURVE.a;
        %                 Amp(k) = CURVE.b;
        %                 GoF(k)= GOF.rsquare;
        %
        %                 Yp(j,k) = CURVE.c;
        %                 Ip(j,k)= CURVE.b;
        %
        %
        
        %FIXME
        % Could easily calculate moment of fitted peak to give an improved
        % estimate of psf width
        
        % FIXME
        % Could easily recurse until xc_err < 0.01 etc.
    end
    
function result = find_edge( Ip, dir, Nback, Navg)
    % Find end of microtubule
    % Inputs - Ip(:,Nimg) -> kymograph with each row corresponding to
    %                        intensity along tube
    % dir -> +1 -> start from j = 1 and find beginning of tubule
    %     -> -1 -> start from j = size(Ip,1) and find end of tubule
    % Nback  - default = 10 -> number of pixels at each end guaranteed to
    %                  just have background.
    % Navg -> +-Number of frames to average growth over.
    
    if (nargin <3)
        Nback = 10; % Number of background points on each end of contour
    end
    
    if (nargin <4)
        Navg = 5; % +-Number of frames to average over.
    end
    
    
    Nimg = size(Ip,2);
    Nl = size(Ip,1);
    
    idxs = 1  + (1-dir)/2*(Nl-1) + dir * ([1:Nback]-1); % Index of background
    kb = 1 + (1-dir)/2*(Nl-1) + dir * (Nback-1); % Starting point
    
    for j=1:Nimg
        
        % Average profile over frames.
        idx_avg = [max(j-Navg,1):min(Nimg, j+Navg)]; % Frames to average
        Iv = sort(Ip(:,idx_avg),2); % Sort values
        Iv = mean(Iv(:,[3:(length(idx_avg)-2)]),2); % Average but ditch max and min.
        
        % Set noise threshold at 5 sigma above background
        Ibs = mean(Iv(idxs));
        Ibss = std(Iv(idxs));
        Ithresh = Ibs + 2 * Ibss; % - that is pretty high change it to 2
        
        % Find foot of tubule
        kf = kb;
        while (kf<(Nl-5))&(kf>5)&(Iv(kf+dir)<Ithresh)
            kf=kf+dir;
        end
        
        % Find top - require that intensity of point be > next 3 points
        ks  = kf+dir;
        while (ks<(Nl-3)) & (ks>3) & (Iv(ks)< max(Iv(ks+dir*[1:3])))
            ks = ks+dir;
        end
        
        
        
        % Now interpolate to find 50% threshold
        % Fixme -> Should fit straight line to all points between 25% and
        % 75% of maximum.
        Imid = (Iv(kf) + Iv(ks))/2;
        kmid = kf;
        while (kmid < (Nl-3)) & (kmid > 2) & ( Iv(kmid+dir)<Imid)
            kmid = kmid + dir;
        end
        
        
        result(j) = kmid + dir * (Imid - Iv(kmid))/(Iv(kmid+dir) - Iv(kmid));
    end
    
    
function [x,y] = contour_to_xy(rb, vn, vp, Y)
    % Go from contour coordinates to xy coordinates
    % rb - column, row of beginnign of contour
    % vn -> normal vector (col,row) of contour
    % vp -> vector perpendicular to contour (col, row)
    % Y -> displacements perpendicular to contour.
    
    idx = [1:size(Y,1)]';
    x  = rb(1) + vn(1) * idx + vp(1) * Y;
    y = rb(2) + vn(2) * idx + vp(2) * Y;
    
    
function result = load_movie(fname)
    % Load tiff movie into memory as 3D array
    
    movie_info = imfinfo(fname);
    
    Nimg = length(movie_info);
    Nrow = movie_info(1).Height;
    Ncol = movie_info(1).Width;
    
    result(Nrow,Ncol,Nimg) = 0;
    
    for j=1:Nimg
        result(:,:,j) = double(imread(fname,j));
    end
    