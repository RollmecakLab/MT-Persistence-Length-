
% simulate a growing dynamic microtubule
%- currently set for 1 value of kappa and 2 differnet backgroudn levels

% - setup simulation parameters !!!!!!!!
delta_t = .050;  % Time steps in seconds  -> most data will be at 50ms exposure so set this to 50!
Nt = 200; % # of time points for now lets test 1000 points


gamma = 3e-3 ;clf

sigma_m =0;
lseg=1.08e-07;%


% setup MT image paramters  for a 50 ms exposure 
MTsig =  700 ; % mean Bg subtracted int of MT but BEFORE convolution, 
MTstd =  500 ; % Std of BG

BGsig = 0;%1000;
BGstd = 0;% 25;


CNR =  ( MTsig - BGsig)/sqrt( MTstd.^2 + BGstd.^2)  ;


PSFstd = 2.2 ; % std of the PSF in pixels measured from data @ 50 ms exposure and 108 nm pixels

ParentDataFolder='/Users/spectorjo/Desktop/DesktopJuly2024/Lp ms 2024 /MATLAB SCRIPTS FOR SUBMISSION/New Folder'


GrowthRate = 2.8 ; % gr in um/min
Lmin = 12;
Lmax = 24 ;

Total_time = ceil((Lmax-Lmin)./GrowthRate)
Nt=Total_time * (60) * (1/delta_t)  % - to get the num of frames





L=[Lmin,Lmax]*1e-6; 

%- here loop over Kappas

% fo rLp = 1 2 3 4 5 7 10 mm
kappas=[ (1e-3)*1.38064852e-23*305.15 (2e-3)*1.38064852e-23*305.15 (3e-3)*1.38064852e-23*305.15 (4e-3)*1.38064852e-23*305.15 (5e-3)*1.38064852e-23*305.15 (6e-3)*1.38064852e-23*305.15  (7e-3)*1.38064852e-23*305.15 (8e-3)*1.38064852e-23*305.15 (9e-3)*1.38064852e-23*305.15 (10e-3)*1.38064852e-23*305.15];


for(k=1:1)% - can change to numbe rof runs here
    kappa=kappas(k)
    for(runNum=1:1) %- can loop over kappas here
        SubDataFolder = [ ParentDataFolder,['/RunNum_',num2str(runNum)],['/kappa_',num2str(kappa)]]; % - make a sub folder for the lengths...
        mkdir(SubDataFolder);
        runNum
        kappa
        % - now simulate the tubes
        
        [x, y, ym, L, a, t] = generate_growing_tubule(L, lseg, delta_t, Nt, kappa, sigma_m);
        
        %-  ok now that we got that let's make some images!
        
        sizeX=500;
        sizeY=200;
        raw_image = zeros(sizeX,sizeY,numel(ym(1,:)));
        
        tubeStart= 100;% where to start the free part
        seedStart = tubeStart-10; % - lets add a seed for fun. let the seed be 1 um

        num_images = numel(ym(1,:)); % - the number of simulated images
        
        
        
        for (image_num = 1 :num_images)
            clear current_tube_image
            clear current_tube_positions
            ym_pixel{image_num}=(ym{image_num}./(108e-9)); % - gets the displacements in pixels 
            x_pixel{image_num}=1:numel(x{image_num}); % to get in to pixel form
            
            current_tube_image(:,1)=x_pixel{image_num};
            current_tube_image(:,2)=ym_pixel{image_num};
            
            current_tube_positions =  tubeStart+current_tube_image ;
            
            for(j=1:numel(x{image_num}))
                ypix = current_tube_positions(j,2); 
                yl  = floor(ypix); 
                yh = ceil(ypix);  
               
               
                raw_image(current_tube_positions(j,1),yl,image_num) = 1 -(ypix-yl);  
                raw_image(current_tube_positions(j,1),yh,image_num) = 1 - (yh-ypix);  
 
            end
            % - now add in the Seed....
            
            raw_image(seedStart:tubeStart,tubeStart,:)=1; % the seed is 10 pixels long here, can change this later if needed..
        end
        
        % - now create the IntScaledImage
        IntScaledImage=raw_image.*normrnd(MTsig,MTstd,sizeX,sizeY,num_images);
        ConvImage=imgaussfilt(IntScaledImage,PSFstd);
        clear IntScaledImage %- to try and save on memory.
 
        %_ add loop with noise here...
        
        BgExpMean = 964;
        bglevel=[0,1,1.5,2];
        for(bg=1:2)
            
            
            bg
            if(bg==1)
                BgStd=0
                Bgmean = 0; 
            elseif(bg==2)
            BgStd = 35.04;
            Bgmean = 964;
            else
                Bgmean=964;
                BgStd=35.04*2;
            end
            BgImageNoise = ConvImage+normrnd(Bgmean,BgStd*1,sizeX,sizeY,size(ConvImage,3));

            BgImageNoise=uint16(BgImageNoise);
            imgname=[SubDataFolder,...
                '/_BG_',num2str(Bgmean),'_BGstd_',num2str(BgStd*1),'_kappa',num2str(kappa),'_run_',num2str(runNum),'.tiff'];
            
            % - Try faster tiff wiritng here...
            fTIF=Fast_Tiff_Write(imgname);
                for(k=1:num_images)
            fTIF.WriteIMG(BgImageNoise(:,:,k));
            fracdone=round((k./num_images)*100)
                end
            fTIF.close;
            
        end
           
    end % - ends loop for a given kapp
end % - ends loop over  different kappa


