% fix the diffusion gradient problem
% add XYZ spoilers to the front of the sequence
% balance the spoilers by adding three spoilers after the fat sup pulse
% Oct 03 2023 Qiang Liu
% qliu30@mgh.harvard.edu

clc;close all;clear all;
seq_file='epidiff_R3_1p5_88sli_A.seq'
% Set system limits
B0field=2.89;% Prisma
lims = mr.opts('MaxGrad',78,'GradUnit','mT/m',...
    'MaxSlew',150,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0',B0field); 
lims1 = mr.opts('MaxGrad',68,'GradUnit','mT/m',...
    'MaxSlew',100,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', B0field);  % 50/130
lims2 = mr.opts('MaxGrad',68,'GradUnit','mT/m',...
    'MaxSlew',150,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', B0field);  % for diffusion gradients 78/120
diff_dir=[pwd '/diffusion_table/'];
table=xlsread([diff_dir 'Book3_30B.xlsx']);
table=table(1:30,:);
table=[0 0 0;0 0 0;0 0 0; table];
table(:,3)=-table(:,3); % Pulseq uses left-handed coordinate

% For IDEA sequences, only b0 uses crushers,  in our
% sequence, we introduce a factor to keep the echo time the same for b0/b1000
crusher_switch=0;

seq=mr.Sequence(lims);     % Create a new sequence object
fov=220e-3; Nx=146; Ny=Nx; % Define FOV and resolution
thickness=1.5e-3;            % slice thinckness
Nslices=88;
bFactor=1000.*ones(1,size(table,1)); % s/mm^2
bFactor(1:3)=0;
TE=58e-3; %
RSegment=1;% multishot factor
R=3;% In-plane accerelation factor
Echotimeshift=0; % for multishot
diffusion_count=size(table,1);
recon=false; % plot the traj for check
% ref 1 works as the dummy scan, and I also turn off the PE gradients to collect the reference
for ref=1:diffusion_count

    if ref ==1
        pe_enable=0;               % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
    else
        pe_enable=1;
    end

    for Nmulti=1:RSegment
        ro_os=1;                   % oversampling factor (in contrast to the product sequence we don't really need it)
        readoutTime=5.9e-4;        % this controls the readout bandwidth
        partFourierFactor=0.75;    % partial Fourier factor: 1: full sampling 0: start with ky=0

        tRFex=3e-3;
        tRFref=3e-3;

        % Create fat-sat pulse
        sat_ppm=-3.45;
        sat_freq=sat_ppm*1e-6*lims.B0*lims.gamma;
        rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims1,'Duration',8e-3,...
            'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
        rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase
        gz_fs = mr.makeTrapezoid('z',lims1,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

        % Spoiler in front of the sequence
        spoiler_amp=3*8*42.58*10e2;
        est_rise=500e-6; % ramp time 280 us
        est_flat=2500e-6; %duration 600 us

        gp_r=mr.makeTrapezoid('x','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gp_p=mr.makeTrapezoid('y','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gp_s=mr.makeTrapezoid('z','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',lims1);

        gn_r=mr.makeTrapezoid('x','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gn_p=mr.makeTrapezoid('y','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);
        gn_s=mr.makeTrapezoid('z','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',lims1);


        % Create 90 degree slice selection pulse and gradient
        [rf, gz, gzReph] = mr.makeSincPulse_QL_1(pi/2,'system',lims,'Duration',tRFex,...
            'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

        % Create 90 degree slice refocusing pulse and gradients
        [rf180, gz180] = mr.makeSincPulse(pi,'system',lims,'Duration',tRFref,...
            'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','refocusing');

        crusher_d=0.95e-3;

        gz180_crusher_1=mr.makeTrapezoid('z',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens
        gz180_crusher_2=mr.makeTrapezoid('y',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens
        gz180_crusher_3=mr.makeTrapezoid('x',lims,'Amplitude',19.05*42.58*10e2,'Duration',crusher_d); % that's what we used for Siemens


        % define the output trigger to play out with every slice excitatuion
        trig=mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'

        % Define other gradients and ADC events
        deltak=1/fov;
        deltaky=RSegment*R*deltak; % Rsegement*R
        kWidth = Nx*deltak;

        % Phase blip in shortest possible time
        blip_dur = ceil(2*sqrt(deltaky/lims.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
        % the split code below fails if this really makes a trpezoid instead of a triangle...
        gy = mr.makeTrapezoid('y',lims,'Area',-deltaky,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center

        % readout gradient is a truncated trapezoid with dead times at the beginnig
        % and at the end each equal to a half of blip_dur
        % the area between the blips should be defined by kWidth
        % we do a two-step calculation: we first increase the area assuming maximum
        % slewrate and then scale down the amlitude to fix the area
        extra_area=blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
        gx = mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
        actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
        gx.amplitude=gx.amplitude/actual_area*kWidth;
        gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
        gx.flatArea = gx.amplitude*gx.flatTime;

        % calculate ADC
        adcDwellNyquist=deltak/gx.amplitude/ro_os;
        adcDwellNyquist=2e-06; % this serves as oversampling, but our primary aim is to match GE's adc
        % round-down dwell time to 100 ns
        adcDwell=floor(adcDwellNyquist*1e7)*1e-7;

        adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
        % MZ: no idea, whether ceil,round or floor is better for the adcSamples...
        adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
        % realign the ADC with respect to the gradient
        time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
        adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us

        % split the blip into two halves and produce a combined synthetic gradient
        gy_parts = mr.splitGradientAt(gy, blip_dur/2, lims);
        [gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
        gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims);

        % pe_enable support
        gy_blipup.waveform=gy_blipup.waveform*pe_enable;
        gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
        gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;

        % phase encoding and partial Fourier

        Ny_pre=round((partFourierFactor-1/2)*Ny-1);  % PE steps prior to ky=0, excluding the central line
        Ny_pre=round(Ny_pre/RSegment/R);
        Ny_post=round(Ny/2+1); % PE lines after the k-space center including the central line
        Ny_post=round(Ny_post/RSegment/R);
        Ny_meas=Ny_pre+Ny_post;

        % Pre-phasing gradients
        gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2);
        gyPre = mr.makeTrapezoid('y',lims,'Area',(Ny_pre*deltaky-(Nmulti-1)*R*deltak));
        [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
        % relax the PE prepahser to reduce stimulation
        gyPre = mr.makeTrapezoid('y',lims,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
        gyPre.amplitude=gyPre.amplitude*pe_enable;

        %% for after the first echo QL

        gxPre_post = mr.makeTrapezoid('x',lims1,'Area',-gx.area/2);
        gyPre_post = mr.makeTrapezoid('y',lims1,'Area',(Ny_post-1)*deltaky);

        [gxPre_post,gyPre_post]=mr.align('right',gxPre_post,'left',gyPre_post);
        % relax the PE prepahser to reduce stimulation
        gyPre_post = mr.makeTrapezoid('y',lims1,'Area',gyPre_post.area,'Duration',mr.calcDuration(gxPre_post,gyPre_post));

        gyPre_post.amplitude=gyPre_post.amplitude*pe_enable;

        if (mod(Ny_meas,2)==0)
            gxPre_post.amplitude=-gxPre_post.amplitude;
        end

        if (ref==3)
            gyPre=mr.scaleGrad(gyPre,-1);
            gy_blipdownup=mr.scaleGrad(gy_blipdownup,-1);
            gy_blipdown=mr.scaleGrad(gy_blipdown,-1);
            gy_blipup=mr.scaleGrad(gy_blipup,-1);
            gyPre_post=mr.scaleGrad(gyPre_post,-1);
        end

        if (ref < 4) % need to adjust delta TE1 and delta TE2 for b0, bcuz of the crusher
            crusher_switch=1;
        else
            crusher_switch=0;
        end

        % Calculate delay times
        durationToCenter = (Ny_pre+0.5)*mr.calcDuration(gx);
        rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);
        rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
        delayTE1=ceil((TE/2 - mr.calcDuration(rf,gz) - mr.calcDuration(gzReph) + rfCenterInclDelay - rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
        delayTE1=delayTE1-crusher_switch*mr.calcDuration(gz180_crusher_1); %QL
        delayTE2tmp=ceil((TE/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
        assert(delayTE1>=0);
        gxPre.delay=0;
        gyPre.delay=0;
        delayTE2=delayTE2tmp-mr.calcDuration(gxPre,gyPre);
        delayTE2=delayTE2-crusher_switch*mr.calcDuration(gz180_crusher_1); %QL
        [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
        assert(delayTE2>=0);

        %% diffusion calculating
        small_delta=delayTE2-ceil(lims2.maxGrad/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
        big_delta=delayTE1+mr.calcDuration(rf180,gz180);
        g=sqrt(bFactor(1,ref)*1e6/bFactCalc(1,small_delta,big_delta)); % QL

        gr=ceil(g/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
        gDiff=mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims2);
  %% split the diffusion gradient to 3 axis
        
        g_x=g.*table(ref,1);
        g_y=g.*table(ref,2);
        g_z=g.*table(ref,3);
        
        if ((sum((table(ref,:)).^2)==0)||(abs(sum(table(ref,:)))==1)) % b=0 or dwi with diffusion gradient on one axis we keep using the older version
            
            g_xr=ceil(abs(g_x)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime; % QL: diffusion Oct 27
            g_yr=ceil(abs(g_y)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            g_zr=ceil(abs(g_z)/lims2.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
            gDiff_x=mr.makeTrapezoid('x','amplitude',g_x,'riseTime',g_xr,'flatTime',small_delta-g_xr,'system',lims2);
            gDiff_y=mr.makeTrapezoid('y','amplitude',g_y,'riseTime',g_yr,'flatTime',small_delta-g_yr,'system',lims2);
            gDiff_z=mr.makeTrapezoid('z','amplitude',g_z,'riseTime',g_zr,'flatTime',small_delta-g_zr,'system',lims2);
            duration_diff=max(max(mr.calcDuration(gDiff_x), mr.calcDuration(gDiff_y)), mr.calcDuration(gDiff_z));
            
        else % 3 axis, 'rotation'
            
            [azimuth,elevation,r] = cart2sph(g_x,g_y,g_z);
            polar= -(pi/2-elevation);
            
            Gr=mr.rotate('z',azimuth,mr.rotate('y',polar,gDiff));
            if size(Gr,2)==3
                gDiff_x=Gr{1,2};
                gDiff_y=Gr{1,3};
                gDiff_z=Gr{1,1};
            else
                if size(Gr,2)==2
                    diffusion_blank=find( table(ref,:)==0);
                    switch diffusion_blank
                        case 2
                            gDiff_x=Gr{1,2};
                            gDiff_z=Gr{1,1};
                            gDiff_y=gDiff; gDiff_y.channel='y'; gDiff_y.amplitude=0; gDiff_y.area=0; gDiff_y.flatArea=0;
                        case 1
                            gDiff_z=Gr{1,1};
                            gDiff_y=Gr{1,2};
                            gDiff_x=gDiff; gDiff_x.amplitude=0; gDiff_x.area=0; gDiff_x.flatArea=0;gDiff_x.channel='x';
                        case 3
                            gDiff_x=Gr{1,2};
                            gDiff_y=Gr{1,1};
                            gDiff_z=gDiff; gDiff_z.amplitude=0; gDiff_z.area=0; gDiff_z.flatArea=0;gDiff_z.channel='z';
                    end
                end
            end
            duration_diff=mr.calcDuration(gDiff);
        end
        assert(duration_diff<=delayTE1);
        assert(duration_diff<=delayTE2);

        %% Calculate the echo time shift for multishot EPI Qiang Liu
        actual_esp=gx.riseTime+gx.flatTime+gx.fallTime;
        TEShift=actual_esp/RSegment;
        TEShift=round(TEShift,5);
        TEShift_before_echo=(Nmulti-1)*TEShift;
        if TEShift_before_echo ==0
            TEShift_before_echo=0.00001; % apply the minimum duration for the no delay case
        end
        TEShift_after_echo=(RSegment-(Nmulti-1))*TEShift;
        dETS_before=mr.makeDelay(TEShift_before_echo);
        dETS_after=mr.makeDelay(TEShift_after_echo);
        %% TR calculation
        TR=11000e-3;
        TR=TR/Nslices;
        delayTR = ceil((TR - mr.calcDuration(gp_r) - mr.calcDuration(gn_r) - TE- Ny_post*mr.calcDuration(gx) - mr.calcDuration(gxPre_post)- rfCenterInclDelay)/seq.gradRasterTime)*seq.gradRasterTime;
        delayTR=round(delayTR,3); % pulseq will allow delay time at power -3
        dTR=mr.makeDelay(delayTR);
        %% interleaved
        for islice_1=1:Nslices
            freqOffset_factor(islice_1)=islice_1-1-(Nslices-1)/2;
        end
        slic_indexS=[2:2:Nslices 1:2:Nslices];
        interleaved_freqOffset_factor=freqOffset_factor(slic_indexS);

        %% Define sequence blocks
        for s=1:Nslices
            seq.addBlock(gp_r,gp_p,gp_s);
            seq.addBlock(rf_fs, gn_r,gn_p,gn_s);
            rf.freqOffset=gz.amplitude*thickness*interleaved_freqOffset_factor(s);
            %     rf.phaseOffset=-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
            rf180.freqOffset=gz180.amplitude*thickness*interleaved_freqOffset_factor(s);
            %     rf180.phaseOffset=pi/2-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase
            rf180.phaseOffset=pi/2;
            seq.addBlock(rf,gz,trig);
            seq.addBlock(gzReph); %rewinder QL

            seq.addBlock(mr.makeDelay(delayTE1),gDiff_x, gDiff_y, gDiff_z);
            if (ref < 4) % 3->4
                if(~recon)
                    seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                end
            end

            seq.addBlock(rf180,gz180);

            if (ref < 4)
                if(~recon)
                    seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                end
            end
            seq.addBlock(mr.makeDelay(delayTE2),gDiff_x, gDiff_y, gDiff_z);

            % QL for no blip
            if ref ==1
                gy_blipup.last=0;
                gy_blipdown.first=0;
                gy_blipdownup.last=0;
                gy_blipdownup.first=0;
            end

            if (Echotimeshift)
                seq.addBlock(dETS_before); % echotimeshift for multishot
            end
            seq.addBlock(gxPre,gyPre);
            for i=1:Ny_meas
                if i==1
                    seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                elseif i==Ny_meas
                    seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                else
                    seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                end
                % QL: for odd pe lines, the k-space trajectory was not correct
                if (mod(Ny_meas,2)==1)
                    if (i<Ny_meas)
                        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    end
                else
                    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                end
            end

            if (Echotimeshift)
                seq.addBlock(dETS_after); % echotimeshift for multishot
            end
            seq.addBlock(gxPre_post,gyPre_post);
            seq.addBlock(dTR); % seperate
        end %slice loop

    end % for multishot
    disp(['diffusion dir:', num2str(ref)])
    clear gDiff gDiff_x gDiff_y gDiff_z
end % EPI REF diffusion LOOP
%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% do some visualizations

% seq.plot();             % Plot sequence waveforms
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
% figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
% hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
% figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
% axis('equal'); % enforce aspect ratio for the correct trajectory display
% hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov thickness*Nslices]);
seq.setDefinition('Name', 'epi-diff');
seq.write(seq_file);


%%
function b=bFactCalc(g, delta, DELTA)
sigma=1
kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end
