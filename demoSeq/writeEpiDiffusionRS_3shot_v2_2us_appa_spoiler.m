% three-shot reference EPI for Pulseq-diffusion
% Qiang Liu
% qliu30@mgh.harvard.edu

clc;close all;clear all;
seq_file='epidiff_3_shot_ref_1p5mm_88sli_spoiler_noshift.seq';
% Set system limits
lims = mr.opts('MaxGrad',78,'GradUnit','mT/m',...
    'MaxSlew',150,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89);
lims1 = mr.opts('MaxGrad',68,'GradUnit','mT/m',...
    'MaxSlew',100,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89);  % gz fat

seq=mr.Sequence(lims);     % Create a new sequence object
fov=220e-3; Nx=146; Ny=Nx; % Define FOV and resolution
thickness=1.5e-3;            % slice thinckness
Nslices=88;
bFactor=0; % s/mm^2
TE=50e-3; % min is 65
recon=false; % plot the traj for check
RSegment=3;% QL for multishot EPI
Echotimeshift=1; % for multishot
ref_count=3; % QL set for reference scan

for ref=1:ref_count % ref 1 works as the dummy scan, and I also turn off the PE gradients to collect the reference

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
        deltaky=RSegment*deltak; %QL for R2 PE
        kWidth = Nx*deltak;

        % Phase blip in shortest possible time
        blip_dur = ceil(2*sqrt(deltaky/lims.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
        % the split code below fails if this really makes a trpezoid instead of a triangle...
        gy = mr.makeTrapezoid('y',lims,'Area',-deltaky,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center
        %gy = mr.makeTrapezoid('y',lims,'amplitude',deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);

        extra_area=blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
        gx = mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
        actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
        gx.amplitude=gx.amplitude/actual_area*kWidth;
        gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
        gx.flatArea = gx.amplitude*gx.flatTime;

        % calculate ADC
        adcDwellNyquist=deltak/gx.amplitude/ro_os;
        adcDwellNyquist=2e-06; % match the GE version

        % round-down dwell time to 100 ns
        adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
        adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
        % MZ: no idea, whether ceil,round or floor is better for the adcSamples...
        adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
        % realign the ADC with respect to the gradient
        time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
        adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
        % this rounding actually makes the sampling points on odd and even readouts
        % to appear misalligned. However, on the real hardware this misalignment is
        % much stronger anyways due to the grdient delays

        % FOV positioning requires alignment to grad. raster... -> TODO

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
        Ny_pre=round(Ny_pre/RSegment); % QL for R=2
        Ny_post=round(Ny/2+1); % PE lines after the k-space center including the central line
        Ny_post=round(Ny_post/RSegment);
        Ny_meas=Ny_pre+Ny_post;

        % Pre-phasing gradients
        gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2);
        gyPre = mr.makeTrapezoid('y',lims,'Area',(Ny_pre*deltaky-(Nmulti-1)*deltak));% QL I believe it should be deltak instead of deltaky
        [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
        % relax the PE prepahser to reduce stimulation
        gyPre = mr.makeTrapezoid('y',lims,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
        % gxPre's duration is always large because we use PF for PE, so for
        % ms EPI we are safe -- we don't need to notice this duration changing --QL
        gyPre.amplitude=gyPre.amplitude*pe_enable;


        if (ref==3)
            gyPre=mr.scaleGrad(gyPre,-1);
            gy_blipdownup=mr.scaleGrad(gy_blipdownup,-1);
            gy_blipdown=mr.scaleGrad(gy_blipdown,-1);
            gy_blipup=mr.scaleGrad(gy_blipup,-1);
        end


        % Calculate delay times
        durationToCenter = (Ny_pre+0.5)*mr.calcDuration(gx);
        rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);
        rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
        delayTE1=ceil((TE/2 - mr.calcDuration(rf,gz) - mr.calcDuration(gzReph) + rfCenterInclDelay - rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
        delayTE1=delayTE1-mr.calcDuration(gz180_crusher_1); %QL

        delayTE2tmp=ceil((TE/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
        assert(delayTE1>=0);
        gxPre.delay=0;
        gyPre.delay=0;
        delayTE2=delayTE2tmp-mr.calcDuration(gxPre,gyPre);
        delayTE2=delayTE2-mr.calcDuration(gz180_crusher_1); %QL
        [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
        assert(delayTE2>=0);

        % diffusion weithting calculation
        % delayTE2 is our window for small_delta
        % delayTE1+delayTE2-delayTE2 is our big delta
        % we anticipate that we will use the maximum gradient amplitude, so we need
        % to shorten delayTE2 by gmax/max_sr to accommodate the ramp down
        small_delta=delayTE2-ceil(lims.maxGrad/lims.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
        big_delta=delayTE1+mr.calcDuration(rf180,gz180);
        % we define bFactCalc function below to eventually calculate time-optimal
        % gradients. for now we just abuse it with g=1 to give us the coefficient
        g=sqrt(bFactor*1e6/bFactCalc(1,small_delta,big_delta)); % for now it looks too large!
        gr=ceil(g/lims.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
        gDiff=mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims);
        assert(mr.calcDuration(gDiff)<=delayTE1);
        assert(mr.calcDuration(gDiff)<=delayTE2);

        % Calculate the echo time shift for multishot EPI Qiang Liu
        actual_esp=gx.riseTime+gx.flatTime+gx.fallTime;
        TEShift=actual_esp./RSegment;
        TEShift=round(TEShift,5);
        dETS=mr.makeDelay(TEShift); % I want to try a very silly way...since we only have 4 shots Let's test first. QL

        %% TR calculation Qiang Liu
        TR=9000e-3; % ql follow Congyu's sequence /s
        delayTR = ceil((TR - Nslices*mr.calcDuration(gp_r) - Nslices*mr.calcDuration(gn_r) - Nslices*TE- Nslices*Ny_post*mr.calcDuration(gx))/seq.gradRasterTime)*seq.gradRasterTime;
        delayTR=delayTR/Nslices; % QL 0807
        delayTR=round(delayTR,3); % pulseq will allow delay time at power -3
        dTR=mr.makeDelay(delayTR);

        for islice_1=1:Nslices
            freqOffset_factor(islice_1)=islice_1-1-(Nslices-1)/2;
        end
        slic_indexS=[2:2:Nslices 1:2:Nslices];
        interleaved_freqOffset_factor=freqOffset_factor(slic_indexS);

        % Define sequence blocks
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
            seq.addBlock(mr.makeDelay(delayTE1),gDiff);
            if(~recon)
                seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
            end
            seq.addBlock(rf180,gz180); % QL: it used to be gz180n
            if(~recon)
                seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
            end

            seq.addBlock(mr.makeDelay(delayTE2),gDiff);
            % QL for no blip
            if ref ==1
                gy_blipup.last=0;
                gy_blipdown.first=0;
                gy_blipdownup.last=0;
                gy_blipdownup.first=0;
            end

            if (Echotimeshift)
                if (Nmulti ==1)
                    % seq.addBlock(dETS); % echotimeshift for multishot
                else
                    if (Nmulti ==2)
                        seq.addBlock(dETS);
                    else
                        if (Nmulti ==3)
                            seq.addBlock(dETS);
                            seq.addBlock(dETS);

                        end
                    end
                end
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
                %  gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                % QL: for odd pe lines, the k-space trajectory was not correct
                if (mod(Ny_meas,2)==1)
                    if (i<Ny_meas)
                        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    end
                else
                    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                end
                % end of this modfication
            end

            if (Echotimeshift)
                if (Nmulti ==1)
                    seq.addBlock(dETS); % echotimeshift for multishot
                    seq.addBlock(dETS);
                    seq.addBlock(dETS);
                else
                    if (Nmulti ==2)
                        seq.addBlock(dETS);
                        seq.addBlock(dETS);
                        %   seq.addBlock(dETS);
                    else
                        if (Nmulti ==3)
                            seq.addBlock(dETS);
                            %     seq.addBlock(dETS);

                        end
                    end
                end
            end

            seq.addBlock(dTR); % seperate
        end %slice loop

    end % for multishot

end % EPI REF LOOP
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
% trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
%
% % plot k-spaces
% figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
% hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
% figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
% axis('equal'); % enforce aspect ratio for the correct trajectory display
% hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
%% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov thickness*Nslices]);
seq.setDefinition('Name', 'epi-diff');
% [pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/MP_GPA_K2309_2250V_951A_AS82.asc'); % TERRA-XR
seq.write(seq_file);
%%
function b=bFactCalc(g, delta, DELTA)
sigma=1;
kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end
