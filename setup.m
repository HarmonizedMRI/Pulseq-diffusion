function setup()
    % Add Pulseq to the MATLAB path
    tmp=pwd;
    addpath(genpath(tmp));
    % Display a message indicating successful setup
    disp('Setup complete. Pulseq added to MATLAB path.');
    cd demoSeq/
    edit writeEpiDiffusionRS_multishot_dti_v2_appa_R3_1003.m
end