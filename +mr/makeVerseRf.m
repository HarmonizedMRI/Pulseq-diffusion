function [rfv, gzv] = makeVerseRf(rf,gz,varargin)
%makeVerseRf Create a VERSE RF pulse based on mintverse().
%
%   See also  mintverse(), http://mrsrl.stanford.edu/~brian/mintverse/

validPulseUses = mr.getSupportedRfUse();

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeArbitraryRf';
    
    % RF params
    addRequired(parser, 'rf', @isstruct);
    addRequired(parser, 'gz', @isstruct);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    %addParamValue(parser, 'freqOffset', 0, @isnumeric);
    %addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    %addParamValue(parser, 'timeBwProduct', 0, @isnumeric);
    %addParamValue(parser, 'bandwidth', 0, @isnumeric);
    % % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    %addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    % % Delay
    %addParamValue(parser, 'delay', 0, @isnumeric);
    %addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value
    % % whether it is a refocusing pulse (for k-space calculation)
    %addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
end
parse(parser, rf, gz, varargin{:});
opt = parser.Results;

% convert the input gz to samples by converting it to an extended gradient first 
[~,gze]=mr.splitGradientAt(gz,rf.delay-opt.system.gradRasterTime);
tos=(0.5:0.5:(round(rf.shape_dur/opt.system.rfRasterTime)-0.5))*opt.system.rfRasterTime;

rfos=[0 interp1(rf.t,rf.signal,tos,'linear') 0];
gzos=[0 interp1(gze.tt+gze.delay,gze.waveform,tos+rf.delay,'linear') 0];

rflim=max(abs(rfos))/2.2; % this is the current VERSE strength control
gzlim=max(abs(gzos))*1.0; % this is another current VERSE strength control 1.1

[rfvos,gzvos] = mintverse(rfos,gzos,opt.system.rfRasterTime/2,rflim,gzlim,opt.system.maxSlew/2,opt.system.rfRasterTime/2);

rfv=rf; 
rfv.signal=rfvos(2:2:end);
rfv.t=([1:length(rfv.signal)]-0.5)*opt.system.rfRasterTime;
rfv.shape_dur=length(rfv.signal)*opt.system.rfRasterTime;

gzvss=gzvos(11:20:end); % TODO FIXME correct subsumpling calculation for arbitrary dwell time ratios
gzvn=ceil(rfv.shape_dur/opt.system.gradRasterTime);
if length(gzvss)<gzvn
    gzvss(length(gzvss)+1:gzvn)=0; % pad with zeros
end

gzv=mr.makeArbitraryGrad(gz.channel, gzvss,'delay',rf.delay,'system',opt.system); 

end
