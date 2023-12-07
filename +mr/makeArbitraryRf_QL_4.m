function [rf, gz,gzr, delay] = makeArbitraryRf_QL(signal,flip, igSlider, varargin)
%makeArbitraryRf Create an RF pulse with the given pulse shape.
%   rf=makeArbitraryRf(singal, flip) Create RF pulse with complex signal
%   and given flip angle (in radians)
%
%   rf=makeArbitraryRf(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create block pulse with frequency offset and phase offset.
%
%   [rf, gz]=makeArbitraryRf(..., 'Bandwidth', bw, 'SliceThickness', st)
%   Create RF pulse and corresponding slice select gradient. The bandwidth
%   of the pulse must be given for the specified shape.
%
%   See also  Sequence.makeSincPulse, Sequence.addBlock

validPulseUses = mr.getSupportedRfUse();

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeArbitraryRf';
    
    % RF params
    addRequired(parser, 'signal', @isnumeric);
    addRequired(parser, 'flipAngle', @isnumeric);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addParamValue(parser, 'timeBwProduct', 0, @isnumeric);
    addParamValue(parser, 'bandwidth', 0, @isnumeric);
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    % Delay
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
end
parse(parser, signal, flip,varargin{:});
opt = parser.Results;

if opt.dwell==0
    opt.dwell=opt.system.rfRasterTime;
end

% Qiang Liu for dwell time:
opt.dwell=2*opt.dwell;

signal = signal./abs(sum(signal.*opt.dwell))*flip/(2*pi); % original!
%    array_scale=[1.068 1.261 1.122 1.261 1.068];
%   array_scale=[1 9/8 9/8 9/8 1]; % from maxim
array_scale=[1 1 1 1 1];
signal=signal*array_scale(igSlider);

%  signal = signal./sum(abs(signal.*opt.dwell))*flip/(2*pi);


N=  length(signal);
duration = N*opt.dwell;
opt.bandwidth= opt.timeBwProduct/duration; % QL
opt.duration=duration; % QL
t = ((1:N)-0.5)*opt.dwell;

rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.shape_dur=duration;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;
rf.delay = opt.delay;
if ~isempty(opt.use)
    rf.use=opt.use;
end
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end

if nargout>1
    assert(opt.sliceThickness > 0, 'SliceThickness must be provided');
    assert(opt.bandwidth > 0, 'Bandwidth of pulse must be provided');
    %  warning('FIXME: there are some potential issues with the bandwidth and related parameters, double check (e-mail communication)');
    if opt.maxGrad > 0
        opt.system.maxGrad = opt.maxGrad;
    end
    if opt.maxSlew > 0
        opt.system.maxSlew = opt.maxSlew;
    end
    
    BW = opt.bandwidth;
    if opt.timeBwProduct > 0
        BW = opt.timeBwProduct/duration;
    end
    
    amplitude = BW/opt.sliceThickness;
    area = amplitude*duration;
    gz = mr.makeTrapezoid('z', opt.system, 'flatTime', opt.duration, 'flatArea', area);
    %opt.centerpos=mr.calcRfCenter(rf); % QL
    opt.centerpos=0.5; % let us try 0.5 at first.
    
    
    %%  ql: The rewinder is changing too fast, too. Oct 15th
    opt.system.maxGrad = 0.8*opt.system.maxGrad;
    opt.system.maxSlew = 0.5*opt.system.maxSlew; %0.9
    gzr= mr.makeTrapezoid('z', opt.system, 'Area', -area*(1-opt.centerpos)-0.5*(gz.area-area)); %QL
    
    if rf.delay > gz.riseTime
        gz.delay = ceil((rf.delay - gz.riseTime)/opt.system.gradRasterTime)*opt.system.gradRasterTime; % round-up to gradient raster
    end
    if rf.delay < (gz.riseTime+gz.delay)
        rf.delay = gz.riseTime+gz.delay; % these are on the grad raster already which is coarser
    end
end

% v1.4 finally eliminates RF zerofilling
% if rf.ringdownTime > 0
%     tFill = (1:round(rf.ringdownTime/1e-6))*1e-6;  % Round to microsecond
%     rf.t = [rf.t rf.t(end)+tFill];
%     rf.signal = [rf.signal, zeros(size(tFill))];
% end
if rf.ringdownTime > 0 && nargout > 2
    delay=mr.makeDelay(mr.calcDuration(rf)+rf.ringdownTime);
end

end