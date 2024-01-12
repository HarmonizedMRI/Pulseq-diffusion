%% init the music variables

%a_init=440/((2^(1/12))^4); % we could retune the note table to avoid mechanical resonances, but it sounds so strange...
mrMusic.init;

%% melody definition

barDurationSeconds=1.8; 
timeSignature = 4/4;

% Root Beer Rag by Billy Joel
% 3-voice MRI version by Jochen Leupold
melody = {
    % intro bar 0
    [o/2 o/4 o/8 o/16 a1/32 h1/32], ...
    [o/2 o/4 o/8 o/16 a1/32 h1/32], ...
    [o]; ...
    % intro bar 1
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c1/2                                            h1/2                                           ]; ...
    % intro bar 2
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [a1/2                                            f1/2                                           ]; ...
    % intro bar 3
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [cb/2                                            hbb/2                                          ]; ...
    % intro bar 4
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [abb/2                                           fbb/2                                          ]; ...
    % bar 1
    [o/16 c1/16 d1/16 c1/16 d1/16      e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16      e1/8 c1/32 cis1/32], ...
    [o/8        g1/16 o/16  g1/32 o/32 g1/8 o/16  g1/16 o/16  g1/16  o/16 g1/32 o/32 g1/8 o/16         ], ...
    [cb/8       o/8         gbb/8        o/8      cb/8        o/8         gbb/8         o/8            ]; ...
    % bar 2
    [d1/16 o/16 ais1/16 o/16 ais1/16 o/16 ais1/16 o/16 a1/16 g1/16 f1/16 o/16 f1/8 f1/16 c1/32 d1/32], ...
    [f1/16 o/16 d2/16   o/16 d2/16   o/16 d2/16   o/16 c2/16 o/16  o/16  o/16 c2/8 c2/16 o/16       ], ...
    [aisb/8     o/8          aisb/8       o/8          fbb/8       o/8        fbb/8 o/8             ]; ...
    % bar 3
    [e1/16 c1/16 d1/16 c1/16 d1/16      e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32], ...
    [g1/16 o/16  g1/16 o/16  g1/32 o/32 g1/8 o/16  g1/16 o/16  g1/16 o/16  g1/32 o/32 g1/8 o/16    ], ...
    [cb/8        o/8         gbb/8        o/8      cb/8        o/8         gbb/8     o/8           ]; ...
    % bar 4
    [d1/16 o/16 d1/16 o/16 d1/16 o/16 d1/16 o/16 hb/16 ab/16 gb/16 o/16 gb/8 gb/16 c1/32 d1/32], ...
    [f1/16 o/16 f1/16 o/16 f1/16 o/16 f1/16 o/16 d1/16 o/16 o/16 o/16 d1/8 d1/16 o/16], ...
    [aisb/8 o/8 aisb/8 o/8 gbb/8 o/8 gbb/8 o/8]; ...
    % bar 5
    [e1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32], ...
    [g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16], ...
    [cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8]; ...
    % bar 6
    [d1/16 o/16 ais1/16 o/16 ais1/16 o/16 ais1/16 o/16 a1/16 g1/16 f1/16 o/16 f1/8 f1/16 a1/32 a1/32], ...
    [f1/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 c2/16 o/16 o/16 o/16 c2/8 c2/16 a1/16], ...
    [aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8]; ...
    % bar 7
    [a1/16 gis1/16 a1/16 h1/16 d2/16 c2/16 h1/16 a1/16 g1/16 fis1/16 g1/16 gis1/16 a1/16 g1/16 f1/16 e1/16], ...
    [a1/16 gis1/16 a1/16 h1/16 d2/16 c2/16 h1/16 a1/16 g1/16 fis1/16 g1/16 gis1/16 a1/16 g1/16 f1/16 e1/16], ...
    [fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8]; ...
    % bar 8
    [d1/16 cis1/16 d1/16 e1/16 g1/16 f1/16 e1/16 d1/16 c1/4 o/4], ...
    [d1/16 cis1/16 d1/16 e1/16 g1/16 f1/16 e1/16 d1/16 c1/4 o/4], ...
    [db/8 o/8 gb/8 o/8 cb/16 o/16 gb/16 o/16 cb/8 o/8]; ...
    % intro bar 1
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c1/2 h1/2]; ...
    % intro bar 2
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [a1/2 f1/2]; ...
    % intro bar 3
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [cb/2 hbb/2]; ...
    % intro bar 4
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [abb/2 fbb/2]; ...
    % bar 9
    [o/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32], ...
    [o/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16], ...
    [cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8]; ...
    % bar 10
    [d1/16 o/16 ais1/16 o/16 ais1/16 o/16 ais1/16 o/16 a1/16 g1/16 f1/16 o/16 f1/8 f1/16 c1/32 d1/32], ...
    [f1/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 c2/16 o/16 o/16 o/16 c2/8 c2/16 o/16], ...
    [aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8]; ...
    % bar 11
    [e1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32], ...
    [g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16], ...
    [cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8]; ...
    % bar 12
    [d1/16 o/16 d1/16 o/16 d1/16 o/16 d1/16 o/16 hb/16 ab/16 gb/16 o/16 gb/8 gb/16 c1/32 d1/32], ...
    [f1/16 o/16 f1/16 o/16 f1/16 o/16 f1/16 o/16 d1/16 o/16 o/16 o/16 d1/8 d1/16 o/16], ...
    [aisb/8 o/8 aisb/8 o/8 gbb/8 o/8 gbb/8 o/8]; ...
    % bar 13
    [e1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32], ...
    [g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16], ...
    [cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8]; ...
    % bar 14
    [d1/16 o/16 ais1/16 o/16 ais1/16 o/16 ais1/16 o/16 a1/16 g1/16 f1/16 o/16 f1/8 f1/16 a1/32 a1/32], ...
    [f1/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 c2/16 o/16 o/16 o/16 c2/8 c2/16 a1/16], ...
    [aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8]; ...
    % bar 15
    [a1/16 gis1/16 a1/16 h1/16 d2/16 c2/16 h1/16 a1/16 g1/16 fis1/16 g1/16 gis1/16 a1/16 g1/16 f1/16 e1/16], ...
    [a1/16 gis1/16 a1/16 h1/16 d2/16 c2/16 h1/16 a1/16 g1/16 fis1/16 g1/16 gis1/16 a1/16 g1/16 f1/16 e1/16], ...
    [fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8]; ...
    % bar 16
    [d1/16 cis1/16 d1/16 e1/16 g1/16 f1/16 e1/16 d1/16 c1/4 o/4], ...
    [d1/16 cis1/16 d1/16 e1/16 g1/16 f1/16 e1/16 d1/16 c1/4 o/4], ...
    [db/8 o/8 gb/8 o/8 cb/16 o/16 gb/16 o/16 cb/8 o/8]; ...
    % bar 17
    [hb/16 o/16 c1/16 o/16 cis1/16 o/16 d1/16 o/16 e1/16 o/16 d1/16 o/16 hb/16 o/16 gb/16 o/16], ...
    [h1/16 o/16 c2/16 o/16 cis2/16 o/16 d2/16 o/16 e2/16 o/16 d2/16 o/16 h1/16 o/16 g1/16 o/16], ...
    [gb/16 o/16 ab/16 o/16 aisb/16 o/16 hb/16 o/16 c1/16 o/16 hb/16 o/16 gb/16 o/16 eb/16 o/16]; ...
    % bar 18
    [e1/16 o/16 f1/16 o/16 fis1/16 o/16 g1/16 o/16 a1/16 o/16 g1/16 o/16 a1/16 o/16 g1/16 o/16], ...
    [e2/16 o/16 f2/16 o/16 fis2/16 o/16 g2/16 o/16 a2/16 o/16 g2/16 o/16 a2/16 o/16 g2/16 o/16], ...
    [cb/16 o/16 db/16 o/16 disb/16 o/16 eb/16 o/16 fb/16 o/16 eb/16 o/16 fb/16 o/16 eb/16 o/16]; ...
    % bar 19
    [hb/16 o/16 c1/16 o/16 cis1/16 o/16 d1/16 o/16 e1/16 o/16 d1/16 o/16 hb/16 o/16 gb/16 o/16], ...
    [h1/16 o/16 c2/16 o/16 cis2/16 o/16 d2/16 o/16 e2/16 o/16 d2/16 o/16 h1/16 o/16 g1/16 o/16], ...
    [gb/16 o/16 ab/16 o/16 aisb/16 o/16 hb/16 o/16 c1/16 o/16 hb/16 o/16 gb/16 o/16 eb/16 o/16]; ...
    % bar 20
    [ab/16 o/16 gb/16 gb/16 gb/16 o/16 f1/16 o/16 e1/16 e1/16 e1/16 o/16 c1/16 c1/16 o/16 o/16], ...
    [a1/16 o/16 g1/16 g1/16 g1/16 o/16 a1/16 o/16 g1/16 g1/16 g1/16 o/16 c2/16 c2/16 o/16 o/16], ...
    [fb/16 o/16 eb/16 eb/16 eb/16 o/16 db/16 o/16 cb/16 cb/16 cb/16 o/16 cb/16 cb/16 o/16 o/16]; ...
    % bar 21
    [o/16 g1/16 o/16 a1/16 o/16 ais1/16 o/16 h1/16 o/16 c2/16 o/16 h1/16 o/16 g1/16 o/16 e1/16], ...
    [o/16 h1/16 o/16 c2/16 o/16 cis2/16 o/16 d2/16 o/16 e2/16 o/16 d2/16 o/16 h1/16 o/16 g1/16], ...
    [gbb/16 hb/16 abb/16 c1/16 aisbb/16 cis1/16 hbb/16 d1/16 cb/16 e1/16 hbb/16 d1/16 gbb/16 hb/16 ebb/16 gb/16]; ...
    % bar 22
    [o/16 c2/16 o/16 d2/16 o/16 dis1/16 o/16 e2/16 o/16 f2/16 o/16 e2/16 o/16 f2/16 o/16 e2/16], ...
    [o/16 e2/16 o/16 f2/16 o/16 fis2/16 o/16 g2/16 o/16 a2/16 o/16 g2/16 o/16 a1/16 o/16 g2/16], ...
    [cb/16 e1/16 db/16 f1/16 disb/16 fis1/16 eb/16 g1/16 fb/16 a1/16 eb/16 g1/16 fb/16 a1/16 eb/16 g1/16]; ...
    % bar 24
    [o/16 a2/16 c1/16 a1/16 a2/16 fis2/16 c2/16 a1/16 g2/16 c2/16 g1/16 cis2/16 cis2/32 o/32 cis2/8 cis2/16], ...
    [o/16 a2/16 c1/16 a1/16 a2/16 fis2/16 c2/16 a1/16 g2/16 c2/16 g1/16 g2/16 g2/32 o/32 g2/8 g2/16], ...
    [fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8]; ...
    % bar X
    [fis1/16 c2/16 e2/16 c2/16 f1/16 c2/16 d2/16 d2/16 e1/4 e1/16 o/16 o/16 a1/32 h1/32], ...
    [e2/16 c2/16 e2/16 c2/16 e2/16 c2/16 d2/16 d2/16 c2/4 c2/16 o/16 o/16 a1/32 h1/32], ...
    [db/8 o/8 gb/8 o/8 cb/16 cb/16 gb/16 o/16 cb/8 o/8]; ...
    % intro bar 1
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c1/2 h1/2]; ...
    % intro bar 2
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [a1/2 f1/2]; ...
    % intro bar 3
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [cb/2 hbb/2]; ...
    % intro bar 4
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [abb/2 fbb/2]; ...
    % S4 bar 1
    [o/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/16 d2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/32 cis2/32], ...
    [o/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16], ...
    [cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8]; ...
    % S4 bar 2
    [d2/16 o/16 ais2/16 o/16 ais2/16 o/16 ais2/16 o/16 a2/16 g2/16 f2/16 o/16 f2/8 f2/16 c2/32 d2/32], ...
    [f2/16 o/16 d3/16 o/16 d3/16 o/16 d3/16 o/16 c3/16 o/16 o/16 o/16 c3/8 c3/16 o/16], ...
    [aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8]; ...
    % S4 bar 3
    [e2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/16 d2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/32 cis2/32], ...
    [g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16], ...
    [cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8]; ...
    % S4 bar 4
    [d2/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 h1/16 a1/16 g1/16 o/16 g1/8 g1/16 c2/32 d2/32], ...
    [f2/16 o/16 f2/16 o/16 f2/16 o/16 f2/16 o/16 d2/16 o/16 o/16 o/16 d2/8 d2/16 o/16], ...
    [aisb/8 o/8 aisb/8 o/8 gbb/8 o/8 gbb/8 o/8]; ...
    % S4 bar 5
    [e2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/16 d2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/32 cis2/32], ...
    [g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16], ...
    [cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8]; ...
    % S4 bar 6
    [d2/16 o/16 ais2/16 o/16 ais2/16 o/16 ais2/16 o/16 a2/16 g2/16 f2/16 o/16 f2/8 f2/16 a2/32 a2/32], ...
    [f2/16 o/16 d3/16 o/16 d3/16 o/16 d3/16 o/16 c3/16 o/16 o/16 o/16 c3/8 c3/16 a2/16], ...
    [aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8]; ...
    % S4 bar 7
    [a2/16 gis2/16 a2/16 h2/16 d3/16 c3/16 h2/16 a2/16 g2/16 fis2/16 g2/16 gis2/16 a2/16 g2/16 f2/16 e2/16], ...
    [a2/16 gis2/16 a2/16 h2/16 d3/16 c3/16 h2/16 a2/16 g2/16 fis2/16 g2/16 gis2/16 a2/16 g2/16 f2/16 e2/16], ...
    [fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8]; ...
    % S4 bar 8
    [d2/16 cis2/16 d2/16 e2/16 g2/16 g1/16 a1/16 h1/16 e1/4 o/4], ...
    [d2/16 cis2/16 d2/16 e2/16 g2/16 g1/16 a1/16 h1/16 c2/4 o/4], ...
    [db/8 o/8 gb/8 o/8 cb/16 o/16 gb/16 o/16 cb/8 o/8]; ...
    % S4 bar 9
    [a1/16 f1/16 c1/16 f1/16 g1/16 cis1/16 ab/16 cis1/16 f1/16 d1/16 ab/16 d1/16 dis1/16 c1/16 gb/16 c1/16], ...
    [a1/16 f1/16 c1/16 f1/16 g1/16 cis1/16 ab/16 cis1/16 f1/16 d1/16 ab/16 d1/16 dis1/16 c1/16 gb/16 c1/16], ...
    [fb/8 o/8 eb/8 o/8 db/8 o/8 cb/8 o/8]; ...
    % S4 bar 10
    [d1/16 aisb/16 fb/16 aisb/16 c1/16 ab/16 fb/16 ab/16 c1/16 hb/16 gb/16 fb/16 gb/16 aisb/16 gb/16 eb/16], ...
    [d1/16 aisb/16 fb/16 aisb/16 c1/16 ab/16 fb/16 ab/16 c1/16 hb/16 gb/16 fb/16 gb/16 ais1/16 gb/16 eb/16], ...
    [aisbb/8 o/8 abb/8 o/8 gbb/8 o/8 cb/8 o/8]; ...
    % S4 bar 11
    [o/16 a1/16 f1/16 c1/16 g1/16 cis1/16 ab/16 cis1/16 f1/16 d1/16 ab/16 d1/16 dis1/16 c1/16 gb/16 c1/16], ...
    [o/16 a1/16 f1/16 c1/16 g1/16 cis1/16 ab/16 cis1/16 f1/16 d1/16 ab/16 d1/16 dis1/16 c1/16 gb/16 c1/16], ...
    [fb/8 o/8 eb/8 o/8 db/8 o/8 cb/8 o/8]; ...
    % S4 bar 12
    [d1/16 aisb/8 d1/16 aisb/8 gb/8 ab/4 o/8 a1/32 h1/32 c2/16], ...
    [d1/16 f1/8 d1/16 f1/8 e1/8 f1/4 o/8 a1/32 h1/32 c2/16], ...
    [aisbb/8 o/8 cb/8 o/8 fb/16 fb/16 cb/16 o/16 fbb/8 o/8]; ...
    % S4 bar 13
    [e2/16 c2/16 g1/16 c2/16 d2/16 gis1/16 e1/16 gis1/16 c2/16 a1/16 e1/16 a1/16 ais1/16 g1/16 e1/16 g1/16], ...
    [e2/16 c2/16 g1/16 c2/16 d2/16 gis1/16 e1/16 gis1/16 c2/16 a1/16 e1/16 a1/16 ais1/16 g1/16 e1/16 g1/16], ...
    [c1/8 o/8 hb/8 o/8 ab/8 o/8 gb/8 o/8]; ...
    % S5 bar 1
    [a1/16 f1/16 c1/16 a1/16 g1/16 e1/16 c1/16 e1/16 fis1/16 d1/16 c1/16 d1/16 f1/16 d1/16 hb/16 d1/16], ...
    [a1/16 c1/16 f1/16 a1/16 g1/16 e1/16 c1/16 e1/16 fis1/16 d1/16 c1/16 d1/16 f1/16 d1/16 hb/16 d1/16], ...
    [fb/8 o/8 eb/8 o/8 db/8 o/8 gb/8 o/8]; ...
    % S5 bar 2
    [o/16 e2/16 c1/16 g1/16 d2/16 gis1/16 e1/16 gis1/16 c2/16 a1/16 e1/16 a1/16 ais1/16 g1/16 e1/16 g1/16], ...
    [o/16 e2/16 c1/16 g1/16 d2/16 gis1/16 e1/16 gis1/16 c2/16 a1/16 e1/16 a1/16 ais1/16 g1/16 e1/16 g1/16], ...
    [c1/8 o/8 hb/8 o/8 ab/8 o/8 gb/8 o/8]; ...
    % S5 bar 3
    [a1/16 c1/8 a1/16 c1/8 hb/8 e1/4 o/8 g1/32 a1/32 ais1/16], ...
    [a1/16 c2/8 a1/16 c2/8 h1/8 c2/4 o/8 g1/32 a1/32 ais1/16], ...
    [fb/8 o/8 gb/8 o/8 cb/16 cb/16 gb/16 o/16 cb/8 o/8]; ...
    % S5 bar 4
    [h1/16 g2/8 h1/16 g2/8 g2/16 h1/16 g1/8 a1/8 g1/8 f1/8], ...
    [h1/16 f2/8 h1/16 f2/8 f2/16 h1/16 g2/8 f2/8 e2/8 d2/8], ...
    [gb/8 o/8 db/8 o/8 gb/8 fb/8 eb/8 db/8]; ...
%     % s5 bar 5
%     [e2/32 g2/32 a2/16 c3/16 ais2/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 c2/16 a1/16 ais1/16], ...
%     [e2/32 g2/32 a2/16 c3/16 ais2/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 c2/16 a1/16 ais1/16], ...
%     [cb/8 o/8 o/4 o/2]; ...
    % S5 bar 5 end version
    [e2/16 a2/16 c3/16 ais2/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 c2/16 a1/16 ais1/16], ...
    [e2/16 a2/16 c3/16 ais2/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 c2/16 a1/16 ais1/16], ...
    [cb/8 o/8 o/4 o/2]; ...
    % S5 bar 6
    [h1/16 g2/8 h1/16 g2/8 g2/16 h1/16 g1/8 a1/8 g1/8 f1/8], ...
    [h1/16 f2/8 h1/16 f2/8 f2/16 h1/16 g2/8 f2/8 e2/8 d2/8], ...
    [gb/8 o/8 db/8 o/8 gb/8 fb/8 eb/8 db/8]; ...
    % S5 bar 7
    [h2/16 c1/8 h2/16 c1/8 h2/16 c1/8 h2/16 c1/8 h2/8 a1/16 ais1/16], ...
    [c3/16 d2/8 c3/16 d2/8 c3/16 d2/8 c3/16 d2/8 c3/8 a1/16 ais1/16], ...
    [cb/8 o/8 o/4 o/2]; ...
    % S5 bar 8
    [h1/16 g2/8 h1/16 g2/8 g2/16 h1/16 g1/8 a1/8 g1/8 f1/8], ...
    [h1/16 f2/8 h1/16 f2/8 f2/16 h1/16 g2/8 f2/8 e2/8 d2/8], ...
    [gb/8 o/8 hbb/8 o/8 db/8 o/8 eb/8 db/8]; ...
%     % S5 bar 8 end version
%     [h1/16 g2/8 h1/16 g2/8 g2/16 g1/16 g2/16 o/16 g2/8 g1/16 o/16 g1/8], ...
%     [h1/16 g2/8 h1/16 d2/8 d2/16 g1/16 d2/16 o/16 d2/8 g1/16 o/16 d2/8], ...
%     [gb/8 o/8 db/8 o/8 gb/8 fb/8 eb/8 db/8]; ...
%     % S5 bar 9
%     [e1/8 ais1/8 a1/8 gis1/8 g1/32 g1/32 g1/16 a1/16 a2/16 a2/16 a1/16 c2/16 a1/16], ...
%     [g1/8 e2/8 dis2/8 d2/8 cis2/32 cis2/32 cis2/16 a1/16 a2/16 a2/16 a1/16 c2/16 a1/16], ...
%     [cb/16 o/16 cb/16 o/16 hbb/8 aisbb/8 abb/8 o/8 o/8 abb/8]; ...
    % S5 bar 9 end version
    [c3/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 g1/16 gis1/16 a1/16 c2/16 d2/16 a1/16], ...
    [c3/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 g1/16 gis1/16 a1/16 c2/16 d2/16 a1/16], ...
    [cb/8 o/8 o/4 o/2]; ...
    % S5 bar 10
    [d2/16 dis2/16 c2/16 a1/16 a2/16 dis2/16 d2/16 c2/16 a1/16 g1/16 c2/16 ais1/16 a1/16 g1/16 e1/16 dis1/16], ...
    [d2/16 dis2/16 c2/16 a1/16 a2/16 dis2/16 d2/16 c2/16 a1/16 g1/16 c2/16 ais1/16 a1/16 g1/16 e1/16 dis1/16], ...
    [fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8]; ...
    % S5 bar 11
    [d1/16 dis1/16 d1/16 c1/16 ab/16 gb/16 fb/8 eb/4 eb/8 o/8], ...
    [d1/16 dis1/16 d1/16 c1/16 ab/16 gb/16 fb/8 gb/4 gb/8 o/8], ...
    [db/8 o/8 gb/8 o/8 cb/16 cb/16 gb/16 o/16 cb/8 o/8]; ...
    % intro bar 1
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c1/2 hb/2]; ...
    % intro bar 2
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [ab/2 fb/2]; ...
    % intro bar 3
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16], ...
    [cb/2 hbb/2]; ...
    % intro bar 4
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16], ...
    [abb/2 fbb/2]; ...
    % closing bar 1
    [c3/16 a2/16 f2/16 c2/16 c3/16 fis2/16 dis2/16 c2/16 g2/16 c2/16 g1/16 cis2/16 cis2/32 o/32 cis2/8 cis2/16], ...
    [c3/16 a2/16 f2/16 c2/16 c3/16 fis2/16 dis2/16 c2/16 g2/16 c2/16 g1/16 g2/16 g2/32 o/32 g2/8 g2/16], ...
    [fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8]; ...
    % closing bar 2
    [e2/16 c2/16 e2/16 c2/16 e2/16 c2/16 h1/8 g1/4 g2/8 o/8], ...
    [e2/16 c2/16 e2/16 c2/16 e2/16 c2/16 d2/8 e1/4 c3/8 o/8], ...
    [db/8 o/8 gb/8 f1/8 c1/16 c1/16 gb/16 o/16 cb/8 o/8]; ...
};    

%     % bar X
%     [], ...
%     [], ...
%     []; ...

% cat rootbeer.cpp | grep 'freq1.*n1\|[Aa]nfang\|[Tt]akt' | sed 's/freq.=\(.*\); n1=\([0-9]*\); for.*/\1\/\2/;s/^0/o/' | tr '\n\r' '  '
% a1/32 h1/32 //anfang1  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 2  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //anfang3  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 4  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //Takt 1  o/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32 //Takt2  d1/16 o/16 ais1/16 o/16 ais1/16 o/16 ais1/16 o/16 a1/16 g1/16 f1/16 o/16 f1/8 f1/16 c1/32 d1/32 //Takt 3  e1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32 //Takt 4  d1/16 o/16 d1/16 o/16 d1/16 o/16 d1/16 o/16 hb/16 ab/16 gb/16 o/16 gb/8 gb/16 c1/32 d1/32 //Takt 5  e1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32 //Takt 6  d1/16 o/16 ais1/16 o/16 ais1/16 o/16 ais1/16 o/16 a1/16 g1/16 f1/16 o/16 f1/8 f1/16 a1/32 a1/32 //Takt 7  a1/16 gis1/16 a1/16 h1/16 d2/16 c2/16 h1/16 a1/16 g1/16 fis1/16 g1/16 gis1/16 a1/16 g1/16 f1/16 e1/16 //Takt 8  d1/16 cis1/16 d1/16 e1/16 g1/16 f1/16 e1/16 d1/16 c1/4 o/4 //anfang1  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 2  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //anfang3  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 4  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //Takt9  o/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32 //Takt10  d1/16 o/16 ais1/16 o/16 ais1/16 o/16 ais1/16 o/16 a1/16 g1/16 f1/16 o/16 f1/8 f1/16 c1/32 d1/32 //Takt 11  e1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32 //Takt 12  d1/16 o/16 d1/16 o/16 d1/16 o/16 d1/16 o/16 hb/16 ab/16 gb/16 o/16 gb/8 gb/16 c1/32 d1/32 //Takt 13  e1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/16 d1/16 c1/16 d1/16 c1/16 d1/16 e1/8 c1/32 cis1/32 //Takt 14  d1/16 o/16 ais1/16 o/16 ais1/16 o/16 ais1/16 o/16 a1/16 g1/16 f1/16 o/16 f1/8 f1/16 a1/32 a1/32 //Takt 15  a1/16 gis1/16 a1/16 h1/16 d2/16 c2/16 h1/16 a1/16 g1/16 fis1/16 g1/16 gis1/16 a1/16 g1/16 f1/16 e1/16 //Takt 16  d1/16 cis1/16 d1/16 e1/16 g1/16 f1/16 e1/16 d1/16 c1/4 o/4 //Takt 17  hb/16 o/16 c1/16 o/16 cis1/16 o/16 d1/16 o/16 e1/16 o/16 d1/16 o/16 hb/16 o/16 gb/16 o/16 //Takt 18  e1/16 o/16 f1/16 o/16 fis1/16 o/16 g1/16 o/16 a1/16 o/16 g1/16 o/16 a1/16 o/16 g1/16 o/16 //Takt19  hb/16 o/16 c1/16 o/16 cis1/16 o/16 d1/16 o/16 e1/16 o/16 d1/16 o/16 hb/16 o/16 gb/16 o/16 //Takt 20  ab/16 o/16 gb/16 gb/16 gb/16 o/16 f1/16 o/16 e1/16 e1/16 e1/16 o/16 c1/16 c1/16 o/16 o/16 // Takt 21  o/16 g1/16 o/16 a1/16 o/16 ais1/16 o/16 h1/16 o/16 c2/16 o/16 h1/16 o/16 g1/16 o/16 e1/16 // Takt 22  o/16 c2/16 o/16 d2/16 o/16 dis1/16 o/16 e2/16 o/16 f2/16 o/16 e2/16 o/16 f2/16 o/16 e2/16 o/16 a2/16 c1/16 a1/16 a2/16 fis2/16 c2/16 a1/16 g2/16 c2/16 g1/16 cis2/16 cis2/32 o/32 cis2/8 cis2/16 //Takt 24  fis1/16 c2/16 e2/16 c2/16 f1/16 c2/16 d2/16 d2/16 e1/4 e1/16 o/16 o/16 a1/32 h1/32 //anfang1  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 2  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //anfang3  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 4  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //S4 Takt 1  o/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/16 d2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/32 cis2/32 //Takt2  d2/16 o/16 ais2/16 o/16 ais2/16 o/16 ais2/16 o/16 a2/16 g2/16 f2/16 o/16 f2/8 f2/16 c2/32 d2/32 //Takt 3  e2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/16 d2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/32 cis2/32 //Takt 4  d2/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 h1/16 a1/16 g1/16 o/16 g1/8 g1/16 c2/32 d2/32 //Takt 5  e2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/16 d2/16 c2/16 d2/16 c2/16 d2/16 e2/8 c2/32 cis2/32 //Takt 6  d2/16 o/16 ais2/16 o/16 ais2/16 o/16 ais2/16 o/16 a2/16 g2/16 f2/16 o/16 f2/8 f2/16 a2/32 a2/32 //Takt 7  a2/16 gis2/16 a2/16 h2/16 d3/16 c3/16 h2/16 a2/16 g2/16 fis2/16 g2/16 gis2/16 a2/16 g2/16 f2/16 e2/16 //Takt 8  d2/16 cis2/16 d2/16 e2/16 g2/16 g1/16 a1/16 h1/16 e1/4 o/4 //S4 Takt 9  a1/16 f1/16 c1/16 f1/16 g1/16 cis1/16 ab/16 cis1/16 f1/16 d1/16 ab/16 d1/16 dis1/16 c1/16 gb/16 c1/16 //S4 Takt 10  d1/16 aisb/16 fb/16 aisb/16 c1/16 ab/16 fb/16 ab/16 c1/16 hb/16 gb/16 fb/16 gb/16 aisb/16 gb/16 eb/16 //S4 Takt 11  o/16 a1/16 f1/16 c1/16 g1/16 cis1/16 ab/16 cis1/16 f1/16 d1/16 ab/16 d1/16 dis1/16 c1/16 gb/16 c1/16 //S4 Takt 12  d1/16 aisb/8 d1/16 aisb/8 gb/8 ab/4 o/8 a1/32 h1/32 c2/16 //S4 Takt 13  e2/16 c2/16 g1/16 c2/16 d2/16 gis1/16 e1/16 gis1/16 c2/16 a1/16 e1/16 a1/16 ais1/16 g1/16 e1/16 g1/16 //S5 Takt 1  a1/16 f1/16 c1/16 a1/16 g1/16 e1/16 c1/16 e1/16 fis1/16 d1/16 c1/16 d1/16 f1/16 d1/16 hb/16 d1/16 o/16 e2/16 c1/16 g1/16 d2/16 gis1/16 e1/16 gis1/16 c2/16 a1/16 e1/16 a1/16 ais1/16 g1/16 e1/16 g1/16 //S5 Takt 3  a1/16 c1/8 a1/16 c1/8 hb/8 e1/4 o/8 g1/32 a1/32 ais1/16 //S5 Takt 4  h1/16 g2/8 h1/16 g2/8 g2/16 h1/16 g1/8 a1/8 g1/8 f1/8 e2/32 g2/32 a2/16 c3/16 ais2/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 c2/16 a1/16 ais1/16 e2/16 //g2/32 a2/16 c3/16 ais2/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 c2/16 a1/16 ais1/16 //S5 Takt 6  h1/16 g2/8 h1/16 g2/8 g2/16 h1/16 g1/8 a1/8 g1/8 f1/8 h2/16 c1/8 h2/16 c1/8 h2/16 c1/8 h2/16 c1/8 h2/8 a1/16 ais1/16 //S5 Takt 8  h1/16 g2/8 h1/16 g2/8 g2/16 h1/16 g1/8 a1/8 g1/8 f1/8 h1/16 g2/8 h1/16 g2/8 g2/16 g1/16 g2/16 o/16 g2/8 g1/16 o/16 g1/8 //S5 Takt 9  /*e1/8 ais1/8 a1/8 gis1/8 g1/32 g1/32 g1/16 a1/16 a2/16 a2/16 a1/16 c2/16 a1/16 c3/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 g1/16 gis1/16 a1/16 c2/16 d2/16 a1/16 //S5 Takt 10  d2/16 dis2/16 c2/16 a1/16 a2/16 dis2/16 d2/16 c2/16 a1/16 g1/16 c2/16 ais1/16 a1/16 g1/16 e1/16 dis1/16 //S5 Takt 11  d1/16 dis1/16 d1/16 c1/16 ab/16 gb/16 fb/8 eb/4 eb/8 o/8 //anfang1  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 2  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //anfang3  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 4  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 c3/16 a2/16 f2/16 c2/16 c3/16 fis2/16 dis2/16 c2/16 g2/16 c2/16 g1/16 cis2/16 cis2/32 o/32 cis2/8 cis2/16 e2/16 c2/16 e2/16 c2/16 e2/16 c2/16 h1/8 //0/16 g1/4 g2/8
% a1/32 h1/32 //anfang1  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 2  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //anfang3  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 4  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //Takt 1  o/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 //Takt2  f1/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 c2/16 o/16 o/16 o/16 c2/8 c2/16 o/16 //Takt 3  g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 //Takt 4  f1/16 o/16 f1/16 o/16 f1/16 o/16 f1/16 o/16 d1/16 o/16 o/16 o/16 d1/8 d1/16 o/16 //Takt 5  g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 //Takt 6  f1/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 c2/16 o/16 o/16 o/16 c2/8 c2/16 a1/16 //Takt 7  a1/16 gis1/16 a1/16 h1/16 d2/16 c2/16 h1/16 a1/16 g1/16 fis1/16 g1/16 gis1/16 a1/16 g1/16 f1/16 e1/16 //Takt 8  d1/16 cis1/16 d1/16 e1/16 g1/16 f1/16 e1/16 d1/16 c1/4 o/4 //anfang1  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 2  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //anfang3  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 4  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //Takt9  o/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 //Takt10  f1/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 c2/16 o/16 o/16 o/16 c2/8 c2/16 o/16 //Takt 11  g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 //Takt 12  f1/16 o/16 f1/16 o/16 f1/16 o/16 f1/16 o/16 d1/16 o/16 o/16 o/16 d1/8 d1/16 o/16 //Takt 13  g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 g1/16 o/16 g1/16 o/16 g1/32 o/32 g1/8 o/16 //Takt 14  f1/16 o/16 d2/16 o/16 d2/16 o/16 d2/16 o/16 c2/16 o/16 o/16 o/16 c2/8 c2/16 a1/16 //Takt 15  a1/16 gis1/16 a1/16 h1/16 d2/16 c2/16 h1/16 a1/16 g1/16 fis1/16 g1/16 gis1/16 a1/16 g1/16 f1/16 e1/16 //Takt 16  d1/16 cis1/16 d1/16 e1/16 g1/16 f1/16 e1/16 d1/16 c1/4 o/4 //Takt 17  h1/16 o/16 c2/16 o/16 cis2/16 o/16 d2/16 o/16 e2/16 o/16 d2/16 o/16 h1/16 o/16 g1/16 o/16 //Takt 18  e2/16 o/16 f2/16 o/16 fis2/16 o/16 g2/16 o/16 a2/16 o/16 g2/16 o/16 a2/16 o/16 g2/16 o/16 //Takt19  h1/16 o/16 c2/16 o/16 cis2/16 o/16 d2/16 o/16 e2/16 o/16 d2/16 o/16 h1/16 o/16 g1/16 o/16 //Takt 20  a1/16 o/16 g1/16 g1/16 g1/16 o/16 a1/16 o/16 g1/16 g1/16 g1/16 o/16 c2/16 c2/16 o/16 o/16 // Takt 21  o/16 h1/16 o/16 c2/16 o/16 cis2/16 o/16 d2/16 o/16 e2/16 o/16 d2/16 o/16 h1/16 o/16 g1/16 // Takt 22  o/16 e2/16 o/16 f2/16 o/16 fis2/16 o/16 g2/16 o/16 a2/16 o/16 g2/16 o/16 a1/16 o/16 g2/16 o/16 a2/16 c1/16 a1/16 a2/16 fis2/16 c2/16 a1/16 g2/16 c2/16 g1/16 g2/16 g2/32 o/32 g2/8 g2/16 //Takt 24  e2/16 c2/16 e2/16 c2/16 e2/16 c2/16 d2/16 d2/16 c2/4 c2/16 o/16 o/16 a1/32 h1/32 //anfang1  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 2  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //anfang3  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 4  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //S4 Takt 1  o/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 //Takt2  f2/16 o/16 d3/16 o/16 d3/16 o/16 d3/16 o/16 c3/16 o/16 o/16 o/16 c3/8 c3/16 o/16 //Takt 3  g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 //Takt 4  f2/16 o/16 f2/16 o/16 f2/16 o/16 f2/16 o/16 d2/16 o/16 o/16 o/16 d2/8 d2/16 o/16 //Takt 5  g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 g2/16 o/16 g2/16 o/16 g2/32 o/32 g2/8 o/16 //Takt 6  f2/16 o/16 d3/16 o/16 d3/16 o/16 d3/16 o/16 c3/16 o/16 o/16 o/16 c3/8 c3/16 a2/16 //Takt 7  a2/16 gis2/16 a2/16 h2/16 d3/16 c3/16 h2/16 a2/16 g2/16 fis2/16 g2/16 gis2/16 a2/16 g2/16 f2/16 e2/16 //Takt 8  d2/16 cis2/16 d2/16 e2/16 g2/16 g1/16 a1/16 h1/16 c2/4 o/4 //S4 Takt 9  a1/16 f1/16 c1/16 f1/16 g1/16 cis1/16 ab/16 cis1/16 f1/16 d1/16 ab/16 d1/16 dis1/16 c1/16 gb/16 c1/16 //S4 Takt 10  d1/16 aisb/16 fb/16 aisb/16 c1/16 ab/16 fb/16 ab/16 c1/16 hb/16 gb/16 fb/16 gb/16 ais1/16 gb/16 eb/16 //S4 Takt 11  o/16 a1/16 f1/16 c1/16 g1/16 cis1/16 ab/16 cis1/16 f1/16 d1/16 ab/16 d1/16 dis1/16 c1/16 gb/16 c1/16 //S4 Takt 12  d1/16 f1/8 d1/16 f1/8 e1/8 f1/4 o/8 a1/32 h1/32 c2/16 //S4 Takt 13  e2/16 c2/16 g1/16 c2/16 d2/16 gis1/16 e1/16 gis1/16 c2/16 a1/16 e1/16 a1/16 ais1/16 g1/16 e1/16 g1/16 //S5 Takt 1  a1/16 c1/16 f1/16 a1/16 g1/16 e1/16 c1/16 e1/16 fis1/16 d1/16 c1/16 d1/16 f1/16 d1/16 hb/16 d1/16 o/16 e2/16 c1/16 g1/16 d2/16 gis1/16 e1/16 gis1/16 c2/16 a1/16 e1/16 a1/16 ais1/16 g1/16 e1/16 g1/16 //S5 Takt 3  a1/16 c2/8 a1/16 c2/8 h1/8 c2/4 o/8 g1/32 a1/32 ais1/16 //S5 Takt 4  h1/16 f2/8 h1/16 f2/8 f2/16 h1/16 g2/8 f2/8 e2/8 d2/8 e2/32 g2/32 a2/16 c3/16 ais2/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 c2/16 a1/16 ais1/16 e2/16 //g2/32 a2/16 c3/16 ais2/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 c2/16 a1/16 ais1/16 //S5 Takt 6  h1/16 f2/8 h1/16 f2/8 f2/16 h1/16 g2/8 f2/8 e2/8 d2/8 c3/16 d2/8 c3/16 d2/8 c3/16 d2/8 c3/16 d2/8 c3/8 a1/16 ais1/16 //S5 Takt 8  h1/16 f2/8 h1/16 f2/8 f2/16 h1/16 g2/8 f2/8 e2/8 d2/8 h1/16 g2/8 h1/16 d2/8 d2/16 g1/16 d2/16 o/16 d2/8 g1/16 o/16 d2/8 //S5 Takt 9  g1/8 e2/8 dis2/8 d2/8 cis2/32 cis2/32 cis2/16 a1/16 a2/16 a2/16 a1/16 c2/16 a1/16 c3/16 a2/16 g2/16 e2/16 dis2/16 d2/16 dis2/16 d2/16 c2/16 a1/16 g1/16 gis1/16 a1/16 c2/16 d2/16 a1/16 //S5 Takt 10  d2/16 dis2/16 c2/16 a1/16 a2/16 dis2/16 d2/16 c2/16 a1/16 g1/16 c2/16 ais1/16 a1/16 g1/16 e1/16 dis1/16 //S5 Takt 11  d1/16 dis1/16 d1/16 c1/16 ab/16 gb/16 fb/8 gb/4 gb/8 o/8 //anfang1  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 2  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 //anfang3  c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 c1/16 g1/16 c2/16 //Anfang 4  c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 f1/16 c2/16 c1/16 c3/16 a2/16 f2/16 c2/16 c3/16 fis2/16 dis2/16 c2/16 g2/16 c2/16 g1/16 g2/16 g2/32 o/32 g2/8 g2/16 e2/16 c2/16 e2/16 c2/16 e2/16 c2/16 d2/8 e1/4 c3/8
% o/32 o/32 //anfang1  c1/2 h1/2 //Anfang 2  a1/2 f1/2 //anfang3  cb/2 hbb/2 //Anfang 4  abb/2 fbb/2 //Takt 1  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt2  aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8 //Takt 3  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt 4  aisb/8 o/8 aisb/8 o/8 gbb/8 o/8 gbb/8 o/8 //Takt 5  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt 6  aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8 //Takt 7  fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8 //Takt 8  db/8 o/8 gb/8 o/8 cb/16 o/16 gb/16 o/16 cb/8 o/8 //anfang1  c1/2 h1/2 //Anfang 2  a1/2 f1/2 //anfang3  cb/2 hbb/2 //Anfang 4  abb/2 fbb/2 //Takt9  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt10  aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8 //Takt 11  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt 12  aisb/8 o/8 aisb/8 o/8 gbb/8 o/8 gbb/8 o/8 //Takt 13  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt 14  aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8 //Takt 15  fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8 //Takt 16  db/8 o/8 gb/8 o/8 cb/16 o/16 gb/16 o/16 cb/8 o/8 //Takt 17  gb/16 o/16 ab/16 o/16 aisb/16 o/16 hb/16 o/16 c1/16 o/16 hb/16 o/16 gb/16 o/16 eb/16 o/16 //Takt 18  cb/16 o/16 db/16 o/16 disb/16 o/16 eb/16 o/16 fb/16 o/16 eb/16 o/16 fb/16 o/16 eb/16 o/16 //Takt19  gb/16 o/16 ab/16 o/16 aisb/16 o/16 hb/16 o/16 c1/16 o/16 hb/16 o/16 gb/16 o/16 eb/16 o/16 //Takt 20  fb/16 o/16 eb/16 eb/16 eb/16 o/16 db/16 o/16 cb/16 cb/16 cb/16 o/16 cb/16 cb/16 o/16 o/16 // Takt 21  gbb/16 hb/16 abb/16 c1/16 aisbb/16 cis1/16 hbb/16 d1/16 cb/16 e1/16 hbb/16 d1/16 gbb/16 hb/16 ebb/16 gb/16 // Takt 22  cb/16 e1/16 db/16 f1/16 disb/16 fis1/16 eb/16 g1/16 fb/16 a1/16 eb/16 g1/16 fb/16 a1/16 eb/16 g1/16 fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8 //Takt 24  db/8 o/8 gb/8 o/8 cb/16 cb/16 gb/16 o/16 cb/8 o/8 //anfang1  c1/2 h1/2 //Anfang 2  a1/2 f1/2 //anfang3  cb/2 hbb/2 //Anfang 4  abb/2 fbb/2 //S4 Takt 1  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt2  aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8 //Takt 3  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt 4  aisb/8 o/8 aisb/8 o/8 gbb/8 o/8 gbb/8 o/8 //Takt 5  cb/8 o/8 gbb/8 o/8 cb/8 o/8 gbb/8 o/8 //Takt 6  aisb/8 o/8 aisb/8 o/8 fbb/8 o/8 fbb/8 o/8 //Takt 7  fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8 //Takt 8  db/8 o/8 gb/8 o/8 cb/16 o/16 gb/16 o/16 cb/8 o/8 //S4 Takt 9  fb/8 o/8 eb/8 o/8 db/8 o/8 cb/8 o/8 //S4 Takt 10  aisbb/8 o/8 abb/8 o/8 gbb/8 o/8 cb/8 o/8 //S4 Takt 11  fb/8 o/8 eb/8 o/8 db/8 o/8 cb/8 o/8 //S4 Takt 12  aisbb/8 o/8 cb/8 o/8 fb/16 fb/16 cb/16 o/16 fbb/8 o/8 //S4 Takt 13  c1/8 o/8 hb/8 o/8 ab/8 o/8 gb/8 o/8 //S5 Takt 1  fb/8 o/8 eb/8 o/8 db/8 o/8 gb/8 o/8 c1/8 o/8 hb/8 o/8 ab/8 o/8 gb/8 o/8 //S5 Takt 3  fb/8 o/8 gb/8 o/8 cb/16 cb/16 gb/16 o/16 cb/8 o/8 //S5 Takt 4  gb/8 o/8 db/8 o/8 gb/8 fb/8 eb/8 db/8 cb/8 o/8 o/4 o/2 cb/8 o/8 o/4 o/2 //S5 Takt 6  gb/8 o/8 db/8 o/8 gb/8 fb/8 eb/8 db/8 cb/8 o/8 o/4 o/2 //S5 Takt 8  gb/8 o/8 hbb/8 o/8 db/8 o/8 eb/8 db/8 gb/8 o/8 db/8 o/8 gb/8 fb/8 eb/8 db/8 //S5 Takt 9  cb/16 o/16 cb/16 o/16 hbb/8 aisbb/8 abb/8 o/8 o/8 abb/8 cb/8 o/8 o/4 o/2 //S5 Takt 10  //fbb/8 //0/8 //fisbb/8 //0/8 //gbb/8 //0/8 //abb/8 //0/8 fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8 //S5 Takt 11  db/8 o/8 gb/8 o/8 cb/16 cb/16 gb/16 o/16 cb/8 o/8 //anfang1  c1/2 hb/2 //Anfang 2  ab/2 fb/2 //anfang3  cb/2 hbb/2 //Anfang 4  abb/2 fbb/2 fb/8 o/8 fisb/8 o/8 gb/8 o/8 ab/8 o/8 db/8 o/8 gb/8 f1/8 c1/16 c1/16 gb/16 o/16 cb/8 o/8

%% convert to the channel-frequency table

[pitches, durations] = mrMusic.melodyToPitchesAndDurations(melody,'timeSignature',timeSignature);

%% Pulseq sequence 

% never use full gradient performance because it is almost impossible to
% avoid the mechanical resonances of the gradient system completely
sys = mr.opts('MaxGrad',18,'GradUnit','mT/m',...
    'MaxSlew',160,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6 ...
);  
seq=mr.Sequence(sys);      % Create a new sequence object

pulseqUseWave=false; % use the "UseWave" option with case as it is really demanding both on the Pulseq environment and the scanner (shape memory)

seq = mrMusic.musicToSequence(seq, pitches, durations, 'barDurationSeconds', barDurationSeconds, 'pulseqUseWave', pulseqUseWave);

%% output
if pulseqUseWave
    seq.setDefinition('Name', 'rootbeer');
    seq.write('rootbeer.seq');
else
    seq.setDefinition('Name', 'rootbeer1');
    seq.write('rootbeer1.seq');
end

return
%% play

seq.sound();

return
%% optional slow step for checking whether we are staying within slew rate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 
