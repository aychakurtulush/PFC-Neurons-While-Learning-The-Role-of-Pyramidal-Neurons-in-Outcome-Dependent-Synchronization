% Load data
% Behavior are 6 colums data:
% 1st: trial start (ms)
% 2nd: Trial end  (ms)
% 3rd: trial rule (1: go right, 2: go the lit arm; 
%      3: go left, 4: go to dark arm)
% 4th: correct (1) or incorrect (0) trial
% 5th: the animal went right (0) or left (1)
% 6th: light location (switched on randomly at the beginning of each trial, 
%      even for left/right tasks). 1 or 0 (right or left)

Behav{1} = load('Data\150628\150628_Behavior.txt');
Behav{2} = load('Data\150629\150629_Behavior.txt');
Behav{3} = load('Data\150630\150630_Behavior.txt');
Behav{4} = load('Data\150701\150701_Behavior.txt');
%Behav{5} = load('Data\150704\150704_Behavior.txt');
Behav{5} = load('Data\150705\150705_Behavior.txt');
Behav{6} = load('Data\150706\150706_Behavior.txt');
Behav{7} = load('Data\150707\150707_Behavior.txt');
% Behav{8} = load('Data\150708\150708_Behavior.txt');
Behav{8} = load('Data\150711\150711_Behavior.txt');
Behav{9} = load('Data\150712\150712_Behavior.txt');
Behav{10} = load('Data\150713\150713_Behavior.txt');
Behav{11} = load('Data\150714\150714_Behavior.txt');
Behav{12} = load('Data\150715\150715_Behavior.txt');

Behav{13} = load('Data\181011\181011_Behavior.txt');
Behav{14} = load('Data\181012\181012_Behavior.txt');
Behav{15} = load('Data\181014\181014_Behavior.txt');
Behav{16} = load('Data\181017\181017_Behavior.txt');
Behav{17} = load('Data\181018\181018_Behavior.txt');
Behav{18} = load('Data\181019\181019_Behavior.txt');
Behav{19} = load('Data\181020\181020_Behavior.txt');
Behav{20} = load('Data\181021\181021_Behavior.txt');
Behav{21} = load('Data\181024\181024_Behavior.txt');
Behav{22} = load('Data\181025\181025_Behavior.txt');
Behav{23} = load('Data\181026\181026_Behavior.txt');
Behav{24} = load('Data\181027\181027_Behavior.txt');
Behav{25} = load('Data\181102\181102_Behavior.txt');
Behav{26} = load('Data\181103\181103_Behavior.txt');
 
Behav{27} = load('Data\190213\190213_Behavior.txt');
Behav{28} = load('Data\190214\190214_Behavior.txt');
Behav{29} = load('Data\190224\190224_Behavior.txt');
Behav{30} = load('Data\190226\190226_Behavior.txt');
Behav{31} = load('Data\190227\190227_Behavior.txt');
Behav{32} = load('Data\190228\190228_Behavior.txt');
Behav{33} = load('Data\190301\190301_Behavior.txt');
Behav{34} = load('Data\190302\190302_Behavior.txt');
Behav{35} = load('Data\190303\190303_Behavior.txt');
Behav{36} = load('Data\190307\190307_Behavior.txt');
Behav{37} = load('Data\190308\190308_Behavior.txt');

Behav{38} = load('Data\201219\201219_Behavior.txt');
Behav{39} = load('Data\201220\201220_Behavior.txt');
Behav{40} = load('Data\201221\201221_Behavior.txt');
Behav{41} = load('Data\201222\201222_Behavior.txt');
Behav{42} = load('Data\201223\201223_Behavior.txt');
Behav{43} = load('Data\201226\201226_Behavior.txt');
Behav{44} = load('Data\201227\201227_Behavior.txt');
Behav{45} = load('Data\201228\201228_Behavior.txt');
Behav{46} = load('Data\201229\201229_Behavior.txt');
Behav{47} = load('Data\201230\201230_Behavior.txt');
Behav{48} = load('Data\200102\200102_Behavior.txt');
Behav{49} = load('Data\200103\200103_Behavior.txt');
% Behav{50} = load('Data\200104\200104_Behavior.txt');
Behav{50} = load('Data\200105\200105_Behavior.txt');
