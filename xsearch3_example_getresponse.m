function output = xsearch3()

% generic visual search architecture
% started 12/06/2005
% current 12/06/2005
% Todd S. Horowitz

global screenRect MainWindow

% this bit should make the code windows/macos9 portable
if strcmp(computer, 'PCWIN')
        warning off MATLAB:DeprecatedLogicalAPI;
end

% get parameter input (optional)
prompt = {'subject''s initials', 'practice trials',  'experimental trials', 'set sizes'};
defaults = {'xxx',  '5', '50', '[4 8 16]'};
answer = inputdlg(prompt, 'Experimental Setup Information', 1, defaults);

% now decode answer
[subject, pTrials, xTrials, setsizes] = deal(answer{:});
pTrials = str2num(pTrials);
xTrials = str2num(xTrials);
setsizes = str2num(setsizes);

% other major variables
formatString = '%s, %d, %d, %d, %d, %d, %d, %d, %d \n'; % this is for data output
% saved variables: subject, block, trial, target, setSize, timeoutFlag, wrongKeyFlag, err,  RT
keys = [KbName('''')  KbName('a')]; % response key assignments
dataFileName = 'xsearch3Data';
jitter = 10; % coordinates of stimuli will vary by up to this much in their grid cells

if exist(dataFileName, 'file') == 0
        dataFile = fopen(dataFileName, 'a');
        count = fprintf(dataFile, '%s \n', 'subject, block, trial, target, setSize, timeoutFlag, wrongKeyFlag, err,  RT');
        fclose('all');
end

% set the state of the random number generator to some random value (based on the clock)
tau = clock;
rand('state',sum(100*tau));
state=rand('state');
seedFileName = ['pptxSeed-', subject];
save (seedFileName, 'tau');

% basic setup stuff
screens = Screen('Screens');
MainScreen = max(screens);
% video setup

% this sets up a display which is the size of the current monitor
% the color is 8 bits, the lowest color resolution, this is for speed and compactness
% if you're using images, you might want to change bitDepth to 32
% on MacOS9, you can also specify a screenRect smaller than the actual monitor

bitDepth = 8;
MainWindow = screen(0, 'OpenWindow', [], [], bitDepth); 
hz = round(screen(MainWindow, 'FrameRate'))
HideCursor;

% define rects
screenRect = screen(0, 'Rect');
objectSize=40;
objectRect = [0 0 objectSize objectSize]; % this is the size of the individual items
stimulusRect = [0 0 500 500]; % this is an area smaller than the full screen in which to present stimuli
fieldRect = CenterRect(stimulusRect, screenRect);
% define coordinate system

cellX = 100; cellY = cellX; % these are the cells of the imaginary stimulus matrix
nColumns = 5;
nRows = 5;
columns = ((1:nColumns) - .5) * cellX;
rows = ((1:nRows) - .5) * cellY;
[X, Y] = meshgrid(columns, rows); % this is the coordinate system for the stimuli, relative to stimulusRect
nPositions = nColumns * nRows;
nSetSizes = length(setsizes);
maxWait = 5; % allow 5 seconds for response
fixationDuration = 30; % this is in refreshes

% now set colors
dacBits = ScreenDacBits(MainWindow);
white = WhiteIndex(MainWindow);
black = BlackIndex(MainWindow);
midGray = round((white+black)/2);
if round(midGray)==white
        midGray=black;
end
darkGray = (midGray+black)/2;
lightGray = (midGray+white)/2;
red = 1;
blue = 2;
green = 3;
clut = screen(MainWindow, 'GetClut');
clut (red + 1, :) = [255 0 0];
clut (blue + 1, :) = [0 0 255];
clut (green + 1, :) = [0 255 0];
if dacBits == 10
        clut = clut*1024/256;
end
LoadClut(MainWindow, clut);

screen(MainWindow, 'FillRect', midGray);
screen(MainWindow, 'TextFont', 'Skia');
screen(MainWindow, 'TextSize', 24);

% more offscreen windows
backgroundColor = midGray;
stimulus = screen(MainWindow, 'OpenOffscreenWindow', midGray, stimulusRect);

% make stimuli
% here I'm drawing the stimuli, but one can also use alphanumeric characters or import jpegs
itemColor = black;
stroke = 10;
topRect = [0 0 objectSize stroke];
bottomRect = [0 objectSize - stroke objectSize, objectSize];
hCrossRect = [0 (objectSize - stroke)/2 objectSize (objectSize + stroke)/2];
vCrossRect = [(objectSize - stroke)/2 0 (objectSize + stroke)/2 objectSize];
leftRect = [0 0 stroke, objectSize];
rightRect = [objectSize - stroke 0 objectSize objectSize];

t(1) = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(t(1), 'FillRect', itemColor, topRect);
screen(t(1), 'FillRect', itemColor, vCrossRect);
t(2) = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(t(2), 'FillRect', itemColor, rightRect);
screen(t(2), 'FillRect', itemColor, hCrossRect);
t(3) = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(t(3), 'FillRect', itemColor, bottomRect);
screen(t(3), 'FillRect', itemColor, vCrossRect);
t(4) = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(t(4), 'FillRect', itemColor, leftRect);
screen(t(4), 'FillRect', itemColor, hCrossRect);

d(1) = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(d(1), 'FillRect', itemColor, leftRect);
screen(d(1), 'FillRect', itemColor, bottomRect);
d(2) = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(d(2), 'FillRect', itemColor, leftRect);
screen(d(2), 'FillRect', itemColor, topRect);
d(3) = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(d(3), 'FillRect', itemColor, rightRect);
screen(d(3), 'FillRect', itemColor, topRect);
d(4) = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(d(4), 'FillRect', itemColor, rightRect);
screen(d(4), 'FillRect', itemColor, bottomRect);

fixation = screen(MainWindow, 'OpenOffscreenWindow', backgroundColor, objectRect);
screen(fixation, 'FillRect', itemColor, vCrossRect);
screen(fixation, 'FillRect', itemColor, hCrossRect);
fixationRect = CenterRect(objectRect,screenRect);

screen(MainWindow, 'FillRect', backgroundColor); % this blanks the screen
instructionString = {
'Your task is to look for a T in the display';
'If you see a T, press the "quote" key with your right hand.'
'If you are sure that there is no T in the display, press the "a" key with your left hand.';
'Please respond as quickly and accurately as possible.';
['There will be ', num2str(pTrials), ' practice trials followed by ', num2str(xTrials), ' experimental trials.'];
};
CenterCellText(MainWindow, instructionString);
SitnWait;
% block routine

for block = 1:2
        if block == 1
                blockMessage = 'Practice ';
                nTrials = pTrials;
        else
                blockMessage = 'Experimental ';
                nTrials = xTrials;
        end

        nTrialString = num2str(nTrials);
        screen(MainWindow, 'FillRect', midGray);
        message = [' Press any key to begin ', nTrialString, ' ', blockMessage, 'trials'];
        CenterText(message);
        SitnWait;
        
        % trial routine
        for trial = 1:nTrials
                screen(MainWindow, 'FillRect', backgroundColor); %start with a blank screen
                screen(stimulus, 'FillRect', backgroundColor); %start with a blank screen
                
                % now select trial properties
                target = randi(2); % 1 = target present, 2 = target absent
                setSize = setsizes(randi(nSetSizes));
                positionIndex = randperm(nPositions); % randomize order of stimulus placement
                selectItem = mod(randperm(setSize ) - 1, 4) + 1; % this ensures that distractor types will be distributed evenly
                
                % put stimuli onto stimulus screen
                for i = 1:setSize
                        if target == 1 & i == 1
                                item = t(selectItem(1));
                        else
                                item = d(selectItem(i));
                        end
                        dx = round(randi(jitter) - jitter/2);
                        dy = round(randi(jitter) - jitter/2);
                        itemRect = CenterRectOnPoint(objectRect, X(positionIndex(i)) + dx, Y(positionIndex(i)) + dy);
                        screen('CopyWindow', item, stimulus, objectRect, itemRect);
                end
        
                % present fixation cross
                screen('CopyWindow', fixation, MainWindow, objectRect, fixationRect);
                screen(MainWindow, 'WaitBlanking', fixationDuration); % synch to the vertical refresh for fixationDuration refreshes
                % now copy display to screen in one fell swoop
                screen('CopyWindow', stimulus, MainWindow, stimulusRect, fieldRect);
                
                % now wait for a response
                initTime = GetSecs;
                response = [];
                
                while isempty(response)&(GetSecs - initTime) < maxWait
                        [keyIsDown, KbTime, keyCode] = KbCheck;
                        if keyIsDown
                                response = find(keyCode);
                                response = response(1);
                                responseTime = KbTime;
                        end
                end
        
                if isempty(response)
                        timeoutFlag = 1;
                        err = 1;
                        feedback = 'Time out!';
                        RT = 0;
                else
                        timeoutFlag = 0;
                        responseKey = find(keys == response);
                        wrongKeyFlag = 0;
                        RT = round((responseTime - initTime) * 1000); % RT in ms
                        if isempty(responseKey)
                                err = 1;
                                feedback = 'Wrong Key! Use "a" for "no" and "quote" for "yes"!';
                                wrongKeyFlag = 1;
                        elseif responseKey == target
                                err = 0;
                                feedback = 'Correct!';
                        else
                                err = 1;
                                feedback = 'Wrong!';
                        end
                end
        
                % save trial data
                % note that you should open and close the file every time you write to it, otherwise you may lose data in crashes
                dataFile = fopen(dataFileName, 'a');
                fprintf(dataFile, formatString, subject, block, trial, target, setSize, timeoutFlag, wrongKeyFlag, err,  RT);
                fclose(dataFile);
                
        
                screen(MainWindow, 'FillRect', midGray);
                CenterText([feedback, ' - RT = ', num2str(RT), ' ms']);
                WaitSecs(1);
                
        end % trial loop
end % block loop

% clean up and go home

screen(MainWindow, 'FillRect', midGray);
CenterText('Thank you for participating');
SitnWait;
clear screen;

function SitnWait()

FlushEvents('keyDown');
GetChar;

function CenterText (message, xoffset, yoffset, color, window)
% print a text string centered on the screen
% syntax [newX, newY] = CenterText (message, [xoffset], [yoffset], [color], [window])
% if you want the text offset from center, use xoffset and yoffset
% if window is not specified, prints to MainWindow, which must be a global in the calling function
% 2/23/2000 accepts color option
% 4/23/2002 can now print to offscreen windows

global screenX
global screenY
global MainWindow

switch nargin
case 1
        xoffset=0;
        yoffset=0;
        color = [];
        window = MainWindow;
case 2
        yoffset=0;
        color = [];
        window = MainWindow;
case 3
        color = [];
        window = MainWindow;
case 4
        window = MainWindow;
end

windowRect = screen(window, 'Rect');
width = screen(window, 'TextWidth', message);
screen(window, 'DrawText', message, ((windowRect(3)/2)-(width/2))+xoffset, (windowRect(4)/2)+yoffset, color);

function CenterCellText(window, messages, spacing, color)

% prints a cell array of text strings centered on the screen
% syntax [newX, newY] = CenterCellText(window, messages, [spacing], [color])
% messages is a cell array of lines of text
% spacing governs the distance between lines of text, in pixels
% default spacing is 20 pixels
% added color option 1/17/2003; only does one color for the whole thing
% fixed number of options bug 01/24/2003

defaultSpacing = 30;
switch nargin
case 0
        error('CenterCellText requires at least two arguments!');
case 1
        error('CenterCellText requires at least two arguments!');
case 2
        spacing = defaultSpacing;
        color = [];
case 3
        color = [];
end

lines = length (messages);

% find the middle line
middleLine = round(lines/2);

yOffset = spacing*((1:lines)-middleLine);
for y = 1:lines
        CenterText(messages{y}, 0, yOffset(y), color, window);
end