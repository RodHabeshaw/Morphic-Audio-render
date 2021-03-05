%%% Quantive Sample Morph - Rod Habeshaw - March 2021
%%% Not well commented - apologies
%%% Will 'morph' one audio loop into another using linear time and fourier
%%% shift, then performans a dual discrete component swap in time and fourier which needs work.
%%% Files are rendered to the set folder.

clear all
%% Input Reference Signal analysis

%[filename, pathname]=uigetfile('.mp3','Choose a start file. Smaller = Faster');
[filename1, pathname]=uigetfile('.wav','Choose a START file.');

filepath1=[pathname filename1]
[y1,Fs] = audioread(filepath1);

[filename2, pathname]=uigetfile('.wav','Choose an END file.');
filepath2=[pathname filename2]
[y2,Fs] = audioread(filepath2);

Y1=fft(y1);
Y2=fft(y2);
%%
steps=32; %% No of loops

%YY=Y1;Y1=Y2;Y2=YY;

t=(0:(steps-1))/(steps-1);
g=cos(2*pi*t)+j*sin(2*pi*t);

if size(Y1,2)~=size(Y2,2); Y1=Y1(:,1);Y2=Y2(:,1);end  

if size(Y1,1)<size(Y2,1); Y1=[Y1' zeros(size(Y2,1)-size(Y1,1),size(Y1,2))']';
else Y2=[Y2' zeros(size(Y1,1)-size(Y2,1),size(Y1,2))']';
end

%% go fourier

D=Y2-Y1;
OUT=zeros(size(Y1,2));
for n=1:steps
    A=real(ifft(Y1+(n-1)*(D/(steps-1))*g(n)));

    OUT=[OUT' A']';
    
end

% USE WAV  --  
audiowrite(['E:\--Audio Store--\Matlab_Shenans\Morpher\' filename1 '_' filename2 '_' num2str(steps) 'steps_merge.wav'],OUT,Fs);


%% Then redo for sample shift in Time Domain only

if size(y1,2)~=size(y2,2); y1=y1(:,1);y2=y2(:,1);end  

if size(y1,1)<size(y2,1); y1=[y1' zeros(size(y2,1)-size(y1,1),size(y1,2))']';
else y2=[y2' zeros(size(y1,1)-size(y2,1),size(y1,2))']';
end

d=y2-y1;
out=zeros(size(y1,2));
for n=1:steps
    a=y1+(n-1)*(d/(steps-1));

    out=[out' a']';
    
end

audiowrite(['E:\--Audio Store--\Matlab_Shenans\Morpher\' filename1 '_' filename2 '_' num2str(steps) 'steps_time.wav'],out,Fs);




%% Run sample/fft shift loops -- needs work

%% Swap to mono 
y1=y1(:,1);y2=y2(:,1);
Y1=Y1(:,1);Y2=Y2(:,1);


start_a=[[1:size(y1,1)]' y1 ];
end_b=[[1:size(y2,1)]' y2];
start_A=[[1:size(Y1,1)]' Y1 angle(Y1)]; %% for stack by fourier angle
end_B=[[1:size(Y2,1)]' Y2 angle(Y2)];


a_t=sortrows(start_a,2);
b_t=sortrows(end_b,2);
A_t=sortrows(start_A,3);
A_t(:,3)=[];
B_t=sortrows(end_B,3);
B_t(:,3)=[];

d_t=b_t-a_t;
D_t=B_t-A_t;

gain_line=hanning(steps);

gain_line=gain_line/sum(gain_line);


out=zeros(1)';OUT=zeros(1)';

%% differences and positions determined - next loops reconstitute audio


for n=1:steps
    n
    loop(1:size(a_t,1))=0;
    LOOP(1:size(A_t,1))=0;
     
   
   gain=sum(gain_line(1:n));
   
    mult=d_t*((n-1)/(steps-1));%*gain;%*(d_t/(steps-1));
    MULT=D_t*((n-1)/(steps-1));%*gain;%*(D_t/(steps-1));
    
    for m=1:size(a_t,1)
        indx=ceil(a_t(m,1)+mult(m,1));
        INDX=ceil(A_t(m,1)+MULT(m,1));
        if indx==0;indx=1;end
        if INDX==0;INDX=1;end
loop(indx)=loop(indx)+a_t(m,2)+mult(m,2);
LOOP(INDX)=LOOP(INDX)+A_t(m,2)+MULT(m,2)*g(n);

    end
    
    if size(loop,2)>size(loop,1);loop=loop';end
    if size(LOOP,2)>size(LOOP,1);LOOP=LOOP';end

    out=[out' loop']';
    
    LOOP=real(ifft(LOOP));
    OUT=[OUT' LOOP']';
    
end

audiowrite(['E:\--Audio Store--\Matlab_Shenans\Morpher\' filename1 '_' filename2 '_' num2str(steps) 'morphin_time.wav'],out,Fs);
audiowrite(['E:\--Audio Store--\Matlab_Shenans\Morpher\' filename1 '_' filename2 '_' num2str(steps) 'morphin_fourier.wav'],out,Fs);


