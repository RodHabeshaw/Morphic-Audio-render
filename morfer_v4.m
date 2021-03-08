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
steps=16; %% No of loops

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


%start_a=[[1:size(y1,1)]' y1 ];
%end_b=[[1:size(y2,1)]' y2];
%start_A=[[1:size(Y1,1)]' Y1 angle(Y1)]; %% for stack by fourier angle
%end_B=[[1:size(Y2,1)]' Y2 angle(Y2)];

[PKS1 LOCS1 W1 P1 ] =findpeaks(smooth(abs(Y1)));
[PKS2 LOCS2 W2 P2 ] =findpeaks(smooth(abs(Y2)));

%sort by peak
A_fp=sortrows([PKS1 LOCS1 W1 P1 ],1);
B_fp=sortrows([PKS2 LOCS2 W2 P2 ],1);
%equal length
if size(A_fp,1)<size(B_fp,1);A_fp=[zeros(size(B_fp,1)-size(A_fp,1),4)' A_fp']';end
if size(B_fp,1)<size(A_fp,1);B_fp=[zeros(size(A_fp,1)-size(B_fp,1),4)' B_fp']';end
%find difference line
D_fp=B_fp-A_fp;
Dr_fp=A_fp-B_fp;

%reconstitute for each iteration
%A as source and B as source

g1=(0:(steps-1))/(steps-1);
g2=1-g1;  %gain lines to merge the source list for difference sets
%a_t=sortrows(start_a,2);
%b_t=sortrows(end_b,2);
%A_t=sortrows(start_A,3);
%A_t(:,3)=[];
%B_t=sortrows(end_B,3);
%B_t(:,3)=[];
%d_t=b_t-a_t;
%D_t=B_t-A_t;
%gain_line=hanning(steps);
%gain_line=gain_line/sum(gain_line);
out=zeros(1)';OUT=zeros(1)';
% differences and positions determined - next loops reconstitute audio


for n=1:steps
    n
    LOOP(1:size(Y1,1))=0;FT_L(1:size(Y1,1))=0;
    %LOOP(1:size(A_t,1))=0;
    %gain=sum(gain_line(1:n));
    %mult=d_t*((n-1)/(steps-1));%*gain;%*(d_t/(steps-1));
    M1=D_fp*((n-1)/(steps-1));%*gain;%*(D_t/(steps-1));
    M2=Dr_fp*(1-((n-1)/(steps-1)));
     
    
    
    
    for m=1:size(D_fp,1)
        
        %%determine bounds
        st_A_fp=ceil(A_fp(m,2)-A_fp(m,3));if st_A_fp<1;st_A_fp=1;end
        ed_A_fp=floor(A_fp(m,2)+A_fp(m,3));if ed_A_fp>size(Y1,1);st_A_fp=size(Y1,1);end
        if ed_A_fp<1;ed_A_fp=1;end
        st_AL_fp=ceil(A_fp(m,2)+M1(m,2)-(A_fp(m,3)+M1(m,3)));if st_AL_fp<1;st_AL_fp=1;end
        ed_AL_fp=floor(A_fp(m,2)+M1(m,2)+(A_fp(m,3)+M1(m,3)));if ed_AL_fp>size(Y1,1);st_AL_fp=size(Y1,1);end
        if ed_AL_fp<1;ed_AL_fp=1;end
        
        
        %resample in to fit out
        %  ed_A_fp-st_A_fp
        %  ed_AL_fp-st_AL_fp
        bit_A=Y1(st_A_fp:ed_A_fp)+M1(m,1);
        
        p=ed_AL_fp-st_AL_fp;if p<1;p=1; end
        q=ed_A_fp-st_A_fp;if q<1;q=1;end
        
        rbit_A=resample(bit_A,p,q);
        
        
        AAA=FT_L(st_AL_fp:st_AL_fp+length(rbit_A)-1);
        BBB=(rbit_A-1);
        
        if size(AAA,1)>size(BBB,1); FT_L(st_AL_fp:st_AL_fp+length(rbit_A)-1)=AAA+BBB';end
        if size(AAA,1)<size(BBB,1); FT_L(st_AL_fp:st_AL_fp+length(rbit_A)-1)=AAA'+BBB;end
        
%        FT_L(st_AL_fp:st_AL_fp+length(rbit_A)-1)=FT_L(st_AL_fp:st_AL_fp+length(rbit_A)-1)+(rbit_A-1);
        
        
        
        st_B_fp=ceil(B_fp(m,2)-B_fp(m,3));if st_B_fp<1;st_B_fp=1;end
        ed_B_fp=floor(B_fp(m,2)+B_fp(m,3));if ed_B_fp>size(Y2,1);st_B_fp=size(Y2,1);end
        if ed_B_fp<1;ed_B_fp=1;end
        st_BL_fp=ceil(B_fp(m,2)+M2(m,2)-(B_fp(m,3)+M2(m,3)));if st_BL_fp<1;st_BL_fp=1;end
        ed_BL_fp=floor(B_fp(m,2)+M2(m,2)+(B_fp(m,3)+M2(m,3)));if ed_BL_fp>size(Y2,1);st_BL_fp=size(Y2,1);end   
        if ed_BL_fp<1;ed_BL_fp=1;end
        
        bit_B=Y2(st_B_fp:ed_B_fp)+M2(m,1);
        
        p=ed_BL_fp-st_BL_fp;if p<1;p=1; end
        q=ed_B_fp-st_B_fp;if q<1;q=1;end   
        
        rbit_B=resample(bit_B,p,q);
        
        AAA=FT_L(st_BL_fp:st_BL_fp+length(rbit_B)-1);
        BBB=(rbit_B-1);
        
        if size(AAA,1)>size(BBB,1); FT_L(st_BL_fp:st_BL_fp+length(rbit_B)-1)=AAA+BBB';end
        if size(AAA,1)<size(BBB,1); FT_L(st_BL_fp:st_BL_fp+length(rbit_B)-1)=AAA'+BBB;end
        
        %FT_L(st_BL_fp:st_BL_fp+length(rbit_B)-1)=AAA+(rbit_B-1)';
        
        % resample blocks to fit
        
                
    %    indx=ceil(a_t(m,1)+mult(m,1));
   %     INDX=ceil(A_t(m,1)+MULT(m,1));
  %      if indx==0;indx=1;end
 %       if INDX==0;INDX=1;end
%loop(indx)=loop(indx)+a_t(m,2)+mult(m,2);
%LOOP(INDX)=LOOP(INDX)+A_t(m,2)+MULT(m,2)*g(n);  %forgot to unfft?
LOOP=FT_L;
    end
    
    %if size(loop,2)>size(loop,1);loop=loop';end
    if size(LOOP,2)>size(LOOP,1);LOOP=LOOP';end

    %out=[out' loop']';
    
    LOOP=real(ifft(LOOP));
    OUT=[OUT' LOOP']';
    
end

%audiowrite(['E:\--Audio Store--\Matlab_Shenans\Morpher\' filename1 '_' filename2 '_' num2str(steps) 'morphin_time.wav'],out,Fs);
audiowrite(['E:\--Audio Store--\Matlab_Shenans\Morpher\' filename1 '_' filename2 '_' num2str(steps) 'morphin_fourier.wav'],OUT,Fs);


