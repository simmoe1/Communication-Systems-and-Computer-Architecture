%LAB 1 - PART A
%EBRAHIM SIMMONS AND BILAL YUSUF 
%400200042 and 400185626

close all;
clear;
clc;

%VARIABLE DECLARATION 
bit = 18; %# of bits in LFSR
period = (2^bit) - 1; %period 
numOfBytes = floor(period/8); %# of bytes 

S = zeros(1, bit); %creating a row of vectors with 22 zeros ranging from S0 to S21
S(1,1) = 1; %set the inital state to 1 
S_initial = S;
DATA_OUT = []; %long vector to hold randomly generated bits 

%---------------------------------PART 1------------------------------------
while true

    %holding values of first and last bit 
    bit1_initial = S(1, 1); % bit value being moved out
    bit22_initial = S(1, bit); % bit value being passed through XOR 
    
    S(1,1:bit-1) = S(1,2:bit); %shift all bits 1 to the left

    %doing XOR with bit 1 and bit 22
    S(1,bit) = bit1_initial; 
    S(1,bit-1) = xor(bit1_initial, bit22_initial);

    DATA_OUT(end+1) = bit1_initial; %append b1_initial to data_out

    if (S == S_initial)
        fprintf("S_intial is equal to S so break; \n");
        break;
    end
end

%---------------------------------PART 3------------------------------------
%convert to sequence of  bytes and write to external file
fid = fopen("my_random_numbers.txt", "w");

[Row,Col] = size(DATA_OUT); %set row col size
byte = []; %array for bytes

temp = 0; %counter
for axis = 1:Col
    byte(end+1) = DATA_OUT(axis);
    NumBytes = mod(axis, 8);

    if NumBytes == 0 %if multiple of 8/remainder is 0 
        if temp == 16 %start new line after 16 integers 
            fprintf(fid, '\n');
            temp = 0;
        end
        fprintf(fid,"%3g, ", bin2dec(sprintf('%d', byte)));
        %%final_final(end+1) = bin2dec(sprintf('%d', byte));
        temp = temp + 1;
        byte = [];
    end
end

%---------------------------------PART 4-6------------------------------------
%CONSECUTIVE 1 RUNS 
DATA_OUT_1 = (DATA_OUT);
DATA_OUT_1 = sprintf('%d',DATA_OUT_1); %formats into a string
t1 = textscan(DATA_OUT_1,'%s','delimiter','0','multipleDelimsAsOne',1); %scans DATA_OUT0 to look for values in between zeros to find one runs

p = t1{:}; %variable to hold previous values 
data = cellfun('length', p); %Finds length of zero runs 

%number_times1 number of times one appears 
%run_length1 is longest one run 
[numberTimes_1, runLength_1] = hist(data, [1:max(data)]);

conditional_probability = [];%empty array for condiontal probability of 1 run
for k = 1:runLength_1(end)
    conditional_probability(end+1) = numberTimes_1(k)/sum(numberTimes_1);%uses formula from lab document to find the conditonal probability 
end

%CONSECUTIVE 0 RUNS 
DATA_OUT_0 = (DATA_OUT); 
DATA_OUT_0 = sprintf('%d',DATA_OUT_0); %formats into a string
t1 = textscan(DATA_OUT_0,'%s','delimiter','1','multipleDelimsAsOne',1); %scans DATA_OUT0 to look for values in between zeros to find zeros runs

p = t1{:};
data = cellfun('length', p); %Finds length of zero runs 

%number_times0 number of times zero appears 
%run_length 0 is longest zero run 
[numberTimes_0, runLength_0] = hist(data, [1:max(data)]); 

%empty array for condiontal probability of zero run
conditional_probability_0 = [];
for k = 1:runLength_0(end)
    %uses formuila from lab document to find the conditonal probability 
    conditional_probability_0(end+1) = numberTimes_0(k)/sum(numberTimes_0);
end

%Unconditional probability -----------------------------------------
unconditional_probability_0 = [];
for k = 1:runLength_0(end)
     %uses formuila from lab document to find the unconditonal probability 
     unconditional_probability_0(end+1) = numberTimes_0(k)/(sum(numberTimes_0) + sum(numberTimes_1));
 end

 unconditional_probability_1 = [];
 for k = 1:runLength_1(end)
    %uses formuila from lab document to find the unconditonal probability 
     unconditional_probability_1(end+1) = numberTimes_1(k)/(sum(numberTimes_0) + sum(numberTimes_1));
 end
 unconditional_probability_0
 unconditional_probability_1

table_0runs = [runLength_0; numberTimes_0; conditional_probability_0;]; %table for 0 runs 
table_1runs = [runLength_1; numberTimes_1; conditional_probability]; %table for 1 runs 
%taking three vectors and putting them into one

%----------------------------DISPLAY----------------------------------------
disp('A1 SOLUTION:')
disp(DATA_OUT)
disp('A5 SOLUTION: consecutive 0 runs')
disp(table_0runs)
disp('A6 SOLUTION: consecutive 1 runs')
disp(table_1runs)