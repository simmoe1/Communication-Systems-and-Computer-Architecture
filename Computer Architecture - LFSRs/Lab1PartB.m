%LAB 1 - PART B
%EBRAHIM SIMMONS AND BILAL YUSUF 
%400200042 and 400185626

clear all;
clc;

%---------------------------------PART 1------------------------------------
%taking in orginal picture
A = imread('Pencil_Image_DM4.jpg');
figure(1);
clf

%printing orginal image for comparisopm
fprintf('Orginal Image Figure 1\n'); %write image to Matlab figure window #1 
image(A);  

%setting values for rows, col, depth corsponding to image i.e. getting size
[rows,cols,depth] = size(A);

% colon (:) selects all available indices in a dimension 
%put here cause sometimes program recognizes as 2d matrix instead of 3
R = A(:,:,1); %copy rows and columns for depth = 1, to R matrix 
G = A(:,:,2); %copy rows and columns for depth = 2, to G matrix 
B = A(:,:,3); %copy rows and columns for depth = 3, to B matrix  

%Initialize matrices for A encrypt
A_encrypted = zeros(rows,cols,depth);
%Initialize matrices for A decrypt
A_decrypted = zeros(rows,cols,depth);
%Initialize matrices for Random matrix encrypt
RAND_matrix = zeros(rows,cols,depth);

%---------------------------------PART 2------------------------------------
%TRYING TO LOAD LSFT VALUES BUT GIVING US ERROR SO WE USED RANDOM VALUES
%INSTEAD LIKE FROM TUTORIAL/LECTURE
% LSFRNo = load("my_random_numbers.m");
% LSFRNo = LSFRNo(:);
% numsIndex = 1;

% we parse through row, columns, and depth
for rows = 1:rows
    for cols = 1:cols
        for depth = 1:depth
            %assigning random matrix values
            RAND_matrix(rows,cols,depth) = round (255*rand);
           % RAND_matrix(rows,cols,depth) = LSFRNo(numsIndex);

        end
    end
end

%---------------------------------PART 3-5------------------------------------
% we parse through row, columns, and depth
for rows = 1:rows
    for cols = 1:cols
        for depth = 1:depth
            
            %First we convert the random number to bit vector as stated on lab
            bitVectorRand = dec2bin(RAND_matrix(rows,cols,depth),8)-48;

            %2nd we convert the random number to bit vector as stated on lab
            bitVectorImg = dec2bin(A(rows,cols,depth),8)-48;
    
            %Xor the two vector and store as decimal 
            A_encrypted(rows,cols,depth) = binaryVectorToDecimal(bitxor(bitVectorRand,bitVectorImg));
            %Covert the encrypyted image as bit vector doing dec2bin 
            bitVectorEncrypt = dec2bin(A_encrypted(rows,cols,depth),8)-48;
            
            %Xor the two vector and store as decimal 
            A_decrypted(rows,cols,depth) = binaryVectorToDecimal(bitxor(bitVectorRand,bitVectorEncrypt));
            
            %value put here to ensure its running and program didnt hang
            fprintf('Running\n'); 
        end
    end
end

%----------------------------DISPLAY----------------------------------------
%Display figures 2 and 3 below 
figure(2);
clf;
fprintf('encrypted Image Figure 2\n')
image (uint8(A_encrypted));
figure(3);
clf;
fprintf('decrypted Image Figure 3\n')
image (uint8(A_decrypted))