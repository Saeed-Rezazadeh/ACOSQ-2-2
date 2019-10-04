%% The script corresponds to the Algorithm 4 with r = (2 2)
clc
clear all
close all
%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon.
% Also, since the ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for a given epsilon and delta.
FileID = fopen ('Results.txt' , 'a') ;

%% Source distributions
sigma = 1 ; % the standard devision of the source
u = 0 ; % the source's mean value
alpha = 300000  ; % the size of the training set to find the initial codebook using the splitting algorithm

%% Channel's cross-over probability epsilon
epsilon = unique ([10 ^ -6 10^-5 : 2 * 10^-5 : 10^-4 , 10^-4 :  10^-4  : 10^-3 , 10^-3 , 0.005 0.01  0.05  0.1]);
% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use a increase-decrease method.

SIZE = length(epsilon) ;
noise =  [1 : SIZE , SIZE : -1 : 1  , 1 : SIZE , SIZE : -1 : 1] ;

% The variable resolution determines the accuracy of the Riemann summation
resolution = 2 ^ 11 ;


%% Initialize parameters
SDR_2 = zeros (SIZE , 1) ;
SDR_4 = zeros (SIZE , 1) ;
Final_SDR_4 = zeros(SIZE , 1) ;


Probability_y_1_y_2 = zeros (4 , 1) ;
[Training_set , T , delta_u] = initialization (sigma , u , alpha , resolution) ;

% Compute the source pdf. We herein consider a zero-mean
% unit-variance Gaussian source distribution.
u = T(: , 1) ;
f =  1 ./ (sqrt (2 .* pi)) .* exp (-u .^ 2 ./ 2) ;
f = f./ (sum(f) * delta_u ) ;

%% Noise correlation
for delta = [ 0 5 10]
    D_4 = zeros (4 , 1) ;
    for k = 1 : length(noise)
        i = noise(k) ;
        Pr = Channel_with_Memory(4 , epsilon(i) , delta) ;
        
        %% COSQ for rate 2
        
        if (k == 1)
            % Using the splitting algorithm to find the initial codebook.
            [~  , codebook] = kmeans (Training_set , 4 , 'MaxIter',1000 , 'OnlinePhase','on') ;
        else
            % We slightly increase the channel's cross-over probability,
            % setting the codebook from the system with small epsilon as
            % the initial state of the system with new epsilon.
            load ('codebook_rate_2.mat' , 'codebook')
        end
        
        % The first step of the ACOSQ described in Section 4.1 and
        % Algorithm 4. In this step a 2-bit COSQ is designed for the source pdf f as computed in line 41.
        [SDR_2(k) ,  D_2 , T , codebook] = COSQ_2(Pr , f , T(: , 1) ,  codebook , delta_u) ;
        
        % save the codebook to initialize the system with the next value of
        % epsilon. In other words, the codebook obtained at every step is used to initialize
        % the quantizer in the SAME step for the new epsilon.
        save ('codebook_rate_2' , 'codebook' ) ;
        fprintf (FileID , '\nrate = 2\n') ;
        %% COSQ for rate 4
        codebook_4 = [] ;
        for y_1 = 1 : 2
            for y_2 = 1 : 2
                y = (y_1 - 1) * 2 + y_2 ;
                % Based on (4.8) compute the source conditional pdf given the received
                % sequence y_1y_2 where y_1y_2 is the channel output corresponding to
                % the transmitted sequence in the first step.
                [f_u_given_y_1_y_2] = generate_pdf_rate_2(y_1 , y_2 , T , Pr , f , delta_u ) ;
                
                for inner_noise = 1 : k
                    i = noise(inner_noise) ;
                    Pr_z = [(1 - epsilon(i) + delta) / (1 + delta)  , epsilon(i) / (1 + delta) ;
                        (1 - epsilon(i)) / (1 + delta)  , (epsilon(i) + delta) / (1 + delta)] ;
                    
                    
                    if (inner_noise == 1 )
                        % Find the initial codebook for the conditional source
                        % pdf using splitting algorithm.
                        codebook = init_codebook(f_u_given_y_1_y_2 , delta_u , T , alpha) ;
                    else
                        % We slightly increase the channel's cross-over probability,
                        % setting the codebook from the system with small epsilon as
                        % the initial state of the system with new epsilon.
                        Data = ['codebook_y_' num2str(y)] ;
                        load (Data) ;
                    end
                    % The second step of the ACOSQ described in Section 4.1 and
                    % Algorithm 4. In this step a 2-bit COSQ is designed for the conditional source pdf f_u_given_y_1_y_2
                    % as computed in line 79.
                    [SDR_4(k) , D_4(y) , hold_T , codebook] = ...
                        COSQ_4(f_u_given_y_1_y_2 , y_1 , y_2 , Pr_z , T(: , 1 : 2) , codebook , delta_u ) ;
                end
                T(: , 2 + y) = hold_T(: , 3) ;
                % Compute the probability of P(Y_1Y_2 = y_1y_2) for y_1y_2 = 00 , 01 , 10 , 11.
                Probability_y_1_y_2(y) = Pr_y_1_y_2 (y_1 , y_2 , f , T(: , [1 2]) , delta_u , Pr) ;
                Data = ['codebook_y_' num2str(y)] ;
                save (Data , 'codebook' ) ;
                codebook_4  = [codebook_4 ; codebook] ;
            end
        end
        % Store the ultimate partition indexes and codebook to compute
        % experimental results.
        Data = ['T\T_k_' num2str(k) '_delta_' num2str(delta)] ;
        save(Data , 'T' , 'codebook_4') ;
        
        Probability_y_1_y_2 = Probability_y_1_y_2 ./ sum(Probability_y_1_y_2) ;
        % As mentioned in the Thesis, the ultimate distortion at every step is the
        % weighted sum of the conditional distortion functions given the
        % received sequence y_1y_2.
        Final_D_4 = sum(D_4 .* Probability_y_1_y_2) ;
        % Compute the SDR value in the following way as the source is zero
        % mean and unit variance.
        Final_SDR_4(k) = 10 * log10(1 / Final_D_4) ;
        fprintf (FileID , 'Final SDR = %4.2f\n' ,  Final_SDR_4(k)) ;
        fprintf (FileID , '\nrate = 4\n') ;
        fprintf (FileID , 'noise = %d\n' ,  i) ;
    end
    % Pick the best SDR value after the end of the so-called
    % increase-decrease method.
    fprintf (FileID , '\nSDR for rate 2\n') ;
    clear D_2
    D_2 = zeros (SIZE , 1) ;
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var = SDR_2(index) ;
        hold_var = hold_var(:) ;
        final_SDR_2(i) = max(hold_var) ;
        fprintf (FileID , '\ni = %d' , i) ;
        D_2(i) = 10 ^(- final_SDR_2(i) / 10) ;
        fprintf (FileID , '\nD = %f' , D_2(i)) ;
        fprintf (FileID , '\nSDR_2 = %7.4f' , final_SDR_2(i)) ;
    end
    
    fprintf (FileID , '\nSDR for rate 4\n') ;
    clear D_4 ;
    D_4 = zeros(SIZE , 1) ;
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var = Final_SDR_4(index) ;
        hold_var = hold_var(:) ;
        final_SDR_4(i) = max(hold_var) ;
        fprintf (FileID , '\ni = %d' , i) ;
        D_4(i) = 10 ^(- final_SDR_4(i) / 10) ;
        fprintf (FileID , '\nD = %f' , D_4(i)) ;
        fprintf (FileID , '\nSDR_4 = %7.4f\n' , final_SDR_4(i)) ;
    end
    %  Save the best SDR values for every channel parameters.
    clear D_4 D_2 ;
    Data = ['ACOSQ_2_2_delta_' num2str(delta)] ;
    save (Data , 'final_SDR_4' , 'Final_SDR_4' , 'final_SDR_2' , 'SDR_2' , 'epsilon') ;
end
