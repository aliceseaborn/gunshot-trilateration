% ECE 301 - PROJECT 3
%   TRILATERATION DATA ANALYSIS
%
%   AUSTIN DIAL
%   THOMAS LANE
%
%   12/03/2018
%

clear; clc;

%% LOAD DATA SETS
%
    
    % Sampling rate of 125 kHz
    fs = 125000;
    
% Import Excel data
    
    %   CH1 TIME | CH1 DATA | CH2 TIME | CH2 DATA | CH3 TIME | CH3 DATA
    DATA = xlsread('C:\Users\adial\Documents\SCHOOL\8 FA18\ECE 301\Projects\Project 3\Analysis\DATA.xlsx');
    
% Adjust the time axis
    
    t = 1/125000:1/125000:1;
    t = t';
    
% Separate the channels and time v data
    
    CH1_time = t;
    CH1_data = DATA(:, 2);
    
    CH2_time = t;
    CH2_data = DATA(:, 4);
    
    CH3_time = t;
    CH3_data = DATA(:, 6);
    
    
%% PLOT ORIGINAL DATA
%

% Plot sampled signal

	figure(1);
	
	subplot(3,1,1);
	plot(CH1_time, CH1_data, 'm');
    grid on;
	xlabel('Time [s]');
	ylabel('Amplitude');
	title('CH1 Data');
	
	subplot(3,1,2);
	plot(CH2_time, CH2_data, 'c');
    grid on;
	xlabel('Time [s]');
	ylabel('Amplitude');
	title('CH2 Data');
    
	subplot(3,1,3);
    plot(CH3_time, CH3_data, 'b');
    grid on;
	xlabel('Time [s]');
	ylabel('Amplitude');
	title('CH3 Data');
    
    
%% ILLUSTRATE PEAKS
%

    figure(2);
    plot(CH3_time(22000:30000), CH3_data(22000:30000), 'm');
    grid on;
	xlabel('Time [s]');
	ylabel('Amplitude');
	title('CH3 Peak Close Up');
    

%% FIND A PEAK AND ILLUSTRATE IT
%

% Find peaks within specified range
    
    [pks,locs] = findpeaks(CH1_data, CH1_time, 'MinPeakDistance', 0.1, 'MinPeakHeight', 3.8, 'SortStr', 'descend');

% Plot data
    
    figure(3);
    plot(CH1_time, CH1_data, 'b');
    grid on;
	xlabel('Time [s]');
	ylabel('Amplitude');
	title('CH1 Data');
    hold on;

% Find peak position in data
    
    indx = find(abs(CH1_data) == pks(1));
    indx = indx(1);

% Add marker for first peak
    
	plot(CH1_time(indx), CH1_data(indx),'r*');
    text(CH1_time(indx), CH1_data(indx), '\leftarrow First peak');


%% FIND AND PLOT ALL PEAKS
%

% Find peaks within specified range
    
    [CH1pks, CH1locs] = findpeaks(CH1_data, CH1_time, 'MinPeakDistance', 0.1, 'MinPeakHeight', 3.8, 'SortStr', 'descend');
    [CH2pks, CH2locs] = findpeaks(CH2_data, CH2_time, 'MinPeakDistance', 0.1, 'MinPeakHeight', 3.8, 'SortStr', 'descend');
    [CH3pks, CH3locs] = findpeaks(CH3_data, CH3_time, 'MinPeakDistance', 0.1, 'MinPeakHeight', 3.8, 'SortStr', 'descend');

% Find peak position in data

    CH1indx = find(abs(CH1_data) == CH1pks(1));
    CH1indx = CH1indx(1);
    
    CH2indx = find(abs(CH2_data) == CH2pks(1));
    CH2indx = CH2indx(1);

    CH3indx = find(abs(CH3_data) == CH3pks(2));
    CH3indx = CH3indx(1);
    
%% PLOT ALL PEAK RESULTS
%
    
    figure(4);
    
    subplot(3,1,1);
	plot(CH1_time, CH1_data, 'm');
    grid on;
	xlabel('Time [s]');
	ylabel('Amplitude');
	title('CH1 Data');
    hold on;
    
    plot(CH1_time(indx), CH1_data(indx),'r*');
    text(CH1_time(indx), CH1_data(indx), '\leftarrow First peak');
	
	subplot(3,1,2);
	plot(CH2_time, CH2_data, 'c');
    grid on;
	xlabel('Time [s]');
	ylabel('Amplitude');
	title('CH2 Data');
    hold on;
    
    plot(CH2_time(indx), CH2_data(indx),'r*');
    text(CH2_time(indx), CH2_data(indx), '\leftarrow First peak');
    
	subplot(3,1,3);
    plot(CH3_time, CH3_data, 'b');
    grid on;
	xlabel('Time [s]');
	ylabel('Amplitude');
	title('CH1 Data');
    hold on;
    
    plot(CH3_time(indx), CH3_data(indx),'r*');
    text(CH3_time(indx), CH3_data(indx), '\leftarrow First peak');
    
% Purge unnecessary data

    clear CH1locs; clear CH2locs; clear CH3locs;
    clear CH1pks; clear CH2pks; clear CH3pks;

    
%% USE PEAKS TO FIND TIME DELAYS
%
    
% Find time of sound impact w/r to sampling
    
    % We subtracted a set constant time from each sample to get the scale
    %   of the system down to where the time differences are meaningful. I
    %   settled on the number 0.1943 s because that is the approximate time
    %   that the first peaks begin. This clipping makes the time zero set
    %   to where the buzzer was likely fired, rather than when the data
    %   recording began.
    
    impact = floor([ CH1indx, CH2indx, CH3indx ] - (0.1943 * fs) ) / fs;
    
    % Calculate distances using the speed of sound
    dist = impact * 343;    % Get distance from time
    dist = dist * 100;      % Convert from M to CM
    
    
%% INITIALIZATION OF TRILATERATION
%
    
% Initialize points and variables 

    N = 3;
    P = zeros(N, 2);    % Initialize 
    P(1, :) = [1.0 0.5];
    P(2, :) = [4.0 39.0];
    P(3, :) = [45.0 1.2];
    
    P_actual = [16.5 20.0];	% Actual position of agent

    
% Load distance data

    rho = dist';        % Distance data
    Drho = zeros(N, 1); % Change in data
    

    P_est = [3 4];      % Random guess for agent position
    

%% ESTIMATE POSITIONS
%

    imax = 100000000;   % Prespecify number of iterations
    i = 1;              % Instantiate iterations
    tol = 0.1;          % Define minimum acceptable tolerance

    % Tolerance:     LOW X, Y             HIGH X, Y
    target = [ (1-tol)*P_actual(1), (1+tol)*P_actual(1), 
           (1-tol)*P_actual(2), (1+tol)*P_actual(2) ];

    while ( ~( (( target(1, 1) <=  P_est(1) ) && ( P_est(1) <=  target(1, 2) )) && (( target(2, 1) <=  P_est(2) ) && ( P_est(2) <=  target(2, 2) )) ) && i < imax )
    
        for m = 1:N
            rho_est(m) = sqrt( (P_est(1)-P(m,1))^2 + (P_est(2)-P(m,2))^2 );
            Drho(m) = rho(m) - rho_est(m);
            H(:,:) = [P_est(1)-P(:,1), P_est(2)-P(:,2)];
        end
    
        deltaP = inv(transpose(H)*H) * transpose(H) * Drho;
        P_est = P_est + transpose(deltaP);
        P_est_old(i,:) = P_est;                  % Previous values of P_est for plotting
    
        % Advance the iterations
        i = i + 1;
    
    end
    

%% PLOT RESULTS
%

    % Plot actual positions
    figure(1);
    plot(P_actual(1),P_actual(2),'rx','markersize',12,'linewidth',2);
    grid on, hold on;
    
    % Add plot of approximations
    plot(P_est_old(:,1),P_est_old(:,2),'.','linewidth',2);
    title('Position Estimate');
    xlabel('X position'),ylabel('Y position');
    
    
    
    
    
    
    
    
    