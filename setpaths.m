function [] = setpaths
% 
% This function is used to set MATLAB search path to current directory and
% all of its subdirectories in order to allow main functions to use 
% auxiliary procedures stored there.
% =========================================================================
%               PLEASE, RUN THIS FUNCTION FIRST 
%     EACH TIME YOU START A NEW SESSION WITH THIS PROGRAM.
% =========================================================================

addpath(genpath(['.' filesep]));

