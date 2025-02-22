Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

# MATLAB code for analyzing microelectrode, scalp EEG and behavioural data

## Description

This repository contains MATLAB analyses of single unit and multiunit electrophysiology recordings from the human subthalamic nucleus. Analyses of scalp electroencephalogram recordings. Analyses of local field potential recordings from the subthalamic nucleus.

## Content

- .m files for data analysis
- license file

## Installation

- Download the .m file of this package and dependent packages (see dependencies below).
- Add directories (with subdirectories) with the downloaded .m files to the Matlab path.
- Run mountcb.m. Name the mounted Cellbase 'PD_CellBase'. When prompted to 'locate CellBase database file', locate and choose the CellBase.mat file from the data. Choose 'CellBase' for the time stamp conversion question.
- Main function: PD_SSRT_MAIN.m
- Typical installation time: 15 minutes

## Instructions for use

- Instructions to run on data: call PD_SSRT_MAIN.m from Matlab command window.
- Expected output: data analysis figures regarding behavioral responsiveness and local synchrony (see https://www.medrxiv.org/content/10.1101/2024.12.02.24318298v1)
- How to run software on your own data: mount your own CellBase in 'PD_CellBase'. Please observe standard CellBase formatting requirements (see https://github.com/hangyabalazs/CellBase).

## Demo

- Create a directory for the saved files.
- Download the example data set.
- We refer to the full path where the data is located on your computer as filesdir.
- We refer to the full path where the results will be saved on your computer as rootdir.
- After installation, call PD_SSRT_MAIN(rootdir,filesdir).
- Typical run time on demo data: 200 minutes.

## Dependencies

- CellBase data base system. MATLAB code for CellBase is available at: www.github.com/hangyabalazs/CellBase.
- MATLAB analysis package: https://github.com/hangyabalazs/Hangya-Matlab-Code.
- EEGLAB toolbox (Delorme A & Makeig S, 2004, 10.1016/j.jneumeth.2003.10.009): https://github.com/sccn/eeglab
- The following plugins of EEGLAB:
        -FileIO: https://github.com/fieldtrip/website/blob/master/development/module/fileio.md
	-ICLabel: https://github.com/sccn/ICLabel
	-firfilt: https://github.com/sccn/firfilt
	-CleanLine(Tim Mullen, 2011): https://github.com/sccn/cleanline
- FastICA toolbox: https://research.ics.aalto.fi/ica/fastica/
- FieldTrip toolbox (Oostenveld, R., 2011, doi:10.1155/2011/156869): https://github.com/fieldtrip/fieldtrip 
- DBSFILT toolbox (Guillaume Lio,2012, https://github.com/guillaumelio/DBSFILT)
- CSD toolbox (Jürgen Kayser, 2009, doi:10.1016/j.clinph.2005.08.034)
- Confidence intervals for regression are derived using the polypredci.m function (Star Strider, https://www.mathworks.com/matlabcentral/fileexchange/57630-polypredci, MATLAB Central File Exchange, retrieved December 30, 2020).

## System requirements  

Windows 10 64bits  
Any Intel or AMD x86-64 processor  
8 GB RAM  
MatlabR2019b or higher  
MATLAB Statistics Toolbox  
MATLAB Image Processing Toolbox
MATLAB Signal Processing Toolbox

Please contact us with any questions, bug reports, and general feedback:

Johanna Petra Szabó and Balázs Hangya  
hangya.balazs@koki.hu
