%Alleskönner Skript
% kann: -einzelne Spektren NOT WORKING YET
%       -Kinetics (shift correction)NOT WORKING YET
%       -Stark shift (in einzelnen Messungen, Normal)NOT WORKING YET
%       -Stark (in Kinetics) NOT WORKING YET
%       -PLE NOT WORKING YET
%       -Saturiation NOT WORKING YET
%       -Polarization NOT WORKING YET
%
% needs:  -GetValueOfFilename function (not version 02)
%         -getWavelengthAxis function
%         -korrelation function


%Wellenlänge automatisch (prüfen, ob alle Dateien selbes Gitter &
%Wellenlänge haben)

close all
clear all
load refSpectra2.mat
refSpectra = PLallShifted;
refSpectra = 0;
PLallShifted = [];

%% set Paths and Parameter
% path = 'Z:\Computer_Backups\Hyazinthenblau(LabP833)\Messungen\211123_CW14\Stark03'; %Z:\Computer_Backups\Hyazinthenblau(LabP833)\...
% measurement = 'CW14_A4_2'; %CW_12_1_01.....
path = 'Z:\Computer_Backups\Hyazinthenblau(LabP833)\Messungen\210705_SH358_Stark\210705_Stark08'; %Z:\Computer_Backups\Hyazinthenblau(LabP833)\...
measurement = 'CW8_D1_3'; %CW_12_1_01.....
% path = 'Z:\Computer_Backups\Hyazinthenblau(LabP833)\Messungen\210708\Dot15_Stark02'; %Z:\Computer_Backups\Hyazinthenblau(LabP833)\...
% measurement = 'CW8_D1_15'; %CW_12_1_01.....
% path = 'Z:\Computer_Backups\Hyazinthenblau(LabP833)\Messungen\210511_CW4_ErsteKondensatormessung\PLE'; %Z:\Computer_Backups\Hyazinthenblau(LabP833)\...
% measurement = 'CW4_01_03'; %CW_12_1_01.....
% path = 'Z:\Computer_Backups\Hyazinthenblau(LabP833)\Messungen\210708\Dot15_Stark01\'; %Z:\Computer_Backups\Hyazinthenblau(LabP833)\...
% measurement = 'CW8_D1_15_01'; %CW_12_1_01.....
% path = 'Z:\Computer_Backups\Hyazinthenblau(LabP833)\Messungen\210708\Dot16_Stark01\'; %Z:\Computer_Backups\Hyazinthenblau(LabP833)\...
% measurement = 'CW8_D1_16_01'; %CW_12_1_01.....

addPathEasy(path);

type = 2; %spec/kin = 1, starkNormal = 2, starkKin = 3, PLE = 4, saturation = 5, polarization = 6

% corrections
shiftRange = (400:1500);
frontCut = 1;
backCut = 0;
doShift='XCORR';
% gap Szie
gapSize = 0.000001;

intRange =  (400:600);

bgCorrection = false; %use recorded BG
baselineCorrection = true; %use bioinformatics tool for baseline correction
crayCorrection = true;
crayThreshold = 500;
loadWavelengthAxis = true;

% range find param
smRange = 30;
relativeShift = [ 100, 120 ];


%% get all relevant files of folder + import

files = dir(path);
%import parameters
for i=length(files):-1:1
    if isempty(strfind(files(i).name,measurement)) %check if "measurement" is part of name. Add more creteria if needed
        files(i)=[];
    else
        %        files(i).UNIXts = floor(posixtime(datetime(files(i).date,'TimeZone','local')));
        files(i).step = GetValueOfFilename(files(i).name,'Step');
        files(i).power = GetValueOfFilename(files(i).name,'uW');
        files(i).wpAngle = GetValueOfFilename(files(i).name,'DetLam');
        files(i).exWL = GetValueOfFilename(files(i).name,'ex');
        files(i).grating = GetValueOfFilename(files(i).name,'_g');
        files(i).centerWL = round(GetValueOfFilename(files(i).name,'_c')); %Mono gibt exakteren Wert aus
        files(i).temp = GetValueOfFilename(files(i).name,'K_');
        files(i).voltage = GetValueOfFilename(files(i).name,'Volt');
        files(i).Bfield = GetValueOfFilename(files(i).name,'B');
        files(i).tint = GetValueOfFilename(files(i).name,'s_');
    end
end
%% Sort files

switch type
    case 1
        [x,idx]=sort([files.UNIXts]);
    case 2
        [x,idx]=sort([files.voltage]);
    case 3
        [x,idx]=sort([files.UNIXts]);
    case 4
        [x,idx]=sort([files.exWL]);
    case 5
        [x,idx]=sort([files.power]);
    case 6
        [x,idx]=sort([files.wpAngle]);
    otherwise
        [x,idx]=sort([files.UNIXts]);
end
files=files(idx);
x=x';

%%Import files and do background correction
%(order remains like in files structure)

wbar = waitbar(0,'Importing Data Files');
wholeDataSet = [];
maxKinDim = 0;
for i=1:length(files)
    waitbar(i/length(files))
    tempData = fitsread([path '\' files(i).name]);
    if bgCorrection
        tempBG = fread(fopen([path '\BG\' files(i).name(1:size(files(i).name,2)-4) 'BINARYBG']),'int16','ieee-be');
        smoothedBG = smooth(tempBG,100);
        tempData = (tempData - smoothedBG')./files(i).tint; %Normalization to counts per second (cps)
    end
    %cosmic ray correction
    if crayCorrection
        for j=size(tempData,1):-1:1
            if max(tempData(j,:))>(mean(tempData(j,:))+crayThreshold)
                tempData(j,:)=[];   %discard single spectra with cosmic ray
            end
        end
    end

    data{i,1}=x(i);
    data{i,2}=tempData;

    nMeas = size(data,1);


    wholeDataSet=vertcat(wholeDataSet,tempData);
    if size(tempData,1)>maxKinDim
        maxKinDim = size(tempData,1);
    end
end
close(wbar)

%Load WavelengthAxis
if loadWavelengthAxis
    %check if grating and centerWavelength is same in all files
    grating = files(1).grating;
    centerWL = files(1).centerWL;
    for i=1:size(files,1)
        if files(i).grating ~= grating
            f = msgbox('Different gratings found!','Warning!')
            return
        elseif files(i).centerWL ~= centerWL
            f = msgbox('Different center wavelength found!','Warning!')
            return
        end
    end
    %load file with corresponding grating/centerWL
    gratingFile = ['Z:\Computer_Backups\Hyazinthenblau(LabP833)\Kalibrierung_' num2str(grating) 'er_Gitter\g' num2str(grating) '_c' num2str(centerWL) '.mat'];
    load(gratingFile);

end

%%baseline correction
if baselineCorrection
    pixels = 1:1600;
    Y = mean(wholeDataSet); %!!! Not normalized on Power
    y_bl_corr = msbackadj(pixels',Y');
    BL = Y-y_bl_corr';
    for i=1:size(data,1)
        data{i,3} = data{i,2}-BL;
    end
else
    for i=1:size(data,1)
        data{i,3} = data{i,2};
    end

end




%%shift correction of single kinetic measurements
%% Shifting kinetic cycles
wholeDataSetShifted=[];
meanShiftvector = zeros( nMeas,1 );
addSpectra = zeros( nMeas,1600 );
if (maxKinDim>1)
    for i=1:size(data,1)
        curSpecMat = data{i,3};
        % Smooth
        smFunc = 1/smRange*ones(1, smRange);
        smSpecMat = conv2(curSpecMat, smFunc, 'same');
        [maxVal, maxInd] = max(smSpecMat, [], 2);
        meanShiftvector(i) = round( mean(maxInd) );

        addSpectra(i, :) = mean(curSpecMat);
%         % diagnosis
%         figure();
%         plot(curSpecMat','b');
%         hold on
%         plot(maxInd,maxVal, 'rX',...
%                 'MarkerSize', 16)
%         xline( [maxInd - relativeShift(1); maxInd + relativeShift(2)] );
 

%         M_cell{1,1} = data{i,3}';
%         [shiftedSpectra, xcorr_sum, shiftVector] = shift(M_cell,doShift,shiftRange,refSpectra);
%         addSpectra(i, :) = mean(M_cell{1,1}');
%         stdDeviation(i)=std(shiftVector{1,1});
%         meanShiftvector(i)= round(mean(shiftVector{1,1}));
%         %data Spalte 3: korrellierte Daten, Spalte 4: Verschiebungsvektor
%         data{i,4} = shiftedSpectra{1,1};
%         data{i,5} = shiftVector{1,1};
%         data{i,6} = mean(data{i,4},2);
%         shiftedKinsCollection(:,i) = data{i,6};
%         wholeDataSetShifted = horzcat(wholeDataSetShifted,data{i,4}(:,:));
    end
%     relativeShift(1) = 50;
%     relativeShift(2) = 100; 
    ylim([0 10])
    %     lineBorders = zeros(size(wholeDataSet));
    arrayposition = 1;
    M_cellWithZeros = zeros(1600,50);
    for i=1:size(data,1)

        M_cell{1,1} = data{i,3}';
        lengthM_cell = size(M_cell{1,1});
        [shiftedSpectra, xcorr_sum, shiftVector] = shift(M_cell,doShift,meanShiftvector(i) - relativeShift(1):1:meanShiftvector(i) + relativeShift(2),refSpectra);
        stdDeviation(i)=std(WavelengthAxis(shiftVector{1,1},2));
        lineMean(arrayposition:arrayposition + lengthM_cell(2)) = meanShiftvector(i);
        lineBorders(arrayposition:arrayposition + lengthM_cell(2),1) = meanShiftvector(i)-relativeShift(1);
        lineBorders(arrayposition:arrayposition + lengthM_cell(2),2) = meanShiftvector(i)+relativeShift(2);
        arrayposition = arrayposition + lengthM_cell(2) ;
        %data Spalte 3: korrellierte Daten, Spalte 4: Verschiebungsvektor
        data{i,4} = shiftedSpectra{1,1};
        data{i,5} = shiftVector{1,1};
        data{i,6} = mean(data{i,4},2);
        shiftedKinsCollection(:,i) = data{i,6};
        M_cellWithZeros(meanShiftvector(i) - relativeShift(1):1:meanShiftvector(i) + relativeShift(2),i) = data{i,6}(meanShiftvector(i) - relativeShift(1):1:meanShiftvector(i) + relativeShift(2));
        wholeDataSetShifted = horzcat(wholeDataSetShifted,data{i,4}(:,:));

    end
    
%     figure
%     title('ShiftedKinsCollection')
%     imagesc(shiftedKinsCollection)
    %mean spectra of each kinetic collected in one Matrix
else
%     figure
%     plot(data{1,3}(1,:))
end

%% shift average spectra of all kinetics
M_cellWithZero{1,1} = M_cellWithZeros;
M_cell{1,1} = shiftedKinsCollection;

[shiftedKinSums,xcorr_sum,shiftVector] = shiftPairWise(M_cell,'XCORR',shiftRange,refSpectra);
allShifted = shiftedKinSums{1,1};
AddedM_cell{1,1} = addSpectra';
[addedShiftedKinSums,addedXcorr_sum,addedShiftVector] = shift(AddedM_cell,'XCORR',shiftRange,refSpectra);
addedAllShifted = addedShiftedKinSums{1,1};
% figure
% imagesc(allShifted)
% title('allShifted');
PLallShifted=mean(shiftedKinSums{1,1},2);

%%integrate intensity in integration range
for i=1:size(data,1)
    data{i,7} = sum(allShifted(intRange,i));
    specIntInt(i)=data{i,7};
end
specIntInt = specIntInt';
%% normalize and plot graphs according to measurement type

if true % show all data
    mainfigure = figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(2,2,1)
    imagesc(wholeDataSet)
    hold on
    scatter(lineBorders, 1:length(lineBorders),'MarkerEdgeColor',[0.7 0.4 0],...
        'MarkerFaceColor',[0.7 0.4 0],...
        'LineWidth',0.05)
    scatter(lineMean', 1:length(lineBorders),'MarkerEdgeColor',[1 0 0],...
        'MarkerFaceColor',[1 0 0],...
        'LineWidth',0.05)
    %     plot(meanShiftvector)
    xlabel('Pixel number')
    ylabel('Kinetic cycle')
    c = colorbar;
    ylabel(c,'Counts/s')
    title('raw dataset')
    line([shiftRange(1),shiftRange(1)], [0.5,size(wholeDataSet,1)+0.5], 'Color', 'r');
    line([shiftRange(size(shiftRange,2)),shiftRange(size(shiftRange,2))], [0.5,size(wholeDataSet,1)+0.5], 'Color', 'r');

    subplot(2,2,2)
    imagesc(allShifted')
    xlabel('Pixel number')
    ylabel('Kinetic cycle')
    c = colorbar;
    ylabel(c,'Counts/s')
    title('all spectra shifted')
    hold on
    line([intRange(1),intRange(1)], [0.5,size(allShifted,2)+0.5], 'Color', 'r');
    line([intRange(size(intRange,2)),intRange(size(intRange,2))], [0.5,size(allShifted,2)+0.5], 'Color', 'r');

    subplot(2,2,3)
    yyaxis left
    plot(PLallShifted)
    xlabel('Pixel number')
    ylabel('PL int. (Counts/s)')
    hold on
    if baselineCorrection
        yyaxis right
        plot(BL)
        ylabel('Baseline')
    end
    title('sum of all shifted spectra')

    subplot(2,2,4)
    plot(specIntInt(:,1))
    xlabel('kinetic cycle')
    ylabel('PL intensity - spec. int, (counts/s)')
    title('PL Intensity during measurement')
end
if type ==1
    mainfigure = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1)
    yyaxis left
    plot(WavelengthAxis(:,1),PLallShifted)
    xlabel('Wavelength (nm)')
    ylabel('PL int. (Counts/s)')
    hold on
    if baselineCorrection
        yyaxis right
        plot(WavelengthAxis(:,1),BL)
        ylabel('Baseline')
    end
    title('sum of all shifted spectra')
end

if type==2 %stark Plots (per measurement version)
    Voltage = x;
    for i=1:size(shiftedKinsCollection,2)
        normShiftedKins(:,i) = shiftedKinsCollection(:,i)./files(i).power;
    end
    mainfigure = figure('units','normalized','outerposition',[0 0 1 1]);
    %Normalize this on Power so intensity change by field is visible
    subplot(2,2,1)
    uimagesc(Voltage,WavelengthAxis(:,1),normShiftedKins)
    xlabel('Voltage')
    ylabel('Wavelength (nm)')
    title('StarkPlot')
    c = colorbar;
    ylabel(c,'Counts/(s*µW)')

    %Plot peakPositions
    peakPos = shiftVector{1,1};
    for i=1:length(peakPos)
        peakWL(i) = WavelengthAxis(peakPos(i),1);
        peakenergy(i) = WavelengthAxis(peakPos(i),2);
    end


    % derivation
    
    sigmaO = 0.000000000000100;
    field = (Voltage/gapSize) / 10000;
    field(5) = [];
    peakenergy(5) = [];
    stdDeviation(5) = [];
    f = field(frontCut:end-backCut);

    funPerfect = @(x1, f) x1(1) - x1(2).*(f-x1(3)).^2;
    funPerfect2 = @(x) x(1) - x(2).*(f-x(3)).^2 - peakenergy(frontCut:end-backCut);

    fitPerfect = lsqcurvefit(funPerfect, [2, 6.200126550657468e-04, 1500],f ,peakenergy(frontCut:end-backCut)');

    fitPerfect2 = lsqnonlin(funPerfect2, [2, 6.200126550657468e-04, -50]);

    fitPerfect3 = nlinfit(f, peakenergy(frontCut:end-backCut)', funPerfect, [2, 6.200126550657468e-01, 1500]);


    subplot(2,2,4) 
    plot(f,funPerfect(fitPerfect3,f),f,funPerfect(fitPerfect,f), f,peakenergy(frontCut:end-backCut)', 'o' )
    hold on
    plot(f,fitPerfect2(1) - fitPerfect2(2).*(f-fitPerfect2(3)).^2)

    dEdF = -2 * fitPerfect3(2).*(f - fitPerfect3(3));
    funSigma = @(sigmaF,f) 2 * fitPerfect3(2) .* f .* sigmaF - 2 * fitPerfect3(2) .* fitPerfect3(3) .* sigmaF;
    funSigmaJustLin = @(sigmaF,f) 2 * fitPerfect3(2) .* f .* sigmaF(1) -  sigmaF(2);
    sigmaFNonLin = nlinfit(f, stdDeviation(frontCut:end-backCut)', funSigma, 10);
    sigmaFNonJustLin = nlinfit(f, stdDeviation(frontCut:end-backCut)', funSigmaJustLin, [1, 1]);
    sigmaF = dEdF\stdDeviation(frontCut:end-backCut)';
%     fitWavelength = polyfit(field(frontCut:end-backCut),peakenergy(frontCut:end-backCut),2);
% 
% 
% 
%     q = polyder(fitWavelength);
%     derivativWavelength = polyval(q,field(frontCut:end-backCut));
%     fitFactors = polyfit(field(frontCut:end-backCut),stdDeviation(frontCut:end-backCut),1);
%     linearFit = polyval(fitFactors,field(frontCut:end-backCut));
%     fun = @(F) sqrt((polyval(q,field(frontCut:end-backCut))*F(1)).^2+sigmaO.^2) + F(2) - stdDeviation(frontCut:end-backCut);
% %     linearFit = field(frontCut:end-backCut)\stdDeviation(frontCut:end-backCut)';
%     x = lsqnonlin(fun,[10, 0.001]);
%     funCurvefit = @(F,derivativWavelength)sqrt((derivativWavelength*F(1)).^2+sigmaO.^2)+F(2); % https://de.mathworks.com/help/optim/ug/lsqcurvefit.html
%     x2 = lsqcurvefit(funCurvefit, [10, 0.001], derivativWavelength, stdDeviation(frontCut:end-backCut)');
%     FLinearFit = fitFactors(1)/q(1);


    peakPos=shiftVector{1,1};
    subplot(2,2,2)
    plot(Voltage,peakWL, 'o')
    xlabel('Voltage')
    ylabel('Wavelength (nm)')

    subplot(2,2,3)
    plot(field,stdDeviation, 'rx')
    hold on
    plot(field(frontCut:end-backCut),dEdF .* sigmaF)
    hold on
    plot(field(frontCut:end-backCut),2 .* fitPerfect3(2) .* f .* sigmaFNonJustLin(1) -sigmaFNonJustLin(2))
    plot(field(frontCut:end-backCut),2 .* fitPerfect3(2) .* f .* sigmaFNonLin - 2 .* fitPerfect3(2) .* fitPerfect3(3) .* sigmaFNonLin)

    %     plot(field(frontCut:end-backCut),sqrt((polyval(q,field(frontCut:end-backCut))*fitFactors).^2)+sigmaO.^2 , 'x')
%     plot(field(frontCut:end-backCut),linearFit , '.')
%     plot(field(frontCut:end-backCut),(sqrt((derivativWavelength*x(1)).^2+sigmaO.^2) + x(2)) , 'o')
%     plot(field(frontCut:end-backCut),(sqrt((derivativWavelength*x2(1)).^2+sigmaO.^2) + x2(2)) , 'x')
    ylabel('standard Deviation (eV)')

end


if type==3 %stark Plots (kinetic version)
    mainfigure = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1)
    kinStarkData = wholeDataSet';
    imagesc(8:-0.5:-8,WavelengthAxis(:,2),kinStarkData)

    for i=i:size(kinStarkData,2)
        plEnergy(i,1) = WavelengthAxis(data{1,5}(i),2);
    end
    subplot(2,2,2)
    plot(8:-0.5:-8,plEnergy,'x')
    xlabel('Voltage')
    ylabel('PL Energy (eV)')
end

if type == 4 %PLE
    for i=1:size(data,1)
        specIntInt(i,2) = specIntInt(i,1)/files(i).power;
    end
    mainfigure = figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(2,2,1)
    plot(x,specIntInt(:,2),'x');
    xlabel('excitation wavelength (nm)')
    ylabel('Intensity (counts/(s*µW))')
    title('PLE (intensity vs. ex. wavelength)')

end

if type==5 %saturation
    mainfigure = figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(2,2,1)
    plot(x,specIntInt(:,1),'x');
    xlabel('excitation power (µW)')
    ylabel('Intensity (counts/(s))')
    title('Saturation (intensity vs. ex. power)')
end

if type==6 %polarization
    %normalize specIntInt to power as well
    for i=1:size(data,1)
        specIntInt(i,2) = specIntInt(i,1)/files(i).power;
    end
    theta = x*2; %polarization turns twice as fast as lambda/2 plate

    mainfigure = figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(2,2,1)
    plot(theta,specIntInt(:,2),'x')
    xlabel('polarization angle (degree)')
    ylabel('Intensity (counts/(s*µW))')
    title('Polarization (intensity vs. angle)')

    subplot(2,2,2)
    polarscatter(deg2rad(theta),specIntInt(:,2))
    hold on
    polarscatter(deg2rad(theta+180),specIntInt(:,2)) %rotate copy
    theta2 = theta+180;
    %polarscatter(deg2rad(theta+180),fliplr(specIntInt(2,:))) %mirror copy
    %     xlabel('polarization angle (degree)')
    %     ylabel('Intensity (counts/(s*µW))')
    %     title('Polarization (intensity vs. angle)')


end
