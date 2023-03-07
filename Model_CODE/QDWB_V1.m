% ============================== BILAN =============================
% A quasi distributed model for calculation of Basin Water Balance
% This version is taylored only for Rokh-Neishabour basin.
clc;
clear all;
close all;
%% ========================== Initialization ========================
% Reading in Data and initialize Variables
LR = load('LR_Dat.txt', '-ascii'); % Laps Rate: 2 columns, 4375 days
Pr = load('Pre_Dat.txt', '-ascii'); % Precipit: 23 columns, 4375 days
Tm = load('Tem_Dat.txt', '-ascii'); % Temperat: 12 columns, 4375 days
% 12 columns = (3 * 4 stations), each 3 : Tmax, Tmin, Tave.
RE = load('RaD_Dat.txt', '-ascii'); % Ra+Eff_D: 3 columns, 366 days
% Ra, ayD, EffectiveDays;
Ra = RE(:, 1); % Now we have only Ra (366 * 1), and ...
eD = RE(:, 3); % Effective days for humidity estimation (366 * 1)
Kc = load('KIP_Dat.txt', '-ascii'); % 14 row, 16 col [LandUses * (12+1+1+1)]
Ir = Kc(:, 13); % Irrigation OR not for 14 LandUse types (14* 1)
In = Kc(:, 14); % Infilteration: Ratio of surface permeabality (14* 1)
Cc = Kc(:, 15); % Crop cover: Ratio of covered surface (14* 1)
RD = Kc(:, 16); % Root Depth for Crop cover (m) (14* 1)
Kc = Kc(:, 1:12); % 12 Monthly Kc values for all 14 LandUse types (14*12)
CT = load('Cr_Tem.txt', '-ascii'); % Tem Inv-Dist Coeff & dH (25*(4*2))
dH = CT(:, 5:8); % dH for 25 row (Regions) * 4 columns (stations)
CT = CT(:, 1:4); % CT for 25 row (Regions) * 4 columns (stations)
CP = load('Cr_Pre.txt', '-ascii'); % Pre Inv-Dist Coeff (25*23)
AW = load('Sol-AW.txt', '-ASCII'); % Top soil AW in % and mm (43* 2)
Annual_PCP = load('MATLAB_INPUT\PCP.txt', '-ASCII');
% ================== 25 Regions * 23 Rain stations ===================
% BasinWide Data;
% Note: Basin-Wide matrixes are (263 * 250) T
% After being masked they become col. vectors named BW?? (36946*1)
Mask = load('MATLAB_INPUT\mask_rn.txt', '-ascii');
Mask(Mask==0) = 1;
Mask(Mask==-9999) = 0;
Mask = logical(Mask);
dim = sum(sum(Mask)); % Numbers of CELLS
Am2 = 500*500; % Area of each CELL (m^2)
[dim1, dim2] = size(Mask);
% Initializing the Top Soil thickness equal to any value in mm
% This is the layer that contributes to Evaporation from bare soil surface
TopL = 50;
% Constructing Basin-Wide Fixed Matrixes ========================
% a) Soil WP & FC
T = load('MATLAB_INPUT\wp_prcnt_rn.txt', '-ascii');
BWWP = T.*(Mask); % Total WP in %
T = load('MATLAB_INPUT\fc_prcnt_rn.txt', '-ascii');
BWFC = T.*(Mask); % Total FC in %
% Error Trap for BAD data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
BWFC_SCLR = BWFC(Mask);
BWWP_SCLR = BWWP(Mask);
if (any(BWFC_SCLR==0) || any((BWFC_SCLR - BWWP_SCLR)<=0))
    disp(' !!!!!!!! ERROR !!!!!!!!')
    if any(BWFC==0)
        disp(find(BWFC==0));
        disp('BWFC==0');
    else
        disp(find(BWFC - BWWP)<=0);
        disp('(BWFC - BWWP)<=0');
    end
    disp(' !!!!!!!! ERROR !!!!!!!!')
    pause;
end
% b) Slope, Corrected CN & S
T = load('MATLAB_INPUT\cn_rn.txt', '-ascii');
BWCN = T.*(Mask);
T = load('MATLAB_INPUT\slp_prcnt_rn.txt', '-ascii');
BWSL = T.*(Mask)/100;
% Correction of CN for steep slopes
BWSL((BWSL<0.05)) = 0.05;
BWCN = BWCN .* ((323.57 + 15.63 .* (BWSL - 0.05)) ...
    ./ ((BWSL - 0.05) + 323.57));
% c) RSA, LandUse, Geo-Permeability
T = load('MATLAB_INPUT\rsa_14k_rn.txt', '-ascii');
BWRS = T.*(Mask); % Runoff Source
RSA_Count = sum(sum(BWRS));
T = load('MATLAB_INPUT\lu_14_rn.txt', '-ascii');
BWLU = T.*(Mask); % Land Use
T = load('MATLAB_INPUT\gp_0_09_rn.txt', '-ascii');
BWGP = T.*(Mask); % Geologigc Permeability
GP_sclr = BWGP(Mask);
% Removing bad data from GP map:
for m = 1:dim1
    for n = 1:dim2
        if BWGP(m,n)>0.3 && BWGP(m,n)<0.9
            BWGP(m,n) = 0.9;
        elseif BWGP(m,n)>0.2 && BWGP(m,n)<0.3
            BWGP(m,n) = 0.2;
        elseif BWGP(m,n)>0.15 && BWGP(m,n)<0.2
            BWGP(m,n) = 0.15;
        elseif BWGP(m,n)>0.1 && BWGP(m,n)<0.15
            BWGP(m,n) = 0.1;
        elseif BWGP(m,n)>0 && BWGP(m,n)<0.1
            BWGP(m,n) = 0;
        end
    end
end
% d) Basin Wide: Kc, Irr & In & Root Depth
BWIR = zeros(dim1, dim2);
BWIN = zeros(dim1, dim2);
BWCC = zeros(dim1, dim2);
BWKC = zeros(dim1, dim2, 12);
temp = zeros(dim1, dim2);
BWRD = zeros(dim1, dim2);
for m=1:dim1
    for n=1:dim2
        for k = 1:14 % 14 different Land Use
            if BWLU(m,n) == k
                BWIR(m,n) = Ir(k,1);
                BWIN(m,n) = In(k,1);
                BWCC(m,n) = Cc(k,1);
            else
                BWIR(m,n)= BWIR(m,n);
                BWIN(m,n)= BWIN(m,n);
                BWCC(m,n)= BWCC(m,n);
            end
            for mm = 1:12
                if BWLU(m,n) == k
                    BWKC(m,n,mm) = Kc(k,mm);
                else
                    BWKC(m,n,mm) =BWKC(m,n,mm);
                end
            end
        end
    end
end
IR_Count = sum(sum(BWIR));
for m=1:dim1
    for n=1:dim2
        for k = 1:14 % 14 different Land Use
            if BWLU(m,n) == k
                BWRD(m,n) = RD(k,1); % Calculating Root Depth for basin
            else
                BWRD(m,n) = BWRD(m,n);
            end
        end
    end
end
BWRD = BWRD .* 1000; % Convert Root Depth m to mm
BWRD_sclr = BWRD(Mask);
% e) Soil Moisture Thresholds
BWAW = BWFC - BWWP; % AW Percent
% Top Bare Soil Moisture Thresholds
TopL_WP = (BWWP .*0.01) .*TopL; % Wilting Point (mm)
TopL_FC = (BWFC .*0.01) .*TopL; % Field Capacity (mm)
TopL_AW = TopL_FC - TopL_WP; % AW (mm)
for m = 1:dim1
    for n = 1:dim2
        if BWRD(m,n) == 0
            TopL_WP(m,n) = 0;
            TopL_FC(m,n) = 0;
        end
    end
end
% Top Crop Covered Soil Moisture Thresholds
BWWP = BWRD .*BWWP .*0.01; % Wilting Point (mm)
BWFC = BWRD .*BWFC .*0.01; % Field Capacity (mm)
% Calculating AW in mm for the TOP SOIl (Contributes to soil evaporation)
% Total available water of this layer (Total Evaporable water): (BWTE)
BWTE = (BWAW./100) .* TopL; % Calculating: AW% *
% Actual water content for this layer (Available Evaporable water): (BWAE)
BWAE = zeros(dim1, dim2); % Initializing for the first day of
calculation
% SCS-CN: Limits for AMC (I and III)
CN3L = 0.9 .* TopL_AW; % TopSoil Moisture above this, indicates AMC III
CN1L = 0.3 .* TopL_AW; % TopSoil Moisture below this, indicates AMC I
%% ========================== Calculations ============================
sumD = 0;
SW = BWWP; % For the starting day, Soil Water is set to WP
BWAE = BWWP;
SW_bs = TopL_WP;
f = zeros(dim1, dim2);
Bilan_dim = zeros(12, 7);
Bilan_Am2 = zeros(12, 7);
for Year = 1998:2009 % ================= START of Year
    % NOTES:
    % 1- Leap Years: 2000, 2004, 2008
    % 2- Last year only has 357 days
    clc;
    disp('********** Processing **********');
    Bil = 1;
    save('Bil.txt', 'Bil', '-ascii', '-tabs');
    switch Year
        case {2000, 2004, 2008}
            lastD = 366;
            mm = [0 31 61 92 123 152 183 213 244 274 305 336 366];
        case 2009
            lastD = 357;
            mm = [0 31 61 92 123 151 182 212 243 273 304 335 365];
        otherwise
            lastD = 365;
            mm = [0 31 61 92 123 151 182 212 243 273 304 335 365];
    end
    % % Initialize SUM variables
    sumP = 0;
    sumR = 0;
    sumI = 0;
    sumETa = 0;
    sumDP = 0;
    sumIR = 0;
    BWsP = zeros(dim1, dim2);
    BWsR = zeros(dim1, dim2);
    BWsI = zeros(dim1, dim2);
    BWsETa = zeros(dim1, dim2);
    BWsDP = zeros(dim1, dim2);
    BWsIR = zeros(dim1, dim2);
    PCP = Annual_PCP((Year-1997),2);
    if PCP>300
        Irr_Rate = 0.7;
        IT_Rate = 0.3;
    elseif PCP>250 && PCP<300
        Irr_Rate = 0.5;
        IT_Rate = 0.25;
    elseif PCP<250
        Irr_Rate = 0.3;
        IT_Rate = 0.2;
    end
    BWIT = BWWP + (BWFC - BWWP) .*IT_Rate; % The THRESHOLD for Irrigation (mm)
    yearString = int2str(Year);
    h = waitbar(0, ['Year : ', yearString]);
    for ayD = 1:lastD % ================= START of ayD
        waitbar(ayD / lastD);
        R = zeros(dim1, dim2);
        I = zeros(dim1, dim2);
        ETa = zeros(dim1, dim2);
        DP = zeros(dim1, dim2);
        IR = zeros(dim1, dim2);
        month = find((mm<ayD), 1, 'last');
        % Preparing Daily Regional Data ==================================
        % Each day for 25 regions: Temp (Max and Ave), Rain & ETo
        dTa = -dH * LR((sumD + ayD), 2); % dTave, 25 * 4
        % Computation of mean of (Tmax - Tave) for each day
        difmax = Tm((sumD + ayD), [1, 4, 7, 10])...
            -Tm((sumD + ayD), [3, 6, 9, 12]);
        avemax = mean(difmax);
        difmin = Tm((sumD + ayD), [3, 6, 9, 12])...
            -Tm((sumD + ayD), [2, 5, 8, 11]);
        avemin = mean(difmin);
        rTa = dTa + repmat(Tm((sumD + ayD), [3, 6, 9, 12]), [25, 1]); % 25 * 4
        rTx = rTa + avemax; % 25 * 4
        rTn = rTa - avemin; % 25 * 4
        aTx = sum(rTx .* CT, 2); % Regional Tmax, 25 * 1
        aTa = sum(rTa .* CT, 2); % Regional Tave, 25 * 1
        aTn = sum(rTn .* CT, 2); % Regional Tmax, 25 * 1
        ETo = 0.0023 .* (aTa + 17.8) .* Ra(ayD) .* (aTx - aTn).^0.5; % 25 * 1
        aPr = CP * Pr((sumD + ayD), :)'; % 25 * 1
        % Constructing Basin-Wide Daily Matrixes =========================
        % a) Rain
        if sum(aPr) ~= 0
            BWPR = Basin_Wide(aPr);
            BWPR = BWPR.*(Mask);
        else
            BWPR = zeros(dim1,dim2);
        end
        % b) ETo
        BWET = Basin_Wide(ETo);
        BWET = BWET.*(Mask);
        %% =============== START of Calculations ================
        % (263 * 250)
        % % SCS-CN: correction for AMC (I or III)
        CN = BWCN;
        for m = 1:dim1
            for n = 1:dim2
                if SW_bs(m,n) >= CN3L(m,n) % Correction for AMC III
                    CN(m,n) = (23 .*CN(m,n)) ./(10 + 0.130 .*CN(m,n));
                elseif (SW_bs(m,n) <= CN1L(m,n)) % Correction for AMC I
                    CN(m,n) = (4.2 .*CN(m,n)) ./(10 - 0.058 .*CN(m,n));
                end
            end
        end
        % % Calculation of S for every where in the Basin
        BWSQ = (25400 ./ CN) - 254; % Calculates S value
        CN_SCLR = CN(Mask);
        if any(CN_SCLR<0.0003)
            disp('!!!!!!!! CN ERROR !!!!!!!!');
            pause;
        end
        T = ((BWPR - 0.2 * BWSQ) > 0); % To set R = 0, when P <= 0.2S (next line)
        R = ((BWPR - 0.2 * BWSQ).^2 ./ (BWPR + 0.8 * BWSQ)) .* BWRS .* T;
        I = (BWPR - R);
        % Correction of I (and R) for infilteration (surface permeability)
        T = I;
        for m = 1:dim1
            for n = 1:dim2
                if BWRS(m,n) == 1
                    I(m,n) = T(m,n) .*BWIN(m,n);
                    R(m,n) = R(m,n) + (T(m,n) .*(1 - BWIN(m,n)));
                end
            end
        end
        % Calculation of ETa
        fT = (SW - BWWP) ./ (BWFC - BWWP); % fT <= 1, Reducing
        Transpiration
        for m = 1:dim1
            for n = 1:dim2
                if BWRD(m,n) == 0
                    fT(m,n) = 0;
                end
            end
        end
        aT = fT .* ((BWKC(:,:, month) .* BWET) .* BWCC);
        kE = (SW_bs - TopL_WP) ./ TopL_AW;
        for m = 1:dim1
            for n = 1:dim2
                if BWRD(m,n) == 0
                    kE(m,n) = 0;
                end
            end
        end
        aE = kE .* (BWET .* (1 - BWCC));
        ETa = aT + aE; % Total actual Evapo-Transpiration: ETa
        T = SW - BWWP; % ETa can use soil moisture above WP. If ETa is
        for m = 1:dim1 % more it should be replaced with the
            correct value.
            for n = 1:dim2
                if ETa(m,n) > T(m,n)
                    ETa(m,n) = T(m,n);
                end
            end
        end
        % Updating SW :
        SW = SW + I - ETa; % Renew the SW for a new day
        for m = 1:dim1
            for n = 1:dim2
                if BWRD(m,n) == 0
                    SW(m,n) = 0;
                end
            end
        end
        T = TopL ./BWRD;
        SW_bs = SW_bs + (T .*I) - aE;
        % Correction for cells that have less moisture than WP
        for m = 1:dim1
            for n = 1:dim2
                if SW_bs(m,n) < TopL_WP(m,n)
                    SW_bs(m,n) = TopL_WP(m,n);
                end
            end
        end
        % Correcction of bare soil Moisture content in no-soil cells
        for m = 1:dim1
            for n = 1:dim2
                if BWRD(m,n) == 0
                    SW_bs(m,n) = 0;
                end
            end
        end
        % Error Trapping:
        T = (SW_bs - TopL_WP);
        Tsclr = T(Mask);
        if any(Tsclr<0)
            disp(find(Tsclr<0,1));
        end
        % Correction of SW :
        for m = 1:dim1
            for n = 1:dim2
                % Case a: Where / When : SW > FC
                if SW(m,n) > BWFC(m,n)
                    DP(m,n) = SW(m,n) - BWFC(m,n); % Above FC is DeepPerculation
                    SW(m,n) = BWFC(m,n); % Then, SW is set to FC
                    % Case b: Where / When : WP < SW < IT
                elseif (SW(m,n)>BWWP(m,n)) && (SW(m,n)<=BWIT(m,n)) && (BWIR(m,n)==1)
                    IR(m,n) = (BWFC(m,n) - BWWP(m,n)) .*Irr_Rate; %
                    Irrigation amount
                    SW(m,n) = SW(m,n) + IR(m,n); % After
                    irrigation SW reaches FC
                    if SW(m,n) > BWFC(m,n)
                        DP(m,n) = SW(m,n) - BWFC(m,n);
                        SW(m,n) = BWFC(m,n);
                    end
                end
            end
        end
        % Correction of DP (and R) for Geological Permeability
        T = DP;
        DP = T .* (1 - BWGP);
        R = R + (T .* BWGP);
        % % Keeping record of generated data : Annual_Lumped_Sum
        sumP = sumP + sum(BWPR(Mask));
        sumI = sumI + sum(I(Mask));
        sumR = sumR + sum(R(Mask));
        ETa_Mask = ETa(Mask);
        sumETa = sumETa + sum(ETa_Mask);
        sumDP = sumDP + sum(DP(Mask));
        sumIR = sumIR + sum(IR(Mask));
        % % Keeping record of generated data : Annual_Distrib_Sum
        BWsP = BWsP + BWPR;
        BWsR = BWsR + R;
        BWsI = BWsI + I;
        BWsETa = BWsETa + ETa;
        BWsDP = BWsDP + DP;
        BWsIR = BWsIR + IR;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T = sum(Pr((sumD + ayD), :))/23;
        Bil = [Year ayD T sum(aPr)/25 sum(BWPR)/36649 sum(ETa)/36649];
        save('Bil.txt', 'Bil', '-ascii', '-tabs', '-append');
    end % ================= END of ayD
    close(h)
    %% Saving OUTPUT files : Annual_Distrib_Sum ++++++++++++
    fName = ['BW_' int2str(Year)];
    T = [BWsP, BWsR, BWsI, BWsETa, BWsDP, BWsIR];
    save (fName, 'T', '-ASCII');
    sumD = sumD + lastD;
    Y = Year - 1997;
    % Creats two line for this year Bilan: 1) in mm & 2) in Mm^3
    Bilan_dim(Y, 2:7) = [sumP/dim, sumR/dim, sumI/dim, sumETa/dim,...
        sumDP/dim, sumIR/dim];
    Bilan_dim(Y, 1) = Year;
    Bilan_Am2(Y, 2:7) = [sumP sumR sumI sumETa sumDP sumIR].*(Am2/10^9);
    Bilan_Am2(Y, 1) = Year;
    %% Map Generation (Visulaizing Output)
    Annual Precipitation
    for m = 1:dim1
        for n = 1:dim2
            if Mask(m,n) == 0
                BWsP(m,n) = NaN;
            end
        end
    end
    ax1 = subplot(2,2,1);
    cnimage1 = imagesc(BWsP);
    figure(gcf); box on;
    hold on;
    colormap(ax1, flipud(winter))
    colorbar('location','EastOutside');
    set(cnimage1,'AlphaData',~isnan(BWsP)) % this line sets NaN values to white
    title('Annual Total Precipitation');
    axis off;
    axis(ax1,'square')

    % Annual Total ETa
    for m = 1:dim1
        for n = 1:dim2
            if Mask(m,n) == 0
                BWsETa(m,n) = NaN;
            end
        end
    end
    ax2 = subplot(2,2,2);
    cnimage2 = imagesc(BWsETa);
    figure(gcf); box on;
    hold on;
    colormap(ax2, flipud(hot))
    colorbar('location','EastOutside');
    set(cnimage2,'AlphaData',~isnan(BWsETa)) % this line sets NaN values to white
    title('Annual Total actual Evapo-Transpiration');
    axis off;
    axis(ax2,'square')

    % Annual Total Runoff
    for m = 1:dim1
        for n = 1:dim2
            if Mask(m,n) == 0
                BWsR(m,n) = NaN;
            end
        end
    end
    ax3 = subplot(2,2,3);
    cnimage3 = imagesc(BWsR);
    figure(gcf); box on;
    hold on;
    colormap(ax3, flipud(copper))
    colorbar('location','EastOutside');
    set(cnimage3,'AlphaData',~isnan(BWsR)) % this line sets NaN values to white
    title('Annual Total Runoff');
    axis off;
    axis(ax3,'square')

    % Annual Total Deep-Percolation
    for m = 1:dim1
        for n = 1:dim2
            if Mask(m,n) == 0
                BWsDP(m,n) = NaN;
            end
        end
    end
    ax4 = subplot(2,2,4);
    cnimage4 = imagesc(BWsDP);
    figure(gcf); box on;
    hold on;
    colormap(ax4, flipud(summer))
    colorbar('location','EastOutside');
    set(cnimage4,'AlphaData',~isnan(BWsDP)) % this line sets NaN values to white
    title('Annual Total Deep-Perculation');
    axis off;
    axis(ax4,'square')

    yearStringPlus = int2str(Year+1);
    suptitle(strcat('Dynamic Irrigation Threshold Scenario,', ' Year: '...
        , yearString, ' - ', yearStringPlus))

    Fig_Name = strcat('Figure_OutPut/QDWB_Dynamic_Pcp_Irr/',...
        yearString, '_QDWB_Dynamic_Pcp_Irr.png');
    print(gcf, Fig_Name,'-dpng','-r600');
end % ================= END of Year
disp('********** END **********'); %%%%%%%%%%%%%%%%%%%
% Output file, with 7 columns that contains:
% [Year Runoff Infilteration ETa DeepPercolation Irrigation]
% DeepPercolation: Groundwater recharge
% Irrigation: Part of Consumptive Use that was applied via irrigation
save('Bilan_dim.txt', 'Bilan_dim', '-ascii', '-tabs');
save('Bilan_Am2.txt', 'Bilan_Am2', '-ascii', '-tabs');
s = xlswrite('Bilan.xls', Bilan_dim );
