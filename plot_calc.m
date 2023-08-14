function plot_calc


    %% read the young objects sample - OLD
    
    % Reading the BTS sample:
    clear all;
    
    T = readtable('bts_sne_early_data.csv');
    T = T(:,{'ztf_name','name','filter','mag','mag_unc','jdobs','phase','disc_jd','limmag'});
    
    Flag = strcmp(T.filter,'b''g''');
    Filter(Flag) = 1; 
    Flag = strcmp(T.filter,'b''r''');
    Nt   = size(T,1);
    Filter(Flag) = 2; 
    Flag = strcmp(T.filter,'b''i''');
    Nt   = size(T,1);
    Filter(Flag) = 3; 
    T.filter = Filter(:);
    
    BTS_Bright = T;
    clear T;
    
    %% read the entire sample - OLD
    
    %Tcell2table(AAA(2:end,:));
    T=readtable('lc_table_bts.csv');
    T = T(:,{'ztf_name','name','filter','mag','mag_unc','jdobs','phase','disc_jd','limmag'});
    
    Nt   = size(T,1);
    Filter = nan(Nt,1);
    Flag = strcmp(T.filter,'g');
    Filter(Flag) = 1; 
    Flag = strcmp(T.filter,'r');
    Nt   = size(T,1);
    Filter(Flag) = 2; 
    Flag = strcmp(T.filter,'i');
    Nt   = size(T,1);
    Filter(Flag) = 3; 
    T.filter = Filter;
    
    BTS_All = T;
    clear T;
    
    %% merge BTS_All and BTS_Bright - OLD
    N_all    = size(BTS_All,1);
    N_bright = size(BTS_Bright,1);
    
    UnB = unique(BTS_Bright.name);
    UnA = unique(BTS_All.name);
    AddB = setdiff(UnB, UnA);  % names to add to the lisr
    Flag = false(size(BTS_Bright.name));
    for Ib=1:1:numel(AddB)
        Flag = Flag | strcmp(AddB{Ib},BTS_Bright.name);
    end
    
    T = [BTS_All; BTS_Bright(Flag,:)];
    N_bright = sum(Flag);
    
    %T = BTS_All;
    
    T = [BTS_All; BTS_Bright];
    %T.Sample = [zeros(N_all,1); ones(N_bright,1)];
    
    %% NEW
    
    %% load all BTS sample
    %cd data
    cd lcs_all_BTS_transients_20220626
    
    Files = dir('lcs*.csv');
    Nfile = numel(Files);
    for Ifile=1:1:Nfile
        Ifile
        %T = readtable('lcs_all_BTS_transients_20220626_0.csv');
        T = readtable(Files(Ifile).name);
        T = renamevars(T, {'ZTFID','IAUID','magerr','filterid'},{'ztf_name','name','mag_unc','filter'});
        T = sortrows(T, 'jdobs');

        if iscell(T.redshift)
            % make sure redshift is numeric
            T.redshift = cellfun(@str2double,T.redshift);
        end
        
        % correct for extinction
        Rv = 3.08;
        A_g = astro.spec.extinction(1,'ZTF','g',Rv);
        A_r = astro.spec.extinction(1,'ZTF','r',Rv);
        Ebv = T.A_V./3.08;
        Corr_g = Ebv.*A_g;
        Corr_r = Ebv.*A_r;
        Flag = T.filter == 1;
        T.mag(Flag) = T.mag(Flag) - Corr_g(Flag);
        Flag = T.filter == 2;
        T.mag(Flag) = T.mag(Flag) - Corr_r(Flag);
        
        if Ifile==1
            All = T;
        else
            All = [All; T];
        end
        size(All)
    end
    T = All;
    clear All;
    
    cd ..
    
    %% load Arcavi 2018
    
    TA = readtable('compiled_uv_optical_ir_ext_corrected_binned.txt');
    
    
    
    %% calc color and derivatives
    % search parameters
    clear Data_g
    clear Data_r
    clear Data_Type
    clear Data_Name
    
    MaxSep = 30./24;
    MinSep = 18./24;
    MaxErr = 0.1; %0.15; %0.05;
    
    
    UniqueNames = unique(T.ztf_name);
    Nun = numel(UniqueNames);
    
    Counter = 0;
    for Iun=1:1:Nun
        [Iun, Nun]
        % for each transient
        
        FlagTran = strcmp(T.ztf_name, UniqueNames{Iun});
        
        Tsn = T(FlagTran,:);
        
        Flag_g = Tsn.filter==1;
        Flag_r = Tsn.filter==2;
        Tsn_g  = Tsn(Flag_g,:);
        Tsn_r  = Tsn(Flag_r,:);
        
        Used_g = false(size(Tsn_g,1),1);
        Used_r = false(size(Tsn_r,1),1);
        
        Nep_g = size(Tsn_g,1);
        for Ig=1:1:Nep_g
            
            Flag_g = (Tsn_g.jdobs - Tsn_g.jdobs(Ig))>(-eps) & (Tsn_g.jdobs - Tsn_g.jdobs(Ig))<MaxSep & Tsn_g.mag_unc<MaxErr & Tsn_g.mag<Tsn_g.limmag;
            Flag_r = (Tsn_r.jdobs - Tsn_g.jdobs(Ig))>(-eps) & (Tsn_r.jdobs - Tsn_g.jdobs(Ig))<MaxSep & Tsn_r.mag_unc<MaxErr & Tsn_r.mag<Tsn_r.limmag;
            
            Range_g = range(Tsn_g.jdobs(Flag_g));
            Range_r = range(Tsn_r.jdobs(Flag_r));
            
            
            if sum(Flag_g)>=3 && sum(Flag_r)>=3 && Range_g>MinSep && Range_r>MinSep && ~any(Used_g(Flag_g)) && ~any(Used_r(Flag_r))
                Used_g(Flag_g) = true;
                Used_r(Flag_r) = true;
                
                % candidate found
                Counter = Counter + 1;
               
                MeanT_g = mean(Tsn_g.jdobs(Flag_g));
                
                % g
                T_g = Tsn_g.jdobs(Flag_g) - MeanT_g;
                %Par = polyfit(T_g , Tsn_g.mag(Flag_g), 1);
                Ng  = sum(Flag_g);
                H = [ones(Ng,1), T_g];
                [Par, ParErr] = lscov(H,Tsn_g.mag(Flag_g), 1./(Tsn_g.mag_unc(Flag_g).^2));
                %Resid_g = Tsn_g.mag(Flag_g) - polyval(Par,T_g);
                gdot     = Par(2);
                g_Mean   = Par(1);
                g_Err    = ParErr(2); %mean(Tsn_g.mag_unc(Flag_g));
                g_T      = MeanT_g;
                gdot_Err = ParErr(1); %std(Resid_g)./std(Tsn_g.jdobs(Flag_g));
                gdot_Err = max(gdot_Err, sqrt(gdot_Err.^2+0.015.^2));
                
                % r
                MeanT_r = MeanT_g; %mean(Tsn_r.jdobs(Flag_r));
                T_r = Tsn_r.jdobs(Flag_r) - MeanT_r;
                %Par = polyfit(T_r , Tsn_r.mag(Flag_r), 1);
                Nr  = sum(Flag_r);
                H = [ones(Nr,1), T_r];
                [Par, ParErr] = lscov(H,Tsn_r.mag(Flag_r), 1./(Tsn_r.mag_unc(Flag_r).^2));
                %Resid_r = Tsn_r.mag(Flag_r) - polyval(Par,T_r);
                rdot     = Par(2);
                r_Mean   = Par(1);
                r_Err    = ParErr(2); % mean(Tsn_r.mag_unc(Flag_r));
                r_T      = MeanT_r;
                rdot_Err = ParErr(1); %std(Resid_r)./std(Tsn_r.jdobs(Flag_r));
                rdot_Err = max(rdot_Err, sqrt(rdot_Err.^2+0.015.^2));

                Data_g(Counter, :) = [g_T, g_Mean, g_Err, gdot, gdot_Err, MeanT_g];
                Data_r(Counter, :) = [r_T, r_Mean, r_Err, rdot, rdot_Err, MeanT_r];
                Data_Name{Counter} = Tsn_g.name{1};
                Data_Type{Counter} = Tsn_g.type{1};
            end
        end
    end
    % number of measurments
   	Counter   
    % unique objects
    numel(unique(Data_Name))
    
    % unique types
    numel(unique(Data_Type))
    
    % for MaxErr=0.1
    % 15639, 588, 27
        
    %% generate [g-r, rdor, gdot]
        ColorMagDot = array2table([Data_g(:,1), Data_g(:,2) - Data_r(:,2), max(0.015, sqrt(Data_g(:,3).^2 + Data_r(:,3).^2)),Data_g(:,4), Data_g(:,5), Data_r(:,4), Data_r(:,5), Data_g(:,6), Data_r(:,6)]);

    ColorMagDot.Properties.VariableNames = {'mjd','gr','grErr','gdot','gdotErr','rdot','rdotErr','gJD','rJD'};
    ColorMagDot.Name  = Data_Name(:);
    ColorMagDot.type  = Data_Type(:);
    ColorMagDot.mjd = ColorMagDot.mjd - 2400000.5;
    ColorMagDot = ColorMagDot(:,{'Name','type','mjd','gr','grErr','gdot','gdotErr','rdot','rdotErr','gJD','rJD'});
    
    %ColorMagDotBigErr = ColorMagDot;
    
    %% extinction of AT2017gfo
    Ebv_GW = astro.spec.sky_ebv(197.450341598./RAD, -23.381467544./RAD);
    Eg_GW  = astro.spec.extinction(Ebv_GW,'SDSS','g');
    Er_GW  = astro.spec.extinction(Ebv_GW,'SDSS','r');
    Eu_GW  = astro.spec.extinction(Ebv_GW,'SDSS','u');
    Ei_GW  = astro.spec.extinction(Ebv_GW,'SDSS','i');
    EJ_GW  = astro.spec.extinction(Ebv_GW,'2MASS','J');
    EH_GW  = astro.spec.extinction(Ebv_GW,'2MASS','H');
    EK_GW  = astro.spec.extinction(Ebv_GW,'2MASS','K');
    
    
    %% generate [g-r, rdor, gdot] for Arcavi 2018
    % already extinction corrected
    Flag = strcmp(TA.filter,'g');
    TAg = TA(Flag,:);
    Flag = strcmp(TA.filter,'r');
    TAr = TA(Flag,:);
    
    % interpolate TAr to TAg
    TAri.mag   = interp1(TAr.phase, TAr.mag, TAg.phase);
    TAri.dmag  = interp1(TAr.phase, TAr.dmag, TAg.phase);
    TAri.phase = TAg.phase;
    
    ColorMagDotA.phase = TAg.phase;
    ColorMagDotA.gr    = TAg.mag - TAri.mag; % - (Eg_GW - Er_GW);
    ColorMagDotA.grErr = sqrt(TAg.dmag.^2 + TAri.dmag.^2);
    ColorMagDotA.gdot  = tools.deriv.numerical_deriv(ColorMagDotA.phase, TAg.mag);
    ColorMagDotA.rdot  = tools.deriv.numerical_deriv(ColorMagDotA.phase, TAri.mag);
    
    
    %% generate Table
    FID = fopen('Table_BTS.txt','w');    
    tools.table.fprintf(FID,'%20s  %10s  %9.3f   %5.2f %4.2f  %6.3f %5.3f  %6.3f %5.3f\n',ColorMagDot(:,1:9));
    
    fclose(FID);    
    
    %% Prep additional objects
    
    % add GW170817 to plot
    load LC.mat
    
    G_t    = [LC.t].';
    G_g    = [LC.g].' - Eg_GW;
    G_r    = [LC.r].' - Er_GW;
    G_u    = [LC.u].' - Eu_GW;
    G_i    = [LC.i].' - Ei_GW;
    G_J    = [LC.J].' - EJ_GW;
    G_H    = [LC.H].' - EH_GW;
    G_K    = [LC.K].' - EK_GW;
    
    G_gr   = G_g - G_r;
    G_gdot = tools.deriv.numerical_deriv(G_t, G_g);
    G_rdot = tools.deriv.numerical_deriv(G_t, G_r);

    G_idot = tools.deriv.numerical_deriv(G_t, [LC.i].');
    G_udot = tools.deriv.numerical_deriv(G_t, [LC.u].');
    G_Jdot = tools.deriv.numerical_deriv(G_t, [LC.J].');
    G_Hdot = tools.deriv.numerical_deriv(G_t, [LC.H].');
    G_Kdot = tools.deriv.numerical_deriv(G_t, [LC.K].');
    
    GW170817 = array2table([G_t, G_gr, G_gdot, G_rdot, G_udot, G_idot, G_Jdot, G_Hdot, G_Kdot,    G_u-G_g, G_r-G_i, G_J-G_H, G_H-G_K  ]);
    GW170817.Properties.VariableNames = {'t', 'gr', 'gdot', 'rdot',  'udot', 'idot','Jdot', 'Hdot', 'Kdot',  'ug','ri', 'JH', 'HK'};
    
    
    AddObjects = array2table([-0.2   0.5;1.9   0; -0.39  -0.65;0.2 -1.3]);
    AddObjects.Properties.VariableNames = {'gr','gdot'};
    AddObjects.Name = {'AT2018lqh','M85OT-1','SN2019ghp','SN2013fs'}';
    
    % print table of GW170817 data
    FID = fopen('Table_GW170817_Colors.txt','w');
    %tools.table.fprintf(FID, '%4.1f  %5.2f %5.2f %5.2f   %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n',GW170817);
    tools.table.fprintf(FID, '$%4.1f$ &  $%5.2f$ & $%5.2f$ & $%5.2f$ &   $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ \n',GW170817);
    fclose(FID);
    
    
    %% add models
    
    Pars.kM    = [0.3 0.8 0.6];
    Pars.ke    = [0.4 0.4 0.4];
    Pars.vM    = [0.15 0.15 0.15];
    Pars.Gamma = [0.6 0 0.2];
    Pars.Alpha = [0.7 0.7 0.5]
    Npar = numel(Pars.kM);
    
    t=logspace(log10(0.1),log10(2),10)'.*86400;
    for Ipar=1:1:Npar
    
        [L,Res]=astro.supernova.kilonova_multizone(t, 'kM',Pars.kM(Ipar), 'vM',Pars.vM(Ipar),'Gamma',Pars.Gamma(Ipar),'ke',Pars.ke(Ipar), 'Alpha',Pars.Alpha(Ipar));
        [Mag_g]=astro.spec.blackbody_mag_c(Res.T,'SDSS','g','AB',Res.rph);
        [Mag_r]=astro.spec.blackbody_mag_c(Res.T,'SDSS','r','AB',Res.rph);
        Mag_gdot = tools.deriv.numerical_deriv(t./86400, Mag_g);
        Mag_rdot = tools.deriv.numerical_deriv(t./86400, Mag_r);
        
        Model(Ipar).Table = array2table([t, Mag_g-Mag_r, Mag_gdot, Mag_rdot]);
        Model(Ipar).Table.Properties.VariableNames = {'t','gr','gdot','rdot'};
    end
    
    
    % plots
    %% select first observation of each transient
    UniqueName = unique(ColorMagDot.Name);
    Nun = numel(UniqueName);
    
    FlagNotSame  = false(size(ColorMagDot,1),1);
    % skip the first name as it is '-'
    for Iun=2:1:Nun
        IndN = find(strcmp(UniqueName{Iun}, ColorMagDot.Name));

        Diff = [1;diff(ColorMagDot.gJD(IndN))];
        
        FlagNotSame(IndN(Diff>0.5)) = true;

        if any(diff(ColorMagDot.gJD(IndN))<0)
            'a'
        end

    end

    ColorMagDot = ColorMagDot(FlagNotSame,:);


    % skip the first name as it is '-'
    
    FlagFirstObs = false(size(ColorMagDot,1),1);
    for Iun=2:1:Nun
        IndN = find(strcmp(UniqueName{Iun}, ColorMagDot.Name));

        [~,MinIjd] = min(ColorMagDot.gJD(IndN));
        %IndN(MinIjd);
        FlagFirstObs(IndN(MinIjd)) = true;
    end


    %% plot gdot histograms

    EdgesGdtot = (-2:0.2:2);
    X    = (EdgesGdtot(2:end) + EdgesGdtot(1:end-1)).*0.5;
    [N]  = histcounts(ColorMagDot.gdot,EdgesGdtot);
    [Nf] = histcounts(ColorMagDot.gdot(FlagFirstObs),EdgesGdtot);
    plot(X,N./sum(N),'k-','LineWidth',2);
    hold on
    plot(X,Nf./sum(Nf),'b-','LineWidth',2)
    set(gca,'YS','log')
    axis([-1 1 3e-4 1])

    H = xlabel('$\dot{g}$ [mag/day]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    H = ylabel('Number per 0.2 mag bin');
    H.FontSize = 18;
    H.Interpreter = 'latex';

    UniqueType = unique(ColorMagDot.type)
    Group{1} = {'AGN','AGN?'}
    Group{2} = {'CV','CV?','nova'}
    Group{3} = {'LBV','SLSN-I','SLSN-II','SN II','SN IIP','SN IIn','SN IIb','SN Ib','SN Ib-pec','SN Ib/c','SN Ibn','SN Ic','SN Ic-BL','TDE'}
    Group{4} = {'SN Ia','SN Ia-91T','SN Ia-CSM','SN Ia-SC','SN Ia-pec','SN Iax'}

    Color = colororder;
    for Igroup=1:1:numel(Group)
        FlagMember = ismember(ColorMagDot.type,Group{Igroup});
        [Nmem]  = histcounts(ColorMagDot.gdot(FlagMember),EdgesGdtot);
        plot(X,Nmem./sum(Nmem),'b--','LineWidth',1,'Color',Color(Igroup,:))
    end

    legend('All','first','AGN','CV','SN','SN Ia')    

    print gdot_hist.eps -depsc2

    %% table of types numbers

    FID = fopen('Table_ObjTypes.txt','w');
    Ntype = numel(UniqueType);
    for Itype=1:1:Ntype
        FlagTT = strcmp(ColorMagDot.type,UniqueType{Itype});
        NuniqueObj = numel(unique(ColorMagDot.Name(FlagTT)));

        fprintf(FID, '%20s &  %5d & %5d\\\\ \n',UniqueType{Itype}, sum(FlagTT), NuniqueObj)
        
    end


    %% sources with g-dot<-0.8 and gdot>0.8
    IndFast = find(ColorMagDot.gdot<-0.8)
    IndFast = find(ColorMagDot.gdot>0.8)
    
    Iobj = find(strcmp(ColorMagDot.Name, 'AT2018ctl'))   % CV
    Iobj = find(strcmp(ColorMagDot.Name, 'AT2019cuy'))   % AGN

    plot(ColorMagDot.gJD(Iobj)-2450000,ColorMagDot.gdot(Iobj),'k.')


    ExtremeObj = {0.8, -1.57, 'AT2019cuy','AGN';...
     -0.02, -1.87, 'AT2018ctl','CV';...
     -0.16,-0.98, 'AT2021hoz','CV';...
     0.48, -0.91, 'AT2018cch','AGN';...
     0.93,-0.81, 'AT2018ief','AGN';...
     0.25,0.88,  'AT2019cmi','AGN'};


    %% g-r vs. gdot
    
    % plot BT sample
    Color = 'k';
    PlotType = 'e1';
    switch PlotType
        case 'e'
            plot.errorxy([ColorMagDot.gr, ColorMagDot.gdot,ColorMagDot.grErr, ColorMagDot.gdotErr],...
                        'ColX',1, 'ColY',2, 'ColXe',3, 'ColYe',4, 'EdgeColor',Color, 'FaceColor',Color, 'MarkSize',4);
        case 'e1'
            %plot.errorxy1([ColorMagDot.gr, ColorMagDot.gdot,ColorMagDot.grErr, ColorMagDot.gdotErr],'MarkerSize',3);
            
            FlagIa = strcmp(Data_Type,'SN Ia');
            FlagIa = FlagFirstObs;
            plot.errorxy1([ColorMagDot.gr(~FlagIa), ColorMagDot.gdot(~FlagIa),ColorMagDot.grErr(~FlagIa), ColorMagDot.gdotErr(~FlagIa)],'MarkerSize',3);
            hold on;
            ColorIa = [1 0.5 0];
            plot.errorxy1([ColorMagDot.gr(FlagIa), ColorMagDot.gdot(FlagIa),ColorMagDot.grErr(FlagIa), ColorMagDot.gdotErr(FlagIa)],...
                'MarkerSize',2,'MarkerFaceColor',ColorIa,'MarkerEdgeColor',ColorIa,'LineColor',ColorIa);

        case 'p'
            H = plot(ColorMagDot.gr, ColorMagDot.gdot,'k.');
            H.Color = [0.8 0.8 0.8];
        otherwise
            error('Unknown PlotType option');
    end
    hold on;
    
    %
    
    % 18gep line
    Igep = find(strcmp(ColorMagDot.Name,'SN2018gep'));
    plot(ColorMagDot.gr(Igep),ColorMagDot.gdot(Igep), 'r-')
    
    % plot extinction arrow
    Egr = astro.spec.extinction(0.2,'ZTF','g',3.08)-astro.spec.extinction(0.2,'ZTF','r',3.08);
    %Hq=quiver(1.5,-2, Egr, 0,'LineWidth',2);
    %Hq.MaxHeadSize=10;
    %Hq.Color = 'm';
    
    axis([-0.6 2 -3 2.5]);

    
    % plot GW170817               
    hold on;
    plot(GW170817.gr, GW170817.gdot,'k-','Color',[0.8 0.8 0.8], 'LineWidth',2);
    scatter(GW170817.gr, GW170817.gdot, 50, GW170817.t, 'filled');
    H = colorbar;
    H.Label.String = 'Age [days]';
    H.Label.Interpreter = 'latex';
    % time labels
    for It=1:1:numel(GW170817.t)
        if GW170817.t(It)<3.6
            Ht = text(GW170817.gr(It)+0.02, GW170817.gdot(It)+0.2, sprintf('%3.1f',GW170817.t(It)));
            Ht.Interpreter = 'latex';
        end
    end
    
    Xa = 1.5;
    Ya = -2;
    [NX,NY]=plot.coo2normalized([Xa, Xa+Egr],[Ya, Ya]);
    Ha=annotation('arrow',NX, NY,'Linewidth',2)

    
    % plot Arcavi 2018
    plot(ColorMagDotA.gr, ColorMagDotA.gdot, 'k-','Color',[0.8 0.8 0.8], 'LineWidth',1);
    scatter(ColorMagDotA.gr, ColorMagDotA.gdot, 50, ColorMagDotA.phase)
    
    % plot models
    %Color = plot.generate_colors(Npar);
    %for Ipar=1:1:Npar
    %    Hmodel = plot(Model(Ipar).Table.gr, Model(Ipar).Table.gdot,'k-');
    %    Hmodel.Color = Color(Ipar,:);
    %end
    
    % plot special objects
    Nobj = size(AddObjects,1);
    Symbol = {'gs','rs','g^','gp'};
    for Iobj=1:1:Nobj
        plot(AddObjects.gr(Iobj), AddObjects.gdot(Iobj),Symbol{Iobj}, 'MarkerSize',10, 'MarkerFaceColor',Symbol{Iobj}(1));
    end
    
    H = xlabel('$g-r$ [mag]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    H = ylabel('$\dot{g}$ [mag/day]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    
    axis([-0.6 2 -3 2.5]);

    
    for Ieo=1:1:size(ExtremeObj,1)
        Htt = text(ExtremeObj{Ieo,1}+0.02,ExtremeObj{Ieo,2}-0.12,ExtremeObj{Ieo,4});
        Htt.FontSize = 8;

    end
    
    print GW_gr_gdot.eps -depsc2
    print GW_gr_gdot.jpg -djpeg100
    
    
    %% gdot vs. rdot
    
    % plot BT sample
    Color = 'k';
    %plot.errorxy([ColorMagDot.rdot, ColorMagDot.gdot,ColorMagDot.rdotErr, ColorMagDot.gdotErr],...
    %                    'ColX',1, 'ColY',2, 'ColXe',3, 'ColYe',4, 'EdgeColor',Color, 'FaceColor',Color, 'MarkSize',4);
    
    %plot.errorxy1([ColorMagDot.rdot, ColorMagDot.gdot,ColorMagDot.rdotErr, ColorMagDot.gdotErr],'MarkerSize',3);
    plot.errorxy1([ColorMagDot.rdot(~FlagIa), ColorMagDot.gdot(~FlagIa),ColorMagDot.rdotErr(~FlagIa), ColorMagDot.gdotErr(~FlagIa)],'MarkerSize',3);
    hold on;
    plot.errorxy1([ColorMagDot.rdot(FlagIa), ColorMagDot.gdot(FlagIa),ColorMagDot.rdotErr(FlagIa), ColorMagDot.gdotErr(FlagIa)],...
        'MarkerSize',2,'MarkerFaceColor',ColorIa,'MarkerEdgeColor',ColorIa,'LineColor',ColorIa);
    
    %plot(ColorMagDot.rdot, ColorMagDot.gdot,'k.');
    
    % 18gep line
    Igep = find(strcmp(ColorMagDot.Name,'SN2018gep'));
    plot(ColorMagDot.rdot(Igep),ColorMagDot.gdot(Igep), 'r-')
    
    % plot GW170817               
    hold on;
    plot(GW170817.rdot, GW170817.gdot,'k-','Color',[0.8 0.8 0.8], 'LineWidth',2);
    scatter(GW170817.rdot, GW170817.gdot, 50, GW170817.t, 'filled');
    H = colorbar;
    H.Label.String = 'Age [days]';
    H.Label.Interpreter = 'latex';
    % time labels
    for It=1:1:numel(GW170817.t)
        if GW170817.t(It)<3.6
            Ht = text(GW170817.rdot(It)+0.02, GW170817.gdot(It)+0.2, sprintf('%3.1f',GW170817.t(It)));
            Ht.Interpreter = 'latex';
        end
    end
    
    % plot Arcavi 2018
    plot(ColorMagDotA.rdot, ColorMagDotA.gdot, 'k-','Color',[0.8 0.8 0.8], 'LineWidth',1);
    scatter(ColorMagDotA.rdot, ColorMagDotA.gdot, 50, ColorMagDotA.phase)
    
    % plot models
    Color = plot.generate_colors(Npar);
    for Ipar=1:1:Npar
        Hmodel = plot(Model(Ipar).Table.rdot, Model(Ipar).Table.gdot,'k-');
        Hmodel.Color = Color(Ipar,:);
    end
    
    H = xlabel('$\dot{r}$ [mag/day]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    H = ylabel('$\dot{g}$ [mag/day]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    
    axis([-3.2 1.5 -3.2 2.5])
        
    
    print GW_rdot_gdot.eps -depsc2
    print GW_rdot_gdot.jpg -djpeg100
    
    %% PCA / g-r vs. rdot vs. gdot
    
    GW_Data = [GW170817.gr, GW170817.gdot, GW170817.rdot];
    
    [Coef, GW_PCA, Latent, TS, Ex, Mu] = pca([GW170817.gr, GW170817.gdot, GW170817.rdot])
    
    
    Data    = [ColorMagDot.gr, ColorMagDot.gdot, ColorMagDot.rdot];
    DataErr = Data + [ColorMagDot.grErr, ColorMagDot.gdotErr, ColorMagDot.rdotErr];
    
    Coef = Coef.';
    
    Proj    = Data*Coef;
    ProjErr = DataErr*Coef;
    ProjErr = abs(ProjErr - Proj);
    
    Color = 'k';
    %plot.errorxy([Proj(:,1), Proj(:,2), ProjErr(:,1), ProjErr(:,2)],...
    %        'ColX',1, 'ColY',2, 'ColXe',3, 'ColYe',4, 'EdgeColor',Color, 'FaceColor',Color, 'MarkSize',4);
        
    %plot.errorxy1([Proj(:,1), Proj(:,2), ProjErr(:,1), ProjErr(:,2)],'MarkerSize',3);
    plot.errorxy1([Proj(~FlagIa,1), Proj(~FlagIa,2), ProjErr(~FlagIa,1), ProjErr(~FlagIa,2)],'MarkerSize',3);
    hold on;
    plot.errorxy1([Proj(FlagIa,1), Proj(FlagIa,2), ProjErr(FlagIa,1), ProjErr(FlagIa,2)],...
        'MarkerSize',2,'MarkerFaceColor',ColorIa,'MarkerEdgeColor',ColorIa,'LineColor',ColorIa);

    %plot(Proj(:,1), Proj(:,2),'k.');
    hold on;
    
    % 18gep line
    Igep = find(strcmp(ColorMagDot.Name,'SN2018gep'));
    plot(Proj(Igep,1),Proj(Igep,2), 'r-')
    
    % plot extinction arrow
    Egr = astro.spec.extinction(0.2,'ZTF','g',3.08)-astro.spec.extinction(0.2,'ZTF','r',3.08);
    DataQ = [0 0 0; Egr 0 0];
    DataQPCA = DataQ*Coef;
    DataQPCA = DataQPCA(2,:) - DataQPCA(1,:);
    
    %Hq=quiver(-3,1.2, DataQPCA(1), DataQPCA(2),'LineWidth',2);
    %Hq.MaxHeadSize=1;
    %Hq.Color = 'm';

    axis([-5 2.3 -1.2 1.6])

    Xa = -3;
    Ya = 1;
    [NX,NY]=plot.coo2normalized([Xa, Xa+DataQPCA(1)],[Ya, Ya+DataQPCA(2)]);
    Ha=annotation('arrow',NX, NY,'Linewidth',2)

    
    
    plot(GW_PCA(:,1), GW_PCA(:,2),'k-','Color',[0.8 0.8 0.8],'LineWidth',2)
    scatter(GW_PCA(:,1), GW_PCA(:,2), 50, GW170817.t, 'filled');
    H = colorbar;
    H.Label.String = 'Age [days]';
    H.Label.Interpreter = 'latex';
    % time labels
    for It=1:1:numel(GW170817.t)
        if GW170817.t(It)<3.6
            Ht = text(GW_PCA(It,1)+0.02, GW_PCA(It,2)+0.05, sprintf('%3.1f',GW170817.t(It)));
            Ht.Interpreter = 'latex';
        end
    end
    
    % plot Arcavi 2018
    ArcaviPCA    = [ColorMagDotA.gr, ColorMagDotA.gdot, ColorMagDotA.rdot]*Coef;
    plot(ArcaviPCA(:,1),ArcaviPCA(:,2), 'k-','Color',[0.8 0.8 0.8], 'LineWidth',1);
    scatter(ArcaviPCA(:,1), ArcaviPCA(:,2), 50, ColorMagDotA.phase)
    
    % plot models
    Color = plot.generate_colors(Npar);
    for Ipar=1:1:Npar
        GW_Model = [Model(Ipar).Table.gr, Model(Ipar).Table.gdot, Model(Ipar).Table.rdot];
        GW_ModelPCA = GW_Model*Coef;
        
        Hmodel = plot(GW_ModelPCA(:,1), GW_ModelPCA(:,2),'k-');
        Hmodel.Color = Color(Ipar,:);
    end
    
    H = xlabel('PC1');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    H = ylabel('PC2');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    
    
    axis([-5 2.3 -1.2 1.6])
    
    print GW_PCA.eps -depsc2
    print GW_PCA.jpg -djpeg100
    
    
    %% probability distribution
    Step_gr   = 0.1;
    Step_gdot = 0.1;
    Step_rdot = 0.1;
    
    Vec_gr   = [-1:Step_gr:2];
    Vec_gdot = [-3.5:Step_gdot:2.5];
    Vec_rdot = [-3.5:Step_rdot:2.5];
    N_gr     = numel(Vec_gr);
    N_gdot   = numel(Vec_gdot);
    N_rdot   = numel(Vec_rdot);
    
    mean(ColorMagDot.grErr)
    mean(ColorMagDot.gdotErr)
    Sigma_gr = 0.1;
    Sigma_gdot = 0.05;
    Sigma_rdot = 0.05;
    
    DataPos      = [ColorMagDot.gr, ColorMagDot.gdot, ColorMagDot.rdot];
    DataPosFirst = [ColorMagDot.gr(FlagFirstObs), ColorMagDot.gdot(FlagFirstObs), ColorMagDot.rdot(FlagFirstObs)];

    CovMat  = diag([Sigma_gr, Sigma_gdot, Sigma_rdot]);
    Npos    = size(DataPos,1);
    
    Cube = zeros(N_gr, N_gdot, N_rdot);
    Coo1 = zeros(N_gr, N_gdot, N_rdot);
    Coo2 = zeros(N_gr, N_gdot, N_rdot);
    Coo3 = zeros(N_gr, N_gdot, N_rdot);
    for Igr=1:1:N_gr
        [Igr, N_gr]
        for Igdot=1:1:N_gdot
            for Irdot=1:1:N_rdot
                gr   = Vec_gr(Igr);
                gdot = Vec_gdot(Igdot);
                rdot = Vec_rdot(Irdot);
                
                AllP      = mvnpdf(DataPos, [gr, gdot, rdot], CovMat);
                AllPfirst = mvnpdf(DataPosFirst, [gr, gdot, rdot], CovMat);

                Cube(Igr, Igdot, Irdot) = sum(AllP);
                CubeFirst(Igr, Igdot, Irdot) = sum(AllPfirst);
                Coo1(Igr, Igdot, Irdot) = gr;
                Coo2(Igr, Igdot, Irdot) = gdot;
                Coo3(Igr, Igdot, Irdot) = rdot;
                
            end
        end
    end
    
    % sanity check - should be 1
    sum(Cube,'all').*prod([Step_gr, Step_gdot, Step_rdot])./Npos
    Cube = Cube.*prod([Step_gr, Step_gdot, Step_rdot])./Npos;
    
    %% generate Table of probabilities in phase space
    
    Max  = max(Cube,[],'all');
    AllP = AllP./Max;

    MaxFirst  = max(CubeFirst,[],'all');
    AllPfirst = AllPfirst./MaxFirst;

    Flag = Cube>(1./1000) | CubeFirst>(1./1000)
    CubeFlag = Cube;
    CubeFlag(~Flag) = 0;
    
    FID = fopen('Table_SparseProb.txt','w');
    fprintf(FID, '%5.1f %5.1f %5.1f  %5.3e  %5.3e\n',[Coo1(Flag(:)), Coo2(Flag(:)), Coo3(Flag(:)), Cube(Flag(:)), CubeFirst(Flag(:))].');
    %fprintf(FID, '$%5.1f$ & $%5.1f$ & $%5.1f$ &  $%5.3e$  $%5.3e$\n',[Coo1(Flag(:)), Coo2(Flag(:)), Coo3(Flag(:)), Cube(Flag(:)), CubeFirst(Flag(:))].');
    fclose(FID)
    
    
    %% gdot rdot for transients types
    UnType = unique(ColorMagDot.type);
    Ntype  = numel(UnType);
    for Itype=1:1:Ntype
        Flag = strcmp(ColorMagDot.type, UnType{Itype});
        plot(ColorMagDot.gdot(Flag), ColorMagDot.gr(Flag), '.');
        hold on;
    end
    
    Flag = strcmp(ColorMagDot.type, 'CV');
    plot(ColorMagDot.gdot(Flag), ColorMagDot.gr(Flag), 'k.','MarkerSize',10);
    
    
    
    
    
    %% END
    