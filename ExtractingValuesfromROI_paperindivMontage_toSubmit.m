clear
close all
clc

MainFolder = uigetdir([], 'choose main directory with all subjects'); %Main directory with all subjects


subjects = {''};

load('targetCoord','targetCoord')%load target coordinates

montageFolders = {'simulations_Func_0071_V2','simulations_FC6_TP8-P8_FC5_TP7-P7','simulations_T7-T8_ring2'};
E_norm = nan(length(subjects),3);
E_normal = nan(length(subjects),3);
FOC50 = nan(length(subjects),3);
FOC75 = nan(length(subjects),3);
E_normal_FuncroiL = nan(length(subjects),3);
E_normal_FuncroiR = nan(length(subjects),3);
E_norm_FuncroiL = nan(length(subjects),3);
E_norm_FuncroiR = nan(length(subjects),3);
E_norm_AnatroiR = nan(length(subjects),3);
E_norm_AnatroiL = nan(length(subjects),3);

for subj=1:length(subjects)
    %loop through montages
    for montageNr = 1:length(montageFolders)
        close all force
        currSubj = subjects{subj};
        disp(currSubj)
        
        %cd to subj dir
        cd(fullfile(MainFolder,currSubj))
        
        %%%%%%%%%%%%----------------------------------------------------------------------------
        % load e-field simulations from the individually optimized montage
        fname =[subjects{subj} '_TDCS_1_scalar_fsavg.msh'];
        fname2 =[subjects{subj} '_TDCS_1_scalar.msh'];
        
        pname = fullfile (MainFolder,subjects{subj},montageFolders{montageNr},'fsavg_overlays');
        pname2 = fullfile (MainFolder,subjects{subj},montageFolders{montageNr});
        disp(['........Working on folder: ' pname])
        
        if isequal(fname,0) || isequal(pname,0); return; end
        
        % Load simulations
        m  = mesh_load_gmsh4(fullfile(pname,fname));
        m2 = mesh_load_gmsh4(fullfile(pname2,fname2));
        
        % Show E-normal and E-norm from surface and volume data
        maxScaleEnorm   = max(m.node_data{2,1}.data);
        maxScaleEnormal = max(m.node_data{3,1}.data);
        
        mesh_show_surface(m,'field_idx','E_norm','scaleLimits',[-1*maxScaleEnorm maxScaleEnorm])
        title([subjects{subj} ' E_norm IndivMontage']);
        
        mesh_show_surface(m,'field_idx','E_normal','scaleLimits',[-1*maxScaleEnormal maxScaleEnormal])
        title([subjects{subj} ' E_normal IndivMontage']);
        
        mesh_show_surface(m2,'field_idx','normE','scaleLimits',[-1*maxScaleEnorm maxScaleEnorm])
        title([subjects{subj} ' E_norm IndivMontage']);
        
        mesh_show_surface(m,'field_idx','E_angle','scaleLimits',[])
        title([subjects{subj} ' E_angle IndivMontage']);
        
        % Note: Focality results are only an estimate, as the fsaverage surface
        % rather the individual surface was used to determine the stimulated areas
        disp(' ')
        disp('whole cortex:')
        
        summary1 = mesh_get_fieldpeaks_and_focality(m,'field_idx','E_norm');
        summary2 = mesh_get_fieldpeaks_and_focality(m,'field_idx','E_normal');
        
        FieldVal.Enorm(subj,montageNr).p95      = summary1.perc_values(1);
        FieldVal.Enorm(subj,montageNr).p99      = summary1.perc_values(2);
        FieldVal.Enorm(subj,montageNr).p999     = summary1.perc_values(3);
        FieldVal.Enorm(subj,montageNr).peak_pos = summary1.XYZ_perc;
        
        FieldVal.FOC(subj,montageNr).c50 = summary1.focality_values(1);
        FieldVal.FOC(subj,montageNr).c75 = summary1.focality_values(2);
        
        FieldVal.Enormal(subj,montageNr).p95     = summary2.perc_values(1);
        FieldVal.Enormal(subj,montageNr).p99     = summary2.perc_values(2);
        FieldVal.Enormal(subj,montageNr).p999    = summary2.perc_values(3);
        FieldVal.EnormalNeg(subj,montageNr).p95  = summary2.perc_neg_values(1);
        FieldVal.EnormalNeg(subj,montageNr).p99  = summary2.perc_neg_values(2);
        FieldVal.EnormalNeg(subj,montageNr).p999 = summary2.perc_neg_values(3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ROI analysis
        %
        % the "summary" structure contains the peak values
        % together with their positions:
        %   "summary.percentiles" lists the tested percentile cutoffs -
        %   the 99.9 percentile is the 3rd entry
        %
        %   "summary.perc_values" lists the corresponding values, as they
        %   are also displayed by mesh_get_fieldpeaks_and_focality
        %
        %   "summary.XYZ_perc" lists the corresponding center positions -
        %   the center position for the 99.9 percentile is the 3rd row:
        %  peak_pos=summary1.XYZ_perc(3,:);
        
        
        %%%%%%%%%%%ONLY FOR VISUALIZATION PURPOSES
        peak_posLo = targetCoord(subj).LH;%[-56 -4 -10];
        peak_posRo = targetCoord(subj).RH;%[56 -4 -10];
        % % distance to peak position
        distRo = sqrt(sum(bsxfun(@minus,m.nodes,peak_posRo).^2,2));
        %
        % % extract nodes closer than 3 mm, and the related triangles and data
        node_idxRo = distRo<7;
        m_ROIRo    = mesh_extract_regions(m, 'node_idx', node_idxRo);
        
        distLo=sqrt(sum(bsxfun(@minus,m.nodes,peak_posLo).^2,2));
        
        node_idxLo = distLo<7;
        m_ROILo    = mesh_extract_regions(m, 'node_idx', node_idxLo);
        try
            % % show the extracted ROI:
            % % show the whole cortex semi-transparent
            mesh_show_surface(m,'showSurface',true,'facealpha',0.3);
            % % add the extracted area to the plot
            mesh_show_surface(m_ROIRo,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca);
            mesh_show_surface(m_ROILo,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca);
        catch
        end
        
        % % extract nodes closer than 10 mm, and the related triangles and data
        node_idxR = distRo<7;
        m_ROIR    = mesh_extract_regions(m, 'node_idx', node_idxR);
        
        node_idxL = distLo<7;
        m_ROIL    = mesh_extract_regions(m, 'node_idx', node_idxL);
        
        % % show the extracted ROI:
        % % show the whole cortex semi-transparent
        mesh_show_surface(m,'showSurface',true,'facealpha',0.3);
        % % add the extracted area to the plot
        mesh_show_surface(m_ROIR,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca);
        mesh_show_surface(m_ROIL,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca);
        
        %
        % % get some key results for the spherical ROI
        disp(' ')
        disp('10 mm spherical ROI around peak positions:')
        
        summary1L = mesh_get_fieldpeaks_and_focality(m_ROIL,'field_idx','E_norm');
        summary2L = mesh_get_fieldpeaks_and_focality(m_ROIL,'field_idx','E_normal');
        
        FieldVal.EnormL(subj,montageNr).p95       = summary1L.perc_values(1);
        FieldVal.EnormL(subj,montageNr).p99       = summary1L.perc_values(2);
        FieldVal.EnormL(subj,montageNr).p999      = summary1L.perc_values(3);
        FieldVal.EnormL(subj,montageNr).peak_pos  = summary1L.XYZ_perc;
        
        FieldVal.EnormalL(subj,montageNr).p95     = summary2L.perc_values(1);
        FieldVal.EnormalL(subj,montageNr).p99     = summary2L.perc_values(2);
        FieldVal.EnormalL(subj,montageNr).p999    = summary2L.perc_values(3);
        try
        FieldVal.EnormalLNeg(subj,montageNr).p95  = summary2L.perc_neg_values(1);
        FieldVal.EnormalLNeg(subj,montageNr).p99  = summary2L.perc_neg_values(2);
        FieldVal.EnormalLNeg(subj,montageNr).p999 = summary2L.perc_neg_values(3);
        catch
        end
        summary1R = mesh_get_fieldpeaks_and_focality(m_ROIR,'field_idx','E_norm');
        summary2R = mesh_get_fieldpeaks_and_focality(m_ROIR,'field_idx','E_normal');
        
        FieldVal.EnormR(subj,montageNr).p95       = summary1R.perc_values(1);
        FieldVal.EnormR(subj,montageNr).p99       = summary1R.perc_values(2);
        FieldVal.EnormR(subj,montageNr).p999      = summary1R.perc_values(3);
        FieldVal.EnormR(subj,montageNr).peak_pos  = summary1R.XYZ_perc;
        
        FieldVal.EnormalR(subj,montageNr).p95     = summary2R.perc_values(1);
        FieldVal.EnormalR(subj,montageNr).p99     = summary2R.perc_values(2);
        FieldVal.EnormalR(subj,montageNr).p999    = summary2R.perc_values(3);
        try
        FieldVal.EnormalRNeg(subj,montageNr).p95  = summary2R.perc_neg_values(1);
        FieldVal.EnormalRNeg(subj,montageNr).p99  = summary2R.perc_neg_values(2);
        FieldVal.EnormalRNeg(subj,montageNr).p999 = summary2R.perc_neg_values(3);
        catch
        end
        %anatomical ROI
        peak_posAnatL = [-42.17 -20.37 7.96];
        peak_posAnatR = [46.24 -18.3 8.6];
        % % distance to peak position
        distAnatR = sqrt(sum(bsxfun(@minus,m.nodes,peak_posAnatR).^2,2));
        %
        % % extract nodes closer than 10 mm, and the related triangles and data
        node_idxAnatR = distAnatR<7;
        m_ROIAnatR    = mesh_extract_regions(m, 'node_idx', node_idxAnatR);
        
        distAnatL     = sqrt(sum(bsxfun(@minus,m.nodes,peak_posAnatL).^2,2));
        %
        % % extract nodes closer than 10 mm, and the related triangles and data
        node_idxAnatL = distAnatL<7;
        m_ROIAnatL    = mesh_extract_regions(m, 'node_idx', node_idxAnatL);
        
        % % show the extracted ROI:
        % % show the whole cortex semi-transparent
        mesh_show_surface(m,'showSurface',true,'facealpha',0.3);
        % % add the extracted area to the plot
        mesh_show_surface(m_ROIAnatR,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca);
        mesh_show_surface(m_ROIAnatL,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca);
        
        title(['10 mm spherical ROI around peak positions Anatomical' subjects{subj}]);
        %
        % % get some key results for the spherical ROI
        disp(' ')
        disp('10 mm spherical ROI around peak positions:')
        
        summary1AnatL = mesh_get_fieldpeaks_and_focality(m_ROIAnatL,'field_idx','E_norm');
        summary2AnatL = mesh_get_fieldpeaks_and_focality(m_ROIAnatL,'field_idx','E_normal');
        
        FieldVal.EnormAnatL(subj,montageNr).p95       = summary1AnatL.perc_values(1);
        FieldVal.EnormAnatL(subj,montageNr).p99       = summary1AnatL.perc_values(2);
        FieldVal.EnormAnatL(subj,montageNr).p999      = summary1AnatL.perc_values(3);
        FieldVal.EnormAnatL(subj,montageNr).peak_pos  = summary1AnatL.XYZ_perc;
        
        FieldVal.EnormalAnatL(subj,montageNr).p95     = summary2AnatL.perc_values(1);
        FieldVal.EnormalAnatL(subj,montageNr).p99     = summary2AnatL.perc_values(2);
        FieldVal.EnormalAnatL(subj,montageNr).p999    = summary2AnatL.perc_values(3);
        try
        FieldVal.EnormalAnatLNeg(subj,montageNr).p95  = summary2AnatL.perc_neg_values(1);
        FieldVal.EnormalAnatLNeg(subj,montageNr).p99  = summary2AnatL.perc_neg_values(2);
        FieldVal.EnormalAnatLNeg(subj,montageNr).p999 = summary2AnatL.perc_neg_values(3);
        catch
        end
        summary1AnatR = mesh_get_fieldpeaks_and_focality(m_ROIAnatR,'field_idx','E_norm');
        summary2AnatR = mesh_get_fieldpeaks_and_focality(m_ROIAnatR,'field_idx','E_normal');
        
        FieldVal.EnormAnatR(subj,montageNr).p95       = summary1AnatR.perc_values(1);
        FieldVal.EnormAnatR(subj,montageNr).p99       = summary1AnatR.perc_values(2);
        FieldVal.EnormAnatR(subj,montageNr).p999      = summary1AnatR.perc_values(3);
        FieldVal.EnormAnatR(subj,montageNr).peak_pos  = summary1AnatR.XYZ_perc;
        
        FieldVal.EnormalAnatR(subj,montageNr).p95     = summary2AnatR.perc_values(1);
        FieldVal.EnormalAnatR(subj,montageNr).p99     = summary2AnatR.perc_values(2);
        FieldVal.EnormalAnatR(subj,montageNr).p999    = summary2AnatR.perc_values(3);
        try
        FieldVal.EnormalAnatRNeg(subj,montageNr).p95  = summary2AnatR.perc_neg_values(1);
        FieldVal.EnormalAnatRNeg(subj,montageNr).p99  = summary2AnatR.perc_neg_values(2);
        FieldVal.EnormalAnatRNeg(subj,montageNr).p999 = summary2AnatR.perc_neg_values(3);
        catch
        end
        %get info from whole right and left hemispheres
         % % extract nodes closer than 10 mm, and the related triangles and data
        node_idxR = m.nodes(:,1)>0;
        m_ROIR    = mesh_extract_regions(m, 'node_idx', node_idxR);
        
        node_idxL = m.nodes(:,1)<0;
        m_ROIL    = mesh_extract_regions(m, 'node_idx', node_idxL);
        
        % % show the extracted ROI:
        % % show the whole cortex semi-transparent
        mesh_show_surface(m,'showSurface',true,'facealpha',0.3);
        % % add the extracted area to the plot
        mesh_show_surface(m_ROIR,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca);
        mesh_show_surface(m_ROIL,'showSurface',true,'surfaceColor',[0 0 1],'haxis',gca);
        
        title(['Individual Hemispheres' subjects{subj}]);
        %
        % % get some key results for the spherical ROI
        disp(' ')
        disp('Individual Hemispheres:')
        
        summary1L = mesh_get_fieldpeaks_and_focality(m_ROIL,'field_idx','E_norm');
        summary2L = mesh_get_fieldpeaks_and_focality(m_ROIL,'field_idx','E_normal');
        
        FieldVal.EnormLH(subj,montageNr).p95       = summary1L.perc_values(1);
        FieldVal.EnormLH(subj,montageNr).p99       = summary1L.perc_values(2);
        FieldVal.EnormLH(subj,montageNr).p999      = summary1L.perc_values(3);
        FieldVal.EnormLH(subj,montageNr).peak_pos  = summary1L.XYZ_perc;
        
        FieldVal.EnormalLH(subj,montageNr).p95     = summary2L.perc_values(1);
        FieldVal.EnormalLH(subj,montageNr).p99     = summary2L.perc_values(2);
        FieldVal.EnormalLH(subj,montageNr).p999    = summary2L.perc_values(3);
        try
        FieldVal.EnormalLHNeg(subj,montageNr).p95  = summary2L.perc_neg_values(1);
        FieldVal.EnormalLHNeg(subj,montageNr).p99  = summary2L.perc_neg_values(2);
        FieldVal.EnormalLHNeg(subj,montageNr).p999 = summary2L.perc_neg_values(3);
        catch
        end
        summary1R = mesh_get_fieldpeaks_and_focality(m_ROIR,'field_idx','E_norm');
        summary2R = mesh_get_fieldpeaks_and_focality(m_ROIR,'field_idx','E_normal');
        
        FieldVal.EnormRH(subj,montageNr).p95       = summary1R.perc_values(1);
        FieldVal.EnormRH(subj,montageNr).p99       = summary1R.perc_values(2);
        FieldVal.EnormRH(subj,montageNr).p999      = summary1R.perc_values(3);
        FieldVal.EnormRH(subj,montageNr).peak_pos  = summary1R.XYZ_perc;
        
        FieldVal.EnormalRH(subj,montageNr).p95     = summary2R.perc_values(1);
        FieldVal.EnormalRH(subj,montageNr).p99     = summary2R.perc_values(2);
        FieldVal.EnormalRH(subj,montageNr).p999    = summary2R.perc_values(3);
        try
        FieldVal.EnormalRHNeg(subj,montageNr).p95  = summary2R.perc_neg_values(1);
        FieldVal.EnormalRHNeg(subj,montageNr).p99  = summary2R.perc_neg_values(2);
        FieldVal.EnormalRHNeg(subj,montageNr).p999 = summary2R.perc_neg_values(3);
        catch
        end
        %SAVING ALL FIGURES
        for f = 1:8
            cName=[MainFolder '/GroupAna_0071_V2/' subjects{subj} '_fig_' num2str(f)];
            savefig(figure(f),cName)
        end
    end
end

save([MainFolder '/GroupAna/GroupEfieldAna_orig_wholeHem_0071_V2'], 'FieldVal');
