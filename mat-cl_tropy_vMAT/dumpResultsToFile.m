function dumpResultsToFile (filename, proof_offset, wsizes, sigmas, offset_hst_XCORR, offset_pks_XCORR, offset_hst_PY, offset_pks_PY)

  afile = fopen (filename, 'w');
    
    fprintf (afile, '\n********************* offset correntropy HISTOGRAM *********************\n');
    
    nsigs = length(sigmas);
    nws = length(offset_hst_PY{1}(1,:));
    
    for ws=1:nws
        fprintf (afile, '\n\nwindow size: %d\n\n', wsizes{ws});
        
        for i=1:nsigs
            fprintf (afile, '\t %.3f', sigmas{i});
        end
        fprintf (afile, '\t proof \n');

        ncomb = length(offset_hst_PY);

        for c=1:ncomb
            if (ws > length(offset_hst_PY{c}(1,:)))
                break;
            end
            fprintf (afile, 'c. %d\t', c);
            for i=1:nsigs
                fprintf (afile, '%d \t ', offset_hst_PY{c}{i,ws});
            end
            fprintf (afile, '%d\n', proof_offset(c));
        end
    end
    
    fprintf (afile, '\n********************* offset correntropy PEAKS *********************\n');
    
    for ws=1:nws
        fprintf (afile, '\n\nwindow size: %d\n\n', wsizes{ws});
        
        for i=1:nsigs
            fprintf (afile, '\t %.3f', sigmas{i});
        end
        fprintf (afile, '\t proof \n');

        ncomb = length(offset_pks_PY);

        for c=1:ncomb
            if (ws > length(offset_pks_PY{c}(1,:)))
                break;
            end
            fprintf (afile, 'c. %d\t', c);
            for i=1:nsigs
                fprintf (afile, '%d \t ', offset_pks_PY{c}{i,ws});
            end
            fprintf (afile, '%d\n', proof_offset(c));
        end
    end
%     
%    fprintf (afile, '\n********************* offset correntropy KL *********************\n');
%     
%     for ws=1:nws
%         fprintf (afile, '\n\nwindow size: %d\n\n', wsizes{ws});
%         
%         for i=1:nsigs
%             fprintf (afile, '\t %.3f', sigmas{i});
%         end
%         fprintf (afile, '\t proof \n');
% 
%         ncomb = length(offset_KL);
% 
%         for c=1:ncomb
%             if (ws > length(offset_KL{c}(1,:)))
%                 break;
%             end
%             fprintf (afile, 'c. %d\t', c);
%             for i=1:nsigs
%                 fprintf (afile, '%d \t ', offset_KL{c}{i,ws});
%             end
%             fprintf (afile, '%d\n', proof_offset(c));
%         end
%     end
    
    fprintf (afile, '\n\n********************* offset xcorr HISTOGRAM *********************\n\n');
    
    for ws=1:nws
            fprintf (afile, '\t %.4d', wsizes{ws});
    end
    
    fprintf (afile, '\t proof \n');
    

    ncomb = length(offset_hst_XCORR);
    
    for c=1:ncomb
        if (ws > length(offset_hst_XCORR{c}(:)))
            break;
        end
        fprintf (afile, 'c. %d \t', c);
        for ws=1:nws
            fprintf (afile, '%.2f \t', offset_hst_XCORR{c}{ws});
        end
        fprintf (afile, '%.2f\n', proof_offset(c));
    end
    
    fprintf (afile, '\n\n********************* offset xcorr PEAKS *********************\n\n');
    
    for ws=1:nws
            fprintf (afile, '\t %.4d', wsizes{ws});
    end
    
    fprintf (afile, '\t proof \n');
    
    ncomb = length(offset_pks_XCORR);
    
    for c=1:ncomb
        if (ws > length(offset_pks_XCORR{c}(:)))
            break;
        end
        fprintf (afile, 'c. %d \t', c);
        for ws=1:nws
            fprintf (afile, '%.2f \t', offset_pks_XCORR{c}{ws});
        end
        fprintf (afile, '%.2f\n', proof_offset(c));
    end
    
    fprintf (afile, '\n\nCombinations offsets: \n\n');
    
    ncomb = length(proof_offset);
    
    for cmb=1:ncomb
        fprintf (afile, 'c. %d - %f\n', cmb, proof_offset(cmb));
    end
    
    fclose(afile);


end