clc;
clear;
// Import support functions to this file
exec('support_functions.sce');

loadmatfile('selected_sen_genes.mat');
loadmatfile('selected_antisen_genes.mat');
fasta_file = 'NZ_CP075108.1.fasta';
selected_sen_gene_count = size(selected_sen_genes, 1);
selected_antisen_gene_count = size(selected_antisen_genes, 1);

/********************** Question 2 *************************/

trainingset_gene_count = 100; // 10 & 100

/*
Pribnow bow search DATA
Start position : 30
End Position : 5
*/

prom_begin = 30;
prom_end = 5;

prom_len = prom_begin - prom_end; 
ppm_length = 10; // length of the position probability matrix

initial_k = 0.01; // initial k value for the PPM
prib_freq_mat = ones(4, ppm_length).* initial_k; //creating initial frequency matrix
pribnow_box = ascii('TATAAT');

// Get the promotor sequences aligned with the pribnow box and calculate the frequency
for gene_no = 725:725+trainingset_gene_count 
    
    seq_start_pos = selected_sen_genes(gene_no,1) - prom_begin;
    seq_end_pos = selected_sen_genes(gene_no,1) - prom_end;
    
    prom_seq = get_fasta_at(fasta_file, seq_start_pos, seq_end_pos, 1); // get the bases in the sequence
    
    // Intact query >> match = 1, mismmatch = -1
    [align_x, align_y, prom_pos] = traceback_prom(prom_seq, pribnow_box, 1, -1, gap_penalty); // W matching
    
    // check whether the pribnow box is aligned in the range of the promoter sequence to get the promoter sequence
    if ((prom_pos > 0) & (prom_len - ppm_length - prom_pos > 0)) then
            
            // Update the pribnow position frequency matrix
            promotor = prom_seq((prom_pos + 1):(prom_pos + ppm_length));

            prib_freq_mat = update_pos_freq_mat(promotor, prib_freq_mat, ppm_length);
     end
end

printf('\n------- Position Frequency Matrix -------');
disp(prib_freq_mat); 

PPMatrix = get_ppm(prib_freq_mat); // get the position probability matrix

printf('\n\n------- Position Probabillity Matrix -------');
disp(PPMatrix); 

/********************** Question 3 *************************/

// Find the entrophy using the equiprobable baseline

baseline_mat = ones(1, 4)*0.25;

[w,entrophy] = ppm_info(PPMatrix, baseline_mat);

S = sum(entrophy, 1); // Get the total entropy for each position of the PPM
prom_pos = 1 : ppm_length; // Positions of PPM

//Plot the entropy distribution
plot(prom_pos, S);
xtitle('Entropy vs Number of positions in Promotor','Promoter Position','Entropy');

entropy_thresh = 0.08; // by observing the plot n=10 >> 0.2 , n=100 >> 0.08

/********************** Question 4 *************************/
// get the ppm_seq with maximum consencus score
ppm_seq = ones(1,ppm_length);
for clmn = 1: ppm_length
    max_base_val = max(PPMatrix(:,clmn));
    for i = 1:4
        if i == 1 then
            base = 'A';
        elseif i == 2 then
            base = 'C';
        elseif i == 3 then
            base = 'G';
        else
            base = 'T';
        end
        
        if (max_base_val <= PPMatrix(i,clmn)) then
            
            ppm_seq(1,clmn) = ascii(base);
        end     
    end
end

disp(ppm_seq);

// Calculate the Concensus Scores
bench_score = stat_align(ppm_seq, PPMatrix); // for the initial length
reduced_bench_score = stat_align_entropy(ppm_seq, PPMatrix, entropy_thresh);  // Without the lower positions than the threshold

printf('\n\n------- Benchmark Consensus Scores -------');
printf('\nInitial PPM Benchmark Consensus Score = %.5f\n' , bench_score);
printf('\nReduced PPM Benchmark Consensus Score = %.5f \n' , reduced_bench_score);

// test set data
//flipped_PPMatrix = flipdim(PPMatrix, 2);  //filp the ppm to align with the anti sense strand genes
for stat_align_thresh = -1:-1:-5
    // Sequence aligning for other genes in the sense strand
    testset_sen_valid_proms = 0;
    testset_sen_reduced_valid_proms = 0;
    testset_sen_start_loc = 1001; // using previous test set data
    
    for gene_index = testset_sen_start_loc:1:selected_sen_gene_count
        
        seq_start_pos = selected_sen_genes(gene_index,1) - prom_begin;
        seq_end_pos = selected_sen_genes(gene_index,1) - prom_end;
        
        prom_seq = get_fasta_at(fasta_file, seq_start_pos, seq_end_pos, 1); // get the bases in the sequence
        
        stat_align_score = stat_align(prom_seq, PPMatrix);
        reduced_stat_align_score = stat_align_entropy(prom_seq, PPMatrix, entropy_thresh);
        
        // check the score lies within the threshold
        max_align_score = max(stat_align_score);
        if ((max_align_score - bench_score) > stat_align_thresh) then
            testset_sen_valid_proms = testset_sen_valid_proms + 1;
        end
        
        max_reduced_align_score = max(reduced_stat_align_score);
        if ((max_reduced_align_score - reduced_bench_score) > stat_align_thresh) then
            testset_sen_reduced_valid_proms = testset_sen_reduced_valid_proms + 1;
        end    
    end
    sen_testset_gene_count = selected_sen_gene_count - 1000;
    printf('\n\n------- Threshold = %d -------', stat_align_thresh);
    printf('\nSense Strand | Test set gene count = %d', sen_testset_gene_count);
    printf('\nSense Strand | Test set VALID gene count = %d', testset_sen_valid_proms);
    printf('\nSense Strand | Test set VALID gene percentage = %.2f', 100*testset_sen_valid_proms/sen_testset_gene_count);
    printf('\nSense Strand | Test set Reduced VALID gene count = %d', testset_sen_reduced_valid_proms);
    printf('\nSense Strand | Test set Reduced VALID gene percentage = %.2f', 100*testset_sen_reduced_valid_proms/sen_testset_gene_count);

    // Sequence aligning for other genes in the anti sense strand
    testset_antisen_valid_proms = 0;
    testset_antisen_reduced_valid_proms = 0;
    testset_antisen_start_loc = 1;
    
    for gene_index = testset_antisen_start_loc:1:selected_antisen_gene_count
        
        seq_start_pos = selected_antisen_genes(gene_index,2) + prom_end;
        seq_end_pos = selected_antisen_genes(gene_index,2) + prom_begin;
        
        prom_seq = get_fasta_at(fasta_file, seq_start_pos, seq_end_pos, 0); // get the bases in the sequence
        
        flipped_prom_seq = flipdim(prom_seq, 2);
        stat_align_score = stat_align(flipped_prom_seq, PPMatrix);
        reduced_stat_align_score = stat_align_entropy(flipped_prom_seq, PPMatrix, entropy_thresh);
        
        // check the score lies within the threshold
        max_align_score = max(stat_align_score);
        if ((max_align_score - bench_score) > stat_align_thresh) then
            testset_antisen_valid_proms = testset_antisen_valid_proms + 1;
        end
        
        max_reduced_align_score = max(reduced_stat_align_score);
        if ((max_reduced_align_score - reduced_bench_score) > stat_align_thresh) then
            testset_antisen_reduced_valid_proms = testset_antisen_reduced_valid_proms + 1;
        end    
    end
    printf('\n\nAnti Sense Strand | Test set gene count = %d', selected_antisen_gene_count);
    printf('\nAnti Sense Strand | Test set VALID gene count = %d', testset_antisen_valid_proms);
    printf('\nAnti Sense Strand | Test set VALID gene percentage = %.2f ', 100*testset_antisen_valid_proms/selected_antisen_gene_count);
    printf('\nAnti Sense Strand | Test set Reduced VALID gene count = %d', testset_antisen_reduced_valid_proms);
    printf('\nAnti Sense Strand | Test set Reduced VALID gene percentage = %.2f', 100*testset_antisen_reduced_valid_proms/selected_antisen_gene_count);

    
    total_test_set_gene_count = sen_testset_gene_count + selected_antisen_gene_count;
    total_valid_gene_count = testset_sen_valid_proms + testset_antisen_valid_proms;
    total_valid_reduced_gene_count = testset_sen_reduced_valid_proms + testset_antisen_reduced_valid_proms;
    
    printf('\n\nTotal | Test set gene count = %d', total_test_set_gene_count);
    printf('\nTotal | Test set VALID gene count = %d', total_valid_gene_count);
    printf('\nTotal | Test set VALID gene percentage = %.2f', 100*total_valid_gene_count/total_test_set_gene_count);
    printf('\nTotal | Test set Reduced VALID gene count = %d', total_valid_reduced_gene_count);
    printf('\nTotal | Test set Reduced VALID gene percentage = %.2f', 100*total_valid_reduced_gene_count/total_test_set_gene_count);
    
    printf('\nTotal | Proportion of genes that do not have promotors = %.2f', 100*(total_test_set_gene_count-total_valid_gene_count)/total_test_set_gene_count);
    printf('\nTotal | Proportion of reduced genes that do not have promotors = %.2f', 100*(total_test_set_gene_count-total_valid_reduced_gene_count)/total_test_set_gene_count);    
end
