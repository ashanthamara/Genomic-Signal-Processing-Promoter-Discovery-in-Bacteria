clc;
clear;

/*Priliminaries*/
// Import support functions to this file
exec('support_functions.sce');
// Define raw data files
fasta_file = 'NZ_CP075108.1.fasta';
protein_table = 'NZ_CP075108.1.txt';

//fasta_file = 'NZ_CP075108.1_thamara.fasta';
//protein_table = 'NZ_CP075108_thamara.1.txt';

//fasta_file = 'sequence_dineth.fasta';
//protein_table = 'proteins_152_1556684_dineth.txt';
// define threshold parameters
upstream_thresh = 50;
downstream_tresh = 3;

// Read the protein table and exact the gene locations of the sense and the anti sense strands
[sen_gene_pos, antisen_gene_pos, sen_ncod, antisen_ncod] = get_protein_pos_array(protein_table);

sen_gen_count = size(sen_gene_pos, 1);
antisen_gen_count = size(antisen_gene_pos, 1);

printf('\n------- Number of genes from the protein table -------');
printf('\nNumber of Genes in the Sense Strand : %d', sen_gen_count); 
printf('\nNumber of Genes in the Anti-Sense Strand: %d\n', antisen_gen_count);

// Question 1 : Remove the genes which are less than 50 bases closer to the other gene.

/******************** SENSE STRAND**********************/
selected_sen_genes = [];

// Methoyonin ascii values in sense strand
ATG = ascii('ATG');
GTG = ascii('GTG');

// Getting the first gene of the sense strand
codon_len = 3; 
start_codon = get_fasta_at(fasta_file, sen_gene_pos(1,1),sen_gene_pos(1,1) + codon_len, 1);
if ((start_codon == ATG)||(start_codon == GTG)) then
    selected_sen_genes(1,:) = sen_gene_pos(1,:);
end

// Selecting the other suitable genes from the sense strand
for gene_index = 1:1:sen_gen_count-1    
    
    next_gene_index = gene_index + 1;
//    disp(gene_index);
    if (sen_gene_pos(next_gene_index, 1) - sen_gene_pos(gene_index, 2)) > upstream_thresh then
        
        start_codon = get_fasta_at(fasta_file, sen_gene_pos(next_gene_index,1), sen_gene_pos(next_gene_index,1) + 3,1);
        
        if ((start_codon == ATG)||(start_codon == GTG)) then
            selected_sen_genes = [selected_sen_genes; sen_gene_pos(next_gene_index,:)];
        end
    end
end

selected_sen_gene_count = size(selected_sen_genes, 1);

/******************** ANTI SENSE STRAND**********************/
selected_antisen_genes = [];

// Methoyonin ascii values in anti sense strand
GTA = ascii('GTA');
GTG = ascii('GTG');

// Getting the first gene of the anti sense strand
start_codon = get_fasta_at(fasta_file, antisen_gene_pos(1,2)- 2, antisen_gene_pos(1,2) + 1, 0);
if ((start_codon == GTA)||(start_codon == GTG)) then
    selected_antisen_genes(1,:) = antisen_gene_pos(1,:);
end

// Selecting the other suitable genes from the anti sense strand
gene_index = 0;
for gene_index = 1:1:antisen_gen_count-1    
    
    next_gene_index = gene_index + 1;
//    disp(gene_index);
    if (antisen_gene_pos(next_gene_index,1) - antisen_gene_pos(gene_index,2)) > upstream_thresh then
     
        start_codon = get_fasta_at(fasta_file, antisen_gene_pos(next_gene_index, 2)-2, antisen_gene_pos(next_gene_index, 2) + 1, 0);
        
        if ((start_codon == GTA)||(start_codon == GTG)) then
            selected_antisen_genes = [selected_antisen_genes; antisen_gene_pos(next_gene_index,:)];
        end
    end
end

selected_antisen_gene_count = size(selected_antisen_genes, 1);

printf('\n------- Number of Valid genes -------');
printf('\nNumber of Valid Genes in the Sense Strand : %d', selected_sen_gene_count); 
printf('\nNumber of Valid Genes in the Anti-Sense Strand: %d\n', selected_antisen_gene_count);

savematfile('selected_sen_genes.mat', 'selected_sen_genes', '-v7.3');
savematfile('selected_antisen_genes.mat', 'selected_antisen_genes','-v7.3');
// Base distribution of bases in non-neglected sense strand and antisense strands
base_dis_sen = zeros(4,1);
base_dis_antisen = zeros(4,1);
// for the sense strand
for count_index = 1:1:selected_sen_gene_count
    gene = get_fasta_at(fasta_file, selected_sen_genes(count_index,1), selected_sen_genes(count_index,2), 1);
    
    for base_index = 1:1:size(gene, 2)
        base_dis_sen = base_inc(gene(base_index), base_dis_sen);
    end
end
printf('\n------- Base Distribution | Sense Strand -------');
printf('\nNumber of Base A: %d', base_dis_sen(1,1)); 
printf('\nNumber of Base C: %d', base_dis_sen(2,1)); 
printf('\nNumber of Base G: %d', base_dis_sen(3,1)); 
printf('\nNumber of Base T: %d\n', base_dis_sen(4,1)); 

// for the anti sense strand
for count_index = 1:1:selected_antisen_gene_count
    gene = get_fasta_at(fasta_file, selected_antisen_genes(count_index,1), selected_antisen_genes(count_index,2),0);
    
    for base_index = 1:1:size(gene,2)
        base_dis_antisen = base_inc(gene(base_index), base_dis_antisen);
    end
end
printf('\n------- Base Distribution | Anti-Sense Strand -------');
printf('\nNumber of Base A: %d', base_dis_antisen(1,1)); 
printf('\nNumber of Base C: %d', base_dis_antisen(2,1)); 
printf('\nNumber of Base G: %d', base_dis_antisen(3,1)); 
printf('\nNumber of Base T: %d\n', base_dis_antisen(4,1)); 


/********************** Question 2 *************************/

trainingset_gene_count = 1000; // As given in the question

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
for gene_no = 1:trainingset_gene_count 
    
    seq_start_pos = selected_sen_genes(gene_no,1) - prom_begin;
    seq_end_pos = selected_sen_genes(gene_no,1) - prom_end;
    
    prom_seq = get_fasta_at(fasta_file, seq_start_pos, seq_end_pos, 1); // get the bases in the sequence
    
    // Intact query >> match = 1, mismmatch = -1
    [align_x, align_y, prom_pos] = traceback_prom(prom_seq, pribnow_box, 1, -1, gap_penalty); // W matching
    
    // check whether the pribnow box is aligned in the range of the promoter sequence to get the promoter sequence
    if ((prom_pos > 0) & (prom_len - ppm_length - prom_pos > 0)) then
            
            // Update the pribnow position frequency matrix
            promotor = prom_seq((prom_pos + 1): (prom_pos + ppm_length));

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

entropy_thresh = 0.04; // by observing the plot

savematfile('PPMatrix.mat', 'PPMatrix', '-v7.3');
save('entropy_thresh.dat','entropy_thresh');

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
    testset_sen_start_loc = trainingset_gene_count + 1;
    
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
    sen_testset_gene_count = selected_sen_gene_count - trainingset_gene_count;
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
        
        flipped_prom_seq = flipdim(prom_seq, 2); // flip the sequence since the strand is anti sense strand
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





