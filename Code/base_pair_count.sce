clc;
clear;

/*Priliminaries*/
// Import support functions to this file
exec('support_functions.sce');
// Define raw data files
fasta_file = 'NZ_CP075108.1.fasta';
protein_table = 'NZ_CP075108.1.txt';


[sen_gene_pos, antisen_gene_pos, sen_ncod, antisen_ncod] = get_protein_pos_array(protein_table);

sen_gen_count = size(sen_gene_pos, 1);
antisen_gen_count = size(antisen_gene_pos, 1);

sen_ncod_count = size(sen_ncod, 1);

printf('\n------- Number of genes from the protein table -------');
printf('\nNumber of Genes in the Sense Strand : %d', sen_gen_count); 
printf('\nNumber of Genes in the Anti-Sense Strand: %d\n', antisen_gen_count);

// Base distribution of bases in non-neglected sense strand and antisense strands
base_dis_sen_cod = zeros(4,1);
base_dis_antisen_cod = zeros(4,1);
// for the sense strand
for count_index = 1:1:sen_gen_count
    gene = get_fasta_at(fasta_file, sen_gene_pos(count_index,1), sen_gene_pos(count_index,2), 1);
    
    for base_index = 1:1:size(gene, 2)
        base_dis_sen_cod = base_inc(gene(base_index), base_dis_sen_cod);
    end
end
printf('\n------- Base Distribution | Sense Strand -------');
printf('\nNumber of Base A: %d', base_dis_sen_cod(1,1)); 
printf('\nNumber of Base C: %d', base_dis_sen_cod(2,1)); 
printf('\nNumber of Base G: %d', base_dis_sen_cod(3,1)); 
printf('\nNumber of Base T: %d\n', base_dis_sen_cod(4,1)); 

// non coding base calculation
base_dis_sen_noncod = zeros(4,1);
// for the sense strand
for count_index = 1:1:sen_gen_count
    ncod_region = get_fasta_at(fasta_file, sen_ncod(count_index,1), sen_ncod(count_index,2), 1);
    
    for base_index = 1:1:size(ncod_region, 2)
        base_dis_sen_noncod = base_inc(ncod_region(base_index), base_dis_sen_noncod);
    end
end
printf('\n------- Base Distribution | Sense Strand non coding region -------');
printf('\nNumber of Base A: %d', base_dis_sen_noncod(1,1)); 
printf('\nNumber of Base C: %d', base_dis_sen_noncod(2,1)); 
printf('\nNumber of Base G: %d', base_dis_sen_noncod(3,1)); 
printf('\nNumber of Base T: %d\n', base_dis_sen_noncod(4,1)); 

total_A = base_dis_sen_cod(1,1) + base_dis_sen_noncod(1,1);
total_C = base_dis_sen_cod(2,1) + base_dis_sen_noncod(2,1);
total_G = base_dis_sen_cod(3,1) + base_dis_sen_noncod(3,1);
total_T = base_dis_sen_cod(4,1) + base_dis_sen_noncod(4,1);

total_bases = total_A + total_C + total_G + total_T;

printf('\n------- Base Distribution | Sense Strand non coding region -------');
printf('\nTotal Number of Base A: %d', total_A); 
printf('\nTotal nNumber of Base C: %d', total_C); 
printf('\nTotal Number of Base G: %d', total_G); 
printf('\nTotal Number of Base T: %d\n', total_T);

printf('\nTotal Number of Base pairs: %d\n', total_bases); 
// for the anti sense strand
for count_index = 1:1:antisen_gen_count
    gene = get_fasta_at(fasta_file, antisen_gene_pos(count_index,1), antisen_gene_pos(count_index,2),0);
    
    for base_index = 1:1:size(gene,2)
        base_dis_antisen_cod = base_inc(gene(base_index), base_dis_antisen_cod);
    end
end
printf('\n------- Base Distribution | Anti-Sense Strand -------');
printf('\nNumber of Base A: %d', base_dis_antisen_cod(1,1)); 
printf('\nNumber of Base C: %d', base_dis_antisen_cod(2,1)); 
printf('\nNumber of Base G: %d', base_dis_antisen_cod(3,1)); 
printf('\nNumber of Base T: %d\n', base_dis_antisen_cod(4,1)); 
