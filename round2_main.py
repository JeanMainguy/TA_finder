# coding: utf-8
import analysis_fct as fct
import find_orf as orf
import sys

print 'Argument List:', str(sys.argv)

metaG_name = sys.argv[1]
input_way = '/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/Programme/data/'+metaG_name+'/'
output_way = 'output/'+metaG_name+'/'
gff_file = input_way+'ICM0104'+metaG_name+'.a.gff'
gff_file_orf = output_way+'/'+metaG_name+'_orf_adjacent.gff'

scaffold_file = input_way+"/ICM0104"+metaG_name+".a.fasta"
fna_file = input_way+'ICM0104'+metaG_name+'.a.fna'
#table_hmm ='data/table_per_domain_hmmsearh_eggpfam_vs_MP0310_36570.txt'
table_hmm =output_way+'merged_table_hmm_'+metaG_name+'.txt'



scaffold = "ICM0104MP0310_1000266"

scaffold_list = sorted(fct.get_list_scaffold(table_hmm))
index = orf.get_index_sca(scaffold_file, scaffold_list)

#scaffold_list = [scaffold]

##Output files:
rasta_group_file = output_way+metaG_name+'group_rasta.txt'
gene_group_file = output_way+metaG_name+'group_gene.txt'
group_rasta_mixt = output_way+metaG_name+'group_rasta_mixt.txt'
csvfile = output_way+metaG_name+'rasta_domain_count.csv'


###PAS TRES ELEGANT MAIS PERMET DE SUPPRIMER LE CONTENU DU FICHIER :-/ 
sep = '\t'
header='n°group'+sep+'contig|gene'+sep*3+'start'+sep+'end'+sep+'strand'+sep+'domains\n'

for fname in [rasta_group_file, gene_group_file, group_rasta_mixt]:
    fl = open(fname, 'w')
    fl.write(header)
    fl.close()


#Open ICM file fna pour recupérer les sequences fasta des gene prédit en nt
id_genTable = 11
dico_seq = {}
fl_fna = open(fna_file, 'r')
line = fl_fna.readline()
dico_seq["fl"]= fl_fna
dico_seq["line"] = line
dico_seq["codon_start"] = orf.getGeneticCode(id_genTable)['start']

#Rasta domain : 
rasta_domain=["NOG.COG1396.clustalo_raw" , "NOG.COG1487.meta_raw" , "NOG.COG1598.clustalo_raw" , "NOG.COG1724.meta_raw" , "NOG.COG1848.meta_raw" , "NOG.COG2002.meta_raw" , "NOG.COG2026.meta_raw" , "NOG.COG2161.meta_raw" , "NOG.COG2336.meta_raw" , "NOG.COG2337.meta_raw" , "NOG.COG3093.meta_raw" , "NOG.COG3549.meta_raw" , "NOG.COG3550.meta_raw" , "NOG.COG3609.meta_raw" , "NOG.COG3654.meta_raw" , "NOG.COG3668.meta_raw" , "NOG.COG4113.meta_raw" , "NOG.COG4118.meta_raw" , "NOG.COG4226.meta_raw" , "NOG.COG4423.meta_raw" , "NOG.COG4456.meta_raw" , "NOG.COG4691.meta_raw" , "NOG.COG5302.meta_raw" , "NOG.COG5499.meta_raw" , "PF01381.20" , "PF01402.19" , "PF01845.15" , "PF01850.19" , "PF02452.15" , "PF04014.16" , "PF04221.10" , "PF05015.11" , "PF05016.13" , "PF05534.10" , "PF07362.10"]

rasta_dico = {}
for d in rasta_domain:
    rasta_dico[d] = {'inside_rasta_group':0, 'lonely_rasta':0, 'inside_group_mixt':0}

composition = {'nonva':{}, 'va':{}, 'mixt':{},'rasta':{},'rasta_mixt':{}, 'lonely_gene':0}

for scaffold in scaffold_list[:10]:
    #print scaffold
    liste = fct.list_line(gff_file, scaffold) #list of the gff line
    liste += fct.list_line(gff_file_orf, scaffold)
    #for l in liste:
        #print l
    genes = fct.get_hmm_genes(scaffold, table_hmm, liste) #list of obj gene

    
    
    genes.check_size(dico_seq)
    #for gene in genes.liste:
        #if True:
            #print gene.gff_line
    genes.set_best_rasta_domain(rasta_domain)
    genes.get_adj() #set the adj gene for each gene which has 
    genes.get_lonely_gene()
    genes.group(composition)
    
    
    
    genes.write_to_file(rasta_group_file)
    genes.write_to_file(gene_group_file, 'group_len_va')
    genes.write_to_file(group_rasta_mixt, 'group_rasta_mixt')
    genes.get_info(composition)
    genes.count_rasta_domain(rasta_dico)

    #ORF, fl, line = orf.get_list_orf(scaffold_file, scaffold, threshold, id_genTable,index, fl_sca, line)
    
    #faire matché les genes seules avec les orf 
    #for g in genes.lonely_gene:
        #ORF.find_the_orf(g)
        
    ##for o in ORF.list_gene_hmm:
        ##print o.gff_obj.gene_number
        
    #ORF.get_adj_orf()
    
    #index_gene = fct.get_gff_obj(liste[-1]).gene_number#index pour pour noter les 
    #ORF.write_faa_gff(index_gene, out_gff_file, out_faa_file)
    #print len(ORF.set_adj_orf)    
    
    #for g in genes.:
        #orf.find_the_orf(dico_orf_list[g.strand], corres[g.strand], g, true_genes_list)
    
    #print listORF

print 'Group with a valid length gene :'
for k in composition['va'].keys():
    print composition['va'][k], 'groups made of', k, 'genes'
print 'Group mixt with at least one rasta domain gene :'
for k in composition['rasta_mixt'].keys():
    print composition['rasta_mixt'][k], 'groups made of', k, 'genes'
print 'Group with rasta gene only :'
for k in composition['rasta'].keys():
    print composition['rasta'][k], 'groups made of', k, 'genes'
print 'number of lonely gene with a rasta domain :', composition['lonely_rasta']
print 'number of lonely gene identified with a domain:', composition['lonely_gene']

fct.write_csv_domain(csvfile, rasta_dico)



