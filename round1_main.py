# coding: utf-8
import analysis_fct as fct
import find_orf as orf
import sys

print 'Argument List:', str(sys.argv)

metaG_name = sys.argv[1]
input_way = '/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/Programme/data/'+metaG_name+'/'
output_way = 'output/'+metaG_name+'/'

gff_file = input_way+'ICM0104'+metaG_name+'.a.gff'
scaffold_file = input_way+"/ICM0104"+metaG_name+".a.fasta"
fna_file = input_way+'ICM0104'+metaG_name+'.a.fna'
#table_hmm ='data/table_per_domain_hmmsearh_eggpfam_vs_MP0310_36570.txt'
table_hmm =input_way+'table_per_domain_hmmsearh_eggpfam_vs_MP0310_36570.txt'

scaffold = "ICM0104MP0310_1000105"

scaffold_list = sorted(fct.get_list_scaffold(table_hmm))


#scaffold_list = [scaffold]
##Output files:
out_faa_file = output_way+metaG_name+'TESTorf_adjacent.faa'
out_gff_file = output_way+metaG_name+'TESTorf_adjacent.gff'

###PAS TRES ELEGANT MAIS PERMET DE SUPPRIMER LE CONTENU DU FICHIER :-/ 
for fname in [out_faa_file,out_gff_file]:
    fl = open(fname, 'w')
    fl.close()

#outpout = 'table_group_lenva.txt'
##reste the file 
#fl = open(outpout,'w')
#sep = '\t'
#header='n°group'+sep+'contig|gene'+sep*3+'start'+sep+'end'+sep+'strand'+sep+'domains\n'
#fl.write(header)
#fl.close()




#Parameter for the orfinder        
threshold = 80 # min(taille)  -10 
id_genTable = 11

##Open of the scaffold file ! Pour le fast fasta function de ORF
fl_sca = open(scaffold_file, 'r')
line = fl_sca.readline()

##Open ICM file fna pour recupérer les sequences fasta des gene prédit en nt
dico_seq = {}
fl_fna = open(fna_file, 'r')
line = fl_fna.readline()
dico_seq["fl"]= fl_fna
dico_seq["line"] = line
dico_seq["codon_start"] = orf.getGeneticCode(id_genTable)['start']

##counter
cnt_lonely_gene = 0
cnt_adj_orf = 0

composition = {'nonva':{}, 'va':{}, 'mixt':{}, 'lonely_gene':0}


for scaffold in scaffold_list[:10]:
    #print scaffold
    liste = fct.list_line(gff_file, scaffold) #list of the gff line
    genes = fct.get_hmm_genes(scaffold, table_hmm, liste) #list of obj gene
    genes.check_size(dico_seq)
    
    genes.get_adj() #set the adj gene for each gene which has 

    genes.group(composition)
    genes.get_lonely_gene()

    ORF, fl, line = orf.get_list_orf(scaffold_file, scaffold, threshold, id_genTable, fl_sca, line)
    
    #faire matché les genes seules avec les orf 
    for g in genes.lonely_gene:
        ORF.find_the_orf(g)
        
    ##for o in ORF.list_gene_hmm:
        ##print o.gff_obj.gene_number
        
    ORF.get_adj_orf()
    
    index_gene = fct.get_gff_obj(liste[-1]).gene_number#index pour pour noter les 
    ORF.write_faa_gff(index_gene, out_gff_file, out_faa_file)
    #print len(ORF.set_adj_orf)    
    ##counter
    cnt_lonely_gene += len(ORF.list_gene_hmm)
    cnt_adj_orf += len(ORF.set_adj_orf)
    #for g in genes.:
        #orf.find_the_orf(dico_orf_list[g.strand], corres[g.strand], g, true_genes_list)
    
    #print listORF
print 'Round 1 :'
print 'group of gene valide :',composition['va']
print 'number of adjacent orf to give to hmm :', cnt_adj_orf
print 'number of lonely hmm gene :',  cnt_lonely_gene

fl_sca.close()





