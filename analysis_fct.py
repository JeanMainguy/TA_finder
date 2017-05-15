# coding: utf-8
import re
from operator import attrgetter
import find_orf as orf

def get_hmm_genes(scaffold, table_hmm, gff_lines):        
    """
    Complete the list_gff_obj with info from the table hmm and 
    """
    
    fl = open(table_hmm, 'r')
    gene_index = set() #to get adjacent gene later. set() --> to get only once the index per gene
    set_TAgene = set() #list of genes with domain found by hmm
    
    genes = ListGene()

    for l in fl:
        if l[:len(scaffold)] == scaffold:
            
            nb_gene, domain = hmmtable_parser(l)

            #print nb_gene
            if nb_gene in genes.numbers:
                genes.appendomain(domain, nb_gene)
                
            else:
                gene = get_gff_obj(gff_lines[nb_gene-1])
                gene.domain.append(domain)
                genes.listappend(gene)
    #print len(genes.liste)
    return genes
     
            
            #number = dico_domain["gene_number"]
            #gene_index.add(number-1) # -1 is there because it is an index of gene number starting with 1 and not 0
            
            #gene = list_gff_obj[number-1]
            
            #gene.domain.append(dico_domain)
            #set_TAgene.add(gene)

    #get_adjacent(list_gff_obj, list(gene_index))
    
    #list_TAgene = list(set_TAgene)
    
    #for g in list_TAgene:
        #g.TA_valid()
    
    #list_of_group, list_gene_solo = group(list_TAgene)

    #return list_of_group, list_gene_solo

def hmmtable_parser(line):
    """
    Regex : parse a line of the hmmtable and return a dico containing the info taken from the line
    """
    pattern = re.compile(r"""
    (?P<scaffold>[^|]+)\|(?P<gene_number>\d{1,4}) #scaffold name and gene_number separated by a |
    \s+-\s+
    (?P<gene_len>\d+) #length of the gene sequence in residu
    \s+
    (?P<domain>[\w.-]+) #domain name
    \s+
    (?P<domain_acc>[\w.+-]+) #domain acc pour pfam
    \s+\d+\s+
    (?P<evalue>[\w.-]+) # E-value of the overall sequence/profile comparison (including all domains).
    \s+
    (?P<score>[\d.]+) # Bit score of the overall sequence/profile comparison (
    \s+[\d.]+\s+
    (?P<domain_number>\d+) # This domain’s number
    \s+
    (?P<total_domain>\d+) # The total number of domains reported in the sequence, ndom.
    \s+[\w.-]+\s+[\w.+-]+\s+[\w.-]+\s+[\w.-]+\s+\d+\s+\d+\s+
    (?P<ali_from>\d+) # The start of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
    \s+
    (?P<ali_to>\d+) # The end of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
    \s+
    (?P<env_from>\d+) # The start of the domain envelope on the sequence, numbered 1..L for a sequence of L residues.  
                      # The envelope defines a subsequence for which their is substantial probability mass supporting a homologous domain, 
                      # whether or not a single discrete alignment can be identified. 
                      # The envelope may extend beyond the endpoints of the MEA alignment, and in fact often does, for weakly scoring domains.
    \s+
    (?P<env_to>\d+) # The end of the domain envelope on the sequence, numbered 1..L for a sequence of L residues.
    \s+[\w.-]+\s+[\w.-]+
    """, re.VERBOSE)
    """
    regex in one line: 
    (?P<scaffold>[^|]+)\|(?P<gene_number>\d{1,4})\s+-\s+(?P<gene_len>\d+)\s+(?P<domain>[\w.-]+)\s+(?P<domain_acc>[\w.+-]+)\s+\d+\s+(?P<evalue>[\w.-]+)\s+(?P<score>[\d.]+)\s+[\d.]+\s+(?P<domain_number>\d+)\s+(?P<total_domain>\d+)\s+[\w.+-]+\s+([\w.+-]+)\s+[\w.+-]+\s+[\w.+-]+\s+\d+\s+\d+\s+(?P<ali_from>\d+)\s+(?P<ali_to>\d+)\s+(?P<env_from>\d+)\s+(?P<env_to>\d+)\s+[\w.-]+\s+[\w.-]+
    
    link : https://regex101.com/r/kKHzsB/1
    """
    
    result = pattern.match(line)
    
    domain = Domain(domain_number=result.group("domain_number"), ali_from=result.group("ali_from"), ali_to=result.group("ali_to"), env_from=result.group("env_from"), env_to=result.group("env_to"), domain_name=result.group("domain"), domain_acc=result.group("domain_acc"), e_value=result.group("evalue"), score=result.group("score"), line=line)
    gene_number = int(result.group("gene_number"))
    return gene_number, domain

def list_line(gff_file, scaffold):
    fl = open(gff_file, 'r')
    list_line_gff = []
    l = fl.readline()
    
    while l[:len(scaffold)] != scaffold and l: #gff file is sorted by scaffold so we use this feature in order to ot scan all the file  
        l = fl.readline()
        
    while l[:len(scaffold)] == scaffold and l:
        list_line_gff.append(l)
        l = fl.readline()
    fl.close()
    return list_line_gff

def get_gff_obj(line):
    result = gff_parser(line)
    gene = Gene(line)
    gene.scaffold = result.group("scaffold")
    gene.feature = result.group("feature")
    gene.start= int(result.group("start"))
    gene.end = int(result.group("end"))
    gene.strand = result.group("strand")
    gene.gene_number = int(result.group("gene_number"))
    gene.len_valide()
    
    return gene

def gff_list_builer(gff_file, scaffold):
    """
    input : name of the gff file and of the scaffold that is currently analysed
    output : list of objet Gene() containing the information of the each line of the file
    """
    fl = open(gff_file, 'r')
    
    list_gff_obj = []
    
    l = fl.readline()
    
    while l[:len(scaffold)] != scaffold: #gff file is sorted by scaffold so we use this feature in order to ot scan all the file  
        l = fl.readline()
        
    while l[:len(scaffold)] == scaffold:
        #print l
        
        result = gff_parser(l)
        gene = Gene()
        gene.gff_line = l
        
        gene.scaffold = result.group("scaffold")
        gene.feature = result.group("feature")
        gene.start= int(result.group("start"))
        gene.end = int(result.group("end"))
        gene.strand = result.group("strand")
        gene.gene_number = int(result.group("gene_number"))
        list_gff_obj.append(gene)
        #print gene.__dict__
        l = fl.readline()
    fl.close()
    #print list_gff_obj
    return list_gff_obj

    
def gff_parser(line):
    pattern = re.compile("^(?P<scaffold>[^\t]*)\t[\w]+\t(?P<feature>\w+)\s(?P<start>\d+)\s+(?P<end>\d+)\s+[\d.-]+\t(?P<strand>[+-])\t[\d.-]+\tID=\w+\|(?P<gene_number>\d+);")
    #print line
    result = pattern.match(line)
    #print result.groupdict()
    return result

class Domain:
    #Threshold overlap !! 
    threshold = 0.1 # 10% of overlaping after that  the domains are considered as overlaping
    
    def __init__(self, domain_number, ali_from, ali_to, env_from, env_to, domain_name, domain_acc, e_value, score, line):
        self.domain_number = int(domain_number)
        self.ali_from = int(ali_from)
        self.ali_to = int(ali_to)
        self.env_from = int(env_from)
        self.env_from = int(env_from)
        self.domain_name = domain_name
        self.domain_acc = domain_acc
        self.e_value = float(e_value)
        self.score = float(score)
        self.line = line
        
    def overlap(self, d):
        if self.ali_from <= d.ali_from <= self.ali_to: #overlap
            if self.ali_to - d.ali_from > (self.ali_to - self.ali_from)*Domain.threshold:
                return True
            else:
                return False
        elif self.ali_from <= d.ali_to <= self.ali_to: #overlap
            if d.ali_to - self.ali_from > (self.ali_to - self.ali_from)*Domain.threshold:
                return True
            else:
                return False
        elif d.ali_from <= self.ali_from and self.ali_to < d.ali_to: #overlap est self anglobé dans d
            return True
        else: #overlap pas 
            return False
    
        
def get_list_scaffold(table_hmm):
    ## input : Name of the table hmm file
    ## output : List of all the scaffold of the table
    
    set_of_scaffold = set()
    fl_hmm = open(table_hmm, 'r')
    for line in fl_hmm:
        if line[0] == "#":
            continue
        set_of_scaffold.add(line[:line.index('|')])
    fl_hmm.close()
    return list(set_of_scaffold)


class Gene:
    #treshold size of gene
    length_min = 80
    length_max = 900
    
    #threshold distance for tandem
    dist_inf = -20
    dist_max = 150
    
    #Allowance of integrity loss of domain (5%) arbitrary number 
    allowance = 0.05
    
    
    def __init__(self, gff_line):

        self.gff_line = gff_line
        
        
        self.scaffold = None
        self.gene_number = None
        self.feature = None
        self.start= None
        self.end = None
        self.strand = None
        self.domain_number = None
        self.domain = []
        self.domain_discarded = []
        self.best_domain = []
        self.rasta_domain = []
        
        #Validity features
        self.len_val = None
        self.post_va = []
        self.prev_va = []
        self.statue = None #statue dit si le gene valide ou non.. il est pas valide dans le cas ou le resize ne marche pas

        self.prev_nonva = []
        self.post_nonva = []
        self.TA_adj = None
        
    def get_best_domain(self):
        #print 'nb domain', len(self.domain)
        #for dd in self.domain:
            #print dd.line
        if len(self.domain) == 1:
            self.best_domain.append(self.domain[0])
            return
        domains = self.domain[:]
        #print domains
        removed = []
        while len(domains) > 0:
            removed = []
            do = domains.pop(0)
            for d in domains: #compare to do
                if do.overlap(d):
                    removed.append(d) #d est éliminé de la liste
                    #print do.domain_number,do.ali_to, 'overlap with ', d.domain_number, d.ali_to 
                    if do.e_value >= d.e_value: #donc si le match de do est de meilleur qualité 
                        #print 'elimination ', do.domain_number,do.ali_to
                        do = d #d devient le nouveau do
                        
                    #else:
                        #print 'elimination ',d.domain_number, d.ali_to 
            #print do.domain_number,do.ali_to  ,'append '
            self.best_domain.append(do)
            for rd in removed:
                domains.remove(rd)
        #print 'nb best domain\n', len(self.best_domain)
        #print self.best_domain
        #print '='*20

    def has_rasta_domain(self, rasta_domain):
        for d in self.domain:
            if d.domain_name in rasta_domain or d.domain_acc in rasta_domain:
                self.rasta_domain.append(d)
                #print 'rasta_domain fouond!'
    
    def resize(self, dico):
        fl = dico['fl']
        line = dico["line"]
        start_codon = dico["codon_start"]
        #print start_codon
        
        scaffold_gnb = self.scaffold+'|'+str(self.gene_number)
        #print 'retreiving info for ', scaffold_gnb
        seq, fl, line = orf.get_fast_fasta(fl,line, scaffold_gnb)
        ##remet les line et fl dans le dico sinon c'est pas cool
        dico['fl'] = fl
        dico['line'] = line
        
        start_po = orf.codon_finder(start_codon, seq["data"])
        if not start_po: #si aucun start trouvé.. dans ce cas là on return False
            return False
            
        for start in start_po:
            if self.length() - start <= Gene.length_max: #on check si la taille du gene avec ce start est ok
                #print self.length() - start , 'est ok reste checker intégrité of the domain !! '
                if self.domain_affected(start) or self.length() - start < Gene.length_min: #si ça affect l'intégrité des domain return True ou si la taille est devenu trop petite
                    #print 'mettre se gene dans un catégorie éliminé '
                    return False 
                else:
                    if self.strand == '+':
                        self.start += start 
                    else:
                        self.end -= start 
                    self.len_valide()
                    if not self.len_val:
                        print self.length()
                    return True #le gène est resizé et validé car au moins un domaine est intègre 
        
        
        #print 'start position', start_po
        #print seq
        return True 
    
    def domain_affected(self, start):
        #print 'len domain', len(self.domain)
        start = start/3 #transformation en start aa
        domain_rm = []
        #print 'nouveau :', start
        for  do in self.domain:
            #print 'Domain',do.ali_from,'-->',do.ali_to
            if start > do.ali_from+(do.ali_to-do.ali_from)*Gene.allowance :
                #print 'domain discard'
                domain_rm.append(do)
            #else:
                #print 'domain is keep'
        for do in domain_rm:
            self.domain_discarded.append(do)
            self.domain.remove(do)
        #print 'len domain', len(self.domain)
        #print 'len domain discard ', len(self.domain_discarded)
        if len(self.domain) ==0:
            #print 'gene va être discard'
            return True #true le gene est affecté par le changement de start, itegrité n'est pas suffisante! 
        else:
            #print 'gene ok'
            return False
        
    def length(self):
        return self.end-self.start+1
    
    def get_adj(self):
        return self.post_va + self.prev_va + self.prev_nonva + self.post_nonva
    
    def add_adj(self, gene_group = [], flag={'len_va':0, 'len_nonva':0}):
        for adj in self.get_adj():
            if adj in gene_group: #si les genes ont deja été traité je les passe ! pour éviter ue boucle infini
                continue
            
            gene_group.append(adj)
            if adj.len_val is True:
                flag['len_va'] += 1
            else:
                flag['len_nonva'] += 1
                
            if adj.sum_link > 1: 
                """ si le gene adjacent a un seul lien alors c'est le gene "self" qu'on test, du coup pas besoin d'aller chercher plus loin on ajoute simplement le gene a la liste adj a ici plus que 1 lien car tous les gènes qu'on test on forcement au moin une connexion du coup on va cherche """
                adj.add_adj(gene_group, flag) 
        return gene_group
    
    
    def len_valide(self):
        if Gene.length_min <= self.length() <= Gene.length_max:
            self.len_val = True
        else:
            self.len_val = False
        
    def tandem_gene(self, gene):
        if Gene.dist_inf<gene.start-self.end<Gene.dist_max:
            #print self.gene_number,'and', gene.gene_number,'are close :', gene.start-self.end , self.strand
            if self.len_val and gene.len_val:
                #print '->>Their length are ok... '
                self.post_va.append(gene)
                gene.prev_va.append(self)
            else:
                self.post_nonva.append(gene)
                gene.prev_nonva.append(self)
            return True
        elif gene.start-self.end>Gene.dist_max:
            return False #false le gene start trop loin 

        
    def sum_link(self): #sum of link with other gene
        return len(self.prev_va) + len(self.post_va) + len(self.prev_nonva) + len(self.post_nonva)
            
    def show(self, attribut=True):
        if attribut is True: #True means here all attribut have to be shown 
            attribut = dir(self)
        for a in attribut:
             print a, getattr(self, a) 
            
            
class ListGene:
    def __init__(self):
        self.liste= []
        self.list_minus = []
        self.list_plus = []
        self.numbers = []
        self.list_pair = []
        self.ggrouped = [] #gene already grouped
        self.ggrouped_rasta = [] #gene already grouped with rasta domain only
        self.ggrouped_rasta_mixt = [] #gene groupe in a group mixt with rasta and not rasta.. 
        self.linked = [] #gene that have a link with another one
        self.group_len_va = [] #list of group so liiste of liste of obj gene ;-) 
        self.group_rasta = [] #liste of group of gene with only rasta
        self.group_rasta_mixt = [] #liste of group with at least one rasta 
        self.lonely_gene = [] #gene that are not adj to any other one
        self.non_valides = [] #liste de gene non valide determiner par check size fonction
        
    def listappend(self, gene):
        self.liste.append(gene)
        if gene.strand == '+':
            self.list_plus.append(gene)
        else:
            self.list_minus.append(gene)
        self.numbers.append(gene.gene_number)
        
    def appendomain(self, domain, nb_gene):
        index = self.numbers.index(nb_gene)
        self.liste[index].domain.append(domain)
        
    def check_size(self, dico_seq):
        #Dico_seq est composé du fl du fichier fna avec les sequence fasta des gene prédit et line qui correspond à la derière ligne 
        #Mise en dico pour pas avoir besoind de le retourner 
        self.liste = sorted(self.liste,  key=attrgetter('gene_number'))
        
        for gene in self.liste:
            if gene.feature == 'ORF':
                continue
            #print gene.gene_number
            if  Gene.length_min <= (gene.end - gene.start) +1 <= Gene.length_max:
                #print 'valide'
                gene.statue = True #statue dit si le gene valide ou non.. il est pas valide dans le cas ou le resize ne marche pas
                continue
            else:
                #print 'resize ', gene.gene_number
                if Gene.length_min > (gene.end - gene.start) +1:
                    #print 'gene trop petit'
                    gene.statue = False
                    self.non_valides.append(gene)
                
                elif gene.resize(dico_seq): #va resizer et si c'est pas possible renvoit un False
                    #print 'gene ', gene.gene_number, 'resize GOOD ! '
                    gene.statue = True
                else:
                    #print 'gene ', gene.gene_number, 'resize FAIL ! '
                    gene.statue = False
                    self.non_valides.append(gene)
        #ATTENTION need to discard the gene in the main liste and also in minus and plus liste !!!!!!!!!!
        for g in self.non_valides:
            #print 'INVALIDE',g.gene_number
            self.liste.remove(g)
            if g.strand == '+':
                self.list_plus.remove(g)
            else:
                self.list_minus.remove(g)
        
    def set_best_rasta_domain(self, list_rasta):
        #### SET the best domain and rasta_domain
        for g in self.liste:
            #print '+++++++'*5
            g.get_best_domain()
            g.has_rasta_domain(list_rasta)
            #print '+++++++'*5
        
    def get_adj(self):
        self.list_plus = sorted(self.list_plus,  key=attrgetter('end'))
        self.list_minus = sorted(self.list_minus,  key=attrgetter('start'))

        
        
        ####
        #TODO Ameliorer mettre un break quand les gene sont de toute facon trop loin 
        
        #Strand plus !
        for gi in range(len(self.list_plus)):
            #print self.list_plus[gi].gene_number,  self.list_plus[gi].len_val
            for gpost in self.list_plus[gi+1:]:

                if self.list_plus[gi].tandem_gene(gpost): #it is not a simple test, tandem_gene() will store the information of prev and post 
                    self.linked.append(self.list_plus[gi])
                
        
        #Strand minus!
        for gi in range(len(self.list_minus)):
            for gpost in self.list_minus[gi+1:]: #small trick with the list to go bac
                if self.list_minus[gi].tandem_gene(gpost): #it is not a simple test, tandem_gene() will store the information of prev and post 
                    self.linked.append(self.list_minus[gi])
               
        
        
    def get_lonely_gene(self):
        for g in self.liste:
            if not g.get_adj(): #si get adj renvoit une liste vide
                self.lonely_gene.append(g)
            
            
    def get_info(self, composition):
             
        try:
            composition["lonely_gene"] += len(self.lonely_gene)
        except KeyError:
            composition["lonely_gene"] = len(self.lonely_gene)
            

        rasta_lonely = 0
        for lonely in self.lonely_gene:
            if len(lonely.rasta_domain) != 0:
                rasta_lonely += 1
        try:
            composition["lonely_rasta"] += rasta_lonely
        except KeyError:
            composition["lonely_rasta"] = rasta_lonely
        
        
        
    def group(self, composition):
        for gene in self.linked:
            if gene in self.ggrouped:
                continue
            #print '\nAaa',gene.gene_number, ':'
            gene_group = []
            flag={'len_va':0, 'len_nonva':0}  #nb de gene avec len valide et non valide permet de savoir si le groupe est mixt ou non

            gene.add_adj(gene_group, flag)
            
            rasta_flag = True
            rasta_flag_group_mixt =False #group avec au moins un gene rasta ! 
            for g in gene_group:
                if len(g.rasta_domain) == 0:
                    rasta_flag = False
                elif len(g.rasta_domain) != 0:
                    rasta_flag_group_mixt = True
                
            self.ggrouped += gene_group #on ajoute les genes deja traité a la liste grouped pour qu'il ne soit pas traité par la suite
            
            if rasta_flag:
                self.ggrouped_rasta += gene_group #gene already grouped with rasta domain only
            elif rasta_flag_group_mixt:
                self.ggrouped_rasta_mixt += gene_group #append gene with are in a group mixt 
            
            #for g in gene_group:
                #print '\n',g.gene_number, g.start, '-->' ,g.end
                #for adj in g.get_adj():
                    #print '\t',adj.gene_number,
            #print '\n',flag
            if rasta_flag is True:
                 #print 'rasta group'
                self.ggrouped_rasta += gene_group
                try:
                    composition['rasta'][len(gene_group)] += 1
                except KeyError:
                    composition['rasta'][len(gene_group)] = 1
            
            #ajout ele groupe a rasta group mixte seulement si il est pas 100% composé de rasta 
            elif rasta_flag_group_mixt:
                #print 'rasta group mixt'
                try:
                    composition['rasta_mixt'][len(gene_group)] += 1
                except KeyError:
                    composition['rasta_mixt'][len(gene_group)] = 1
                
                
                
            if flag['len_va'] == 0:
                
                try:
                    composition['nonva'][len(gene_group)] += 1
                except KeyError:
                    composition['nonva'][len(gene_group)] = 1
                    
            elif flag['len_nonva'] == 0: #donc group homogene len valide
                #Append the group to list of group  of gene to write then in a file
                self.group_len_va.append(gene_group)
                if rasta_flag:
                    self.group_rasta.append(gene_group)
                if rasta_flag_group_mixt:
                    self.group_rasta_mixt.append(gene_group)
                    
                try:
                    composition['va'][len(gene_group)] += 1
                except KeyError:
                    composition['va'][len(gene_group)] = 1
                    
            elif flag['len_va'] >0 and flag['len_nonva']>0:
                try:
                    composition['mixt'][len(gene_group)] += 1
                except KeyError:
                    composition['mixt'][len(gene_group)] = 1
                    
    def write_to_file(self, file_name, liste_name='group_rasta'):
        #write les group valide pour le moment
        fw = open(file_name, 'a')
        sep = '\t'
        if getattr(self, liste_name) :
            #print self.ggrouped_rasta
            fw.write(self.liste[0].scaffold+'\n')
            getattr(self, liste_name)
        for i, group in enumerate(getattr(self, liste_name)):
            fw.write('Group '+str(i+1)+'\n')
            #print group
            for gene in group:
                line = sep+gene.scaffold+'|'+str(gene.gene_number)+sep+str(gene.start)+sep+str(gene.end)+sep+gene.strand+sep+'Rasta:'

                
                for d in gene.rasta_domain:
                    line += d.domain_name+' '+d.domain_acc+','
                if not  gene.rasta_domain:
                    line += str(None)
                line += sep+'|Best:'
                for d in gene.best_domain:
                    line += d.domain_name+' '+d.domain_acc+','
                line += '\n'
                fw.write(line)
            fw.write('\n')
        fw.close()
    
    def count_rasta_domain(self, dico_rasta):
        #count the number of each rasta domain 
        
        for g in self.lonely_gene:
            for dr in g.rasta_domain:
                try:
                    dico_rasta[dr.domain_name]['lonely_rasta'] += 1
                except KeyError:
                    dico_rasta[dr.domain_acc]['lonely_rasta']+= 1 
        for g in self.ggrouped_rasta:
            for dr in g.rasta_domain:
                try:
                    dico_rasta[dr.domain_name]['inside_rasta_group'] += 1
                except KeyError:
                    dico_rasta[dr.domain_acc]['inside_rasta_group'] += 1 
        for g in self.ggrouped_rasta_mixt:
           for dr in g.rasta_domain: 
                try:
                    dico_rasta[dr.domain_name]['inside_group_mixt'] += 1
                except KeyError:
                    dico_rasta[dr.domain_acc]['inside_group_mixt'] += 1 
                
                
                
                
def write_csv_domain(csvfile, dico_domain):
    sep=';'
    name = ''
    lonely = ''
    inside_group_mixt = ''
    inside_rasta_group = ''
    total = ''
    for domain_name in sorted(dico_domain.keys()):
        name += sep+domain_name
        lonely += sep+str(dico_domain[domain_name]['lonely_rasta'])
        inside_rasta_group += sep+str(dico_domain[domain_name]['inside_rasta_group'])
        inside_group_mixt += sep+str(dico_domain[domain_name]['inside_group_mixt'])
        total += sep+str(dico_domain[domain_name]['lonely_rasta']+dico_domain[domain_name]['inside_rasta_group']+dico_domain[domain_name]['inside_group_mixt'])
    
    fl = open(csvfile, 'w')
    fl.write('domaine_name'+name+'\n')
    fl.write('domain_in_rasta_group'+inside_rasta_group+'\n')
    fl.write('domain_in_group_mixt'+inside_group_mixt+'\n')
    fl.write('domain_lonely'+lonely+'\n')
    fl.write('total'+total)
    fl.close()
        
        
        
        
        
        
                
                
                