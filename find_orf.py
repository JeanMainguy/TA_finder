# coding: utf-8
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio import SeqIO
from operator import attrgetter

def get_list_orf(scaffold_file, scaffold, threshold, id_genTable, fl, line):  
    #seq = get_fasta(scaffold_file, scaffold)
    #seq = get_ifasta(scaffold_file, scaffold, index)
    seq, fl, line = get_fast_fasta(fl, line, scaffold)
    

    #seq = {}
    #handle = open(scaffold_file, "r")
    #for record in SeqIO.parse(handle, "fasta"):
        #if record.id == scaffold:
            ##print record.seq
            
            #seq['data'] = str(record.seq)
            #seq['description'] = record.id
            #handle.close()
            #break
    if not seq:
        print scaffold, 'not found'
    #print scaffold, '\n ->longueur',len(seq['data'])
    return findORF(seq, threshold,id_genTable), fl, line
    #return seq 
    
def get_index_sca(scaffold_file, sca_list):
    
    fl = open(scaffold_file, 'r')
    line = fl.readline()
    i=0
    isca=0
    index={}
    while line and isca < len(sca_list):
        if line[0]=='>':
            if line[1:-1] == sca_list[isca]:
                index[line[1:-1]] = i
                isca +=1
        i +=1
        line = fl.readline()
    fl.close()
    if len(index) != len(sca_list):
        print 'PROBLEM in scaffold index'
    return index

def get_ifasta(scaffold_file, scaffold, index):
    fl = open(scaffold_file, 'r')
    isca = index[scaffold]
    #print isca
    i = 0
    line = fl.readline()
    seq = {'description':scaffold, 'data':''}
    while i < isca:
        i+= 1
        fl.readline()
    line = fl.readline() #première ligne de sequence.. 
    while line and line[0] != '>':
        seq['data'] += line.rstrip()
        line = fl.readline()
    fl.close()
    return seq
def get_fast_fasta(fl,line, scaffold):
    #line est la dernière ligne lu dans le fichier et doit correspondre à un entête > 
    #la liste de scaffold est trié donc normalement y'a pas de problem.. Et les scaffold cherché devrait être dans le bonne ordre
    seq = {'description':scaffold, 'data':''}
    #print scaffold+'='+line
    while line and line[:len(scaffold)+1] != '>'+scaffold:
        line = fl.readline()
    line = fl.readline()
    while line:
        if line[0] == '>':
            return seq, fl, line
        seq['data'] += line.rstrip()
        line = fl.readline()
        
def get_fasta(scaffold_file, scaffold):
    fl = open(scaffold_file, 'r')
    line = fl.readline()
    seq = {'description':scaffold, 'data':''}
    while line[:len(scaffold)+1] != '>'+scaffold:
        line = fl.readline()
    if line:
        
        line = fl.readline()
        while line and line[0] != '>':
            seq['data'] += line.rstrip()
            line = fl.readline()
    else:
        print 'scafold not found !'
    fl.close()
    return seq

def findORF(seq,threshold,id_genTable):
    """This function determine the orf of a nucleic sequence 
    This function call the orf_by_frame() function  for each frame and give it to it different parameter especially the frame and the information if the sequence is reverse or not and finally gather the result in a dictionary

    Args:
        seq : sequence where the orf will be determine. SEqObjet from Biopython
        threshold : the minimum size of orf
       # nub_genTable : id Genetic code Table 
    Returns:
        A list of ORFs as a list of dictionaries
    
    """
    #print seq['data']
    genTable = getGeneticCode(id_genTable)
    #print genTable['start']
    
    orf_zero = Orf()
    orf_zero.set_contig(seq)
    
    seq_reverse = dict(seq)
    seq_reverse = complement_reverse(seq_reverse) 
    
    listORF = ListGene()

    for frame in [1,2,3]:

        #print('The program is looking at orf located on frame ', frame)
        orf_by_frame(listORF, seq, genTable, frame, threshold)
        #print(len(listORF), 'orfs found') 
        #print('The program is looking at orf located on frame -', frame, '')
        orf_by_frame(listORF, seq_reverse, genTable, frame, threshold, rev=-1)
    #if listORF:
    #listORF[0].set_contig(seq)
    return listORF


def orf_by_frame(listORF, seq, genTable,  frame, threshold, rev=1):
    """This function find the orf in a sequence for a given frame
    
    First the function determine all the stop codon of the sequence in frame by calling codon_finder(). 
    Then it will find all the start codon flan by two stop codon by calling codon_finder 
    And finally will build the dictionaries with the iformation relative to each orf by calling the remplissage() function
    This function is written by Theo Falgarone. 

    Args:
        seq : sequence where the orf will be determine
        frame : give the frame of the search
        threshold : the minimum size of orf
        genTable : Genetic code Table 
        rev : 1 or -1 to know if the sequence has been completed reverted in order to fill correctly the frame. 
        
    Returns:
        The list of the Orf of the frame as a list of dictionaries
	
	"""
    
    length_seq = len(seq['data'])
    pos_stop =codon_finder(genTable['stop'], seq['data'], frame=frame)
    
    if not pos_stop: #au cas ou il y a pas de codon stop dans la sequence... on sait jamais 
        if length_seq > threshold:
            pos_start = []
            pos_start = codon_finder(genTable['start'], seq['data'], frame=frame) #research of starts in the whole sequence ! 
            # ATTENTION ajout du debut de la sequence aux start... start imaginaire donc, mais necessaire pour coller au gene predit ! :-(
            pos_start = [frame-1] + pos_start
            if pos_start:
                #On determine le stop
                modulo_start = (frame-1)%3
                #TODO completement améliorable !! 
                #print 'frame', frame, ' modulo start =', modulo_start
                if modulo_start == (length_seq-3)%3:
                    stop = length_seq-3
                elif modulo_start == (length_seq-3-1)%3:
                    stop = length_seq-4
                elif modulo_start == (length_seq-3-2)%3:
                    stop = length_seq-5
                    
                remplissage(seq, listORF, genTable, pos_start, stop+2, frame, rev, False, border=True)
        return listORF
    
   
    #research of potential orf between 0 position and the first stop -->  the program doesn't manage circular DNA 
    if pos_stop[0]+3 > threshold: #check if the length from 0 to the first stop is bigger than the trhresold
        pos_start = []
        pos_start = codon_finder(genTable['start'], seq['data'], frame=frame, inf=0, sup=(pos_stop[0]+3-threshold)) #research of starts from 0 to (first stop - thresold) 
        # ATTENTION ajout du debut de la sequence aux start... start imaginaire donc, mais necessaire pour coller au gene predit ! :-(
        pos_start = [frame-1] + pos_start
        if pos_start:
            listORF = remplissage(seq, listORF, genTable, pos_start, pos_stop[0]+2, frame, rev, border=True)
    
    for i in range(len(pos_stop)-1): #recherche entre le stop i et i+1 donc on s'arrete a un stop avant la fin 
        
        if (pos_stop[i+1]+3-(pos_stop[i]+3)) < threshold: #premier +3 pour inclure le codon stop dans la longueur final comme dans ncbi orf finder... 
            continue
        
        else:
            pos_start = []
            pos_start = codon_finder(genTable['start'], seq['data'], inf=pos_stop[i]+3, sup=(pos_stop[i+1]+3-threshold))
            if pos_start:
                listORF = remplissage(seq, listORF, genTable, pos_start, pos_stop[i+1]+2, frame, rev)
    
    #Search of orf from the last stop codon to the end of the seq
    if length_seq-pos_stop[-1]+3 > threshold:
        pos_start = []
        pos_start = codon_finder(genTable['start'], seq['data'], inf=pos_stop[-1]+3, sup=(length_seq-threshold)) #research of starts from 0 to (first stop - thresold) 
        
        if pos_start:
            #looking for the appropriate stop:
            i = 1
            modulo_start = pos_start[-1]%3
            #TODO completement améliorable !! 
            #print 'frame', frame, ' modulo start =', modulo_start
            if modulo_start == (length_seq-3)%3:
                stop = length_seq-3
            elif modulo_start == (length_seq-3-1)%3:
                stop = length_seq-4
            elif modulo_start == (length_seq-3-2)%3:
                stop = length_seq-5
            else:
                exit('NO STOP HAVE BEEN FOUND !!! ERROR exit()')
            
            ##while (pos_start[-1] + 3*i)<= length_seq-2: #-2 car on veut déterminer la position du codon stop qui doit être dans la même frame que le start. Ensuite on ajoute 2 pour avoir la position du dernier codon lorsqu'on rempli.. 
                ##i +=1
            ##stop = pos_start[-1] + 3*(i)
           
            remplissage(seq, listORF, genTable, pos_start, stop+2, frame, rev, False, True)
            #print 'ncbi end', listORF[-1].ncbi_end
            #print 'ncbi start', listORF[-1].ncbi_start
      
    
    
    return listORF

def codon_finder(liste, seq, frame=1, inf=0, sup='Not defined'):
    """Extract position of given codons in a sequence and in a given frame
    
    By default this function search codon of the given list in the whole sequence. However it possible to give the position in the sequence where codon has to be found with sup and inf parameter. 

    Args:
        liste   : liste of codons 
        seq     : DNA sequence
        frame   : give the frame of the search
        inf     : position in sequence where the function has to start to search
        sup     : position in sequence where the function has to stop to search
        
    Returns:
        Return a list of the position of the given codon in the sequence
    """
    
    position = []
    if sup=='Not defined':
        sup = len(seq)-2
    for i in range(inf+frame-1, sup, 3):
        if seq[i:i+3] in liste:
            position.append(i)
    return position

        
def remplissage(seq, listORF, genTable, pos_start, stop, frame, rev, complet=True, border=False):
    """This function fill the orf dictionary 
    
        This function for each start between two stop codon will create a dictionary and add all relevant information

    Args:
        seq : sequence where the orf will be determine
        listORF : list of ORF of the frame
        frame : give the frame of the search
        threshold : the minimum size of orf
        genTable : Genetic code Table 
        pos_start : list of start position between two stop codon
        stop : list of two following stop codon
        rev : 1 or -1 to know if the sequence has been completed reverted in order to fill correctly the frame. 
        
    Returns:
        ListORF : The completed list of Orf as a list of dictionaries
        """
    
    orf = Orf(frame=frame*rev, start=pos_start,end=stop, contig=seq, complet=complet, border=border)
    
    
    listORF.listappend(orf)

        
    return listORF   
    
class Orf:
    #treshold size of gene
    minn = 80
    maxx = 900
    contig = None
    def __init__(self, frame=1, start=[0],end=0, contig={}, complet=True, NCBI_id=11, border=False):
        
        self.frame = frame
        self.contig = contig
        self.start = start
        self.end = end
        self.first_start = start[0]
        
        #Conversion to get start and stop found with official orf finder
        if self.frame < 0:
            len_contig = len(contig['data'])
            self.ncbi_end = len_contig-end
            self.ncbi_start = [len_contig-x for x in start]
            self.strand = '-'
            self.gff_start = self.ncbi_end
        else:
            self.ncbi_start = [x+1 for x in start]
            self.ncbi_end = end+1
            self.strand = '+'
            self.gff_start = self.ncbi_start[0]

        self.complet = complet
        self.border = border
        self.NCBI_id = NCBI_id #May change but here we only deal with code 11 
        self.gff_obj = None
        self.annotation = False
        self.true_gene = False #ATTENTION True if the annotation of the gene is strong !!!
        
    def first_start(self):
        return self.start[0]
    
    def resize(self):
        #TODO Amelioration pour que ça elimine direct quand le premier start donne une sequence trop petite 
        #print 'resize pour', self.end
        if  Orf.minn <=  self.data_len() <= Orf.maxx:
            #print '\ttaille initial ok !'
            return True #C'est bon l'orf est ok 
        start_rm = []
        for i in range(len(self.start)):
            if not Orf.minn <= self.data_len(i) <= Orf.maxx:
                start_rm.append(self.start[i])
        if len(start_rm) == len(self.start):
            print start_rm ,'==' ,self.start
            print self.data_len() ,'tous les starts sont mauvais'
            return False #meaning that every start will be removed then we don't want this orf
        else:
            for start in start_rm:
                self.start.remove(start)
            #reset the start .. 
            self.first_start = self.start[0]
            return True
    
    def set_contig(self, seq):
        Orf.contig = seq
        
    def data_nt(self, index_start = 0):
        return self.contig['data'][self.start[index_start]:self.end+1] ####!!!### Verifier !
    def data_len(self, index_start=0):
        return self.end - self.start[index_start]
    def prot_len(self, index_start=0):
        return self.data_len(index_start) /3
    
    def protein(self):
        seq = {'data':self.data_nt()[3:]}
        return 'M'+ translate(seq)['data']##ATTENTION every prot seq start with a M... 
    
    def write_gff(self, fl, nb):
        #Exemple of gff_line
        #ICM0104MP0310_1000150	GFPMX	CDS	9184	9960	.	-	0	ID=ICM0104MP0310_1000150|9
        if self.strand == '+':
            start = str(self.start[0]+1)
            end = str(self.ncbi_end)
        else:
            start = str(self.ncbi_end)
            end = str(len(self.contig['data'])-self.start[0])
        scaffold = Orf.contig['description']
        line = scaffold+'\tCDS\tORF\t'+start+'\t'+end+'\t.\t-\t0\tID='+scaffold+'|'+str(nb)+';partial='+str(self.complet)+';'+str(len(self.data_nt()))+'\n'
        fl.write(line)
        
        
   
    def write_faa(self,fl, nb):
        #to write seq prot in a fastafile fl is handle file et nb is the number of the gene
        #seq = self.protein()
        seq = self.protein()
        fl.write('>'+Orf.contig['description']+'|'+str(nb)+'\n')
        
        if len(seq)<=60:
            fl.write(seq+'\n')
            return
        for i in range(len(seq)/60):
            fl.write(seq[i*60:(i+1)*60]+'\n')
        if len(seq)%60 != 0:
             fl.write(seq[(i+1)*60:]+'\n')
        
    def show(self, attribut=True):
        if attribut is True: #True means here all attribut have to be shown 
            attribut = dir(self)
        for a in attribut:
             print a, getattr(self, a) 
             
    def get_gff_info(self, g):
        self.gff_obj = g
        #self.annotation = g.annotation
        #self.true_gene = g.true_gene
        #if g.true_gene:
            #self.ann_start = g.start 
        
    def set_index(self, i):
        self.index=i
    def get_gff_end(self):
        if self.frame<0:
            return self.ncbi_start
        return self.ncbi_end
        #record = SeqRecord(self.protein(),self.contig['description']+'|'+nb, '', '')

        #SeqIO.write(mito_frags, fl, "fasta")
        
class ListGene:
    ###THRESHOLD pour recupérer les orf adjacent 
    dist_inf = -20
    dist_max = 150
    
    def __init__(self):
        self.liste= []
        self.list_minus = []
        self.list_plus = []
        self.list_pair = []
        self.list_gene_hmm = [] #orf correspondant au gene hmm qui ne sont pas dans un groupe
        self.set_adj_orf = set() #on utilise un set au ca il y a un orf soit adj de lusieur hmm gene il soit pas ajouté plusieur fois
        
    def write_faa_gff(self, index, gff_file, faa_file): #index_gene is the number of gene already in the scaffold to then be able to give a number to the adjacent orf that ientify them among
        gffl = open(gff_file, 'a')
        faafl = open(faa_file, 'a')
        for orf in list(self.set_adj_orf):
            index += 1
            orf.write_faa(faafl, index)
            orf.write_gff(gffl, index)
        gffl.close()  
        faafl.close()
        
    def listappend(self, gene):
        self.liste.append(gene)
        if gene.strand == '+':
            self.list_plus.append(gene)
        else:
            self.list_minus.append(gene)
        

            
            
    def find_the_orf(self,g):
        """
        Arg : 
        """
        flag_stop = False
        flag_ok = False
        if g.strand == '+':
            list_orf = self.list_plus
            corres  = {'ncbi_end':'end', 'ncbi_start':'start'}
    
        else:
            corres = {'ncbi_end':'start', 'ncbi_start':'end'}
            list_orf = self.list_minus
            
        for orf in list_orf:
            #Trouver le stop associe a g
            if getattr(g, corres['ncbi_end']) == orf.ncbi_end:
                #print 'stop trouvé'
                flag_stop = True
                #trouver les start associe au stop trouve precedement
                if getattr(g, corres['ncbi_start']) in orf.ncbi_start:
                    flag_ok =True
                    orf.get_gff_info(g) #rtreive info of the gff obj to the orf !!!
                    self.list_gene_hmm.append(orf)
                    
                    break
                    #print 'start trouvé'
                else: # on a le stop mais pas le start
                    
                    #print g.gff_line
                    #print 'g start',getattr(g, corres['ncbi_start']), '== longueur contig', len(list_orf[0].contig['data']), 'or == ',1
                    if getattr(g, corres['ncbi_start']) == len(list_orf[0].contig['data']) or getattr(g, corres['ncbi_start']) == 1:
                        print 'PROBLEME DE MODULO ... ', orf.border
                        
                        flag_ok = True
                        orf.get_gff_info(g) #rtreive info of the gff obj to the orf obj!!!
                        self.list_gene_hmm.append(orf)
                        break
                    else:
                        #print 'start pas trouvé donc potentiellement c un mauvais stop qui a ete trouve !<!<!<!!<!<!!<!<!!<!<!<!!<!<!<!!<!<!<!<!<!!<'
                        #print orf.ncbi_start,'-->', orf.ncbi_end , 'border',orf.border,'complet', orf.complet, 
                        flag_stop = False
        if not flag_stop:
            if getattr(g, corres['ncbi_end']) == len(list_orf[0].contig['data']) or getattr(g, corres['ncbi_end']) == 1:
                for ob in list_orf:
                    if ob.border:
                        if getattr(g, corres['ncbi_start']) in ob.ncbi_start:
                            #print 'ça colle au niveau du start de nb ', g.gene_number
                            #print 'orf correspondant', ob.ncbi_start, '-->', ob.ncbi_end
                            flag_ok = True
                            orf.get_gff_info(g) #rtreive info of the gff obj to the orf obj !!!
                            self.list_gene_hmm.append(orf)

        if not flag_ok:
            print 'STOP PAS TROUVEEEE \o/' 
            print g.gff_line
            #exit('BOUUUUUUUHH !')
            fl = open('Error.out', 'a')
            fl.write(g.gff_line)

        
        
    def get_adj_orf(self):
        #J'ai décidé d'utilisr les vrai start et stop des deux brins et pas les conversions. Donc chaque brin a sa propre métrique
        #SORT ! 
        
        plus_by_end = sorted(self.list_plus,  key=attrgetter('end'))
        minus_by_end = sorted(self.list_minus,  key=attrgetter('end'))        
      
        for orf in self.list_gene_hmm:
            if orf.strand == '+':
                l_by_end = plus_by_end 
            else:
                l_by_end = minus_by_end 
          
            ##upstream
            index = l_by_end.index(orf)
            i = 1
            #print '\nORF :',orf.first_start, '-->', orf.end, orf.strand 
            #print 'UPSTREAM'
            while index +i < len(l_by_end):
                if not l_by_end[index+i].resize(): #c'est pas juste un test ça va resizé et si jamais y'a pas de start valide ça return False...
                    i += 1 # on passe donc au tour suivant
                    continue
                    
                #print 'after',l_by_end[index+i].data_len() 
                
                if ListGene.dist_inf <= l_by_end[index+i].first_start - orf.end <= ListGene.dist_max:
                    
                    self.set_adj_orf.add(l_by_end[index+i])
                    #print l_by_start[index+i].first_start,'-->', l_by_start[index+i].end, l_by_start[index+i].strand
                    #print 'distance ok',  l_by_end[index+i].first_start, '-' ,orf.end,'=',l_by_end[index+i].first_start - orf.end,'\n'
                    
                elif l_by_end[index+i].first_start - orf.end > ListGene.dist_max: #si l'orf adj start trop loin du end ww
                    break #on peut faire ça car les listes sont sorted :-) 
                i += 1
            ##downstream 
            i = 1
            #print 'DOWNSTREAM'
            while index -i >= 0:
                if not l_by_end[index-i].resize(): #c'est pas juste un test ça va resizé et si jamais y'a pas de start valide ça return False...
                    i += 1 # on passe donc au suivant et on vérifie qu'on soit pas sortie des bordure de la liste
                    continue
                if ListGene.dist_inf <= orf.first_start - l_by_end[index-i].end <= ListGene.dist_max:
                    #print l_by_end[index-i].first_start,'-->', l_by_end[index-i].end
                    #print orf.first_start,' - ',l_by_end[index-i].end,'=', orf.first_start - l_by_end[index-i].end 
                    self.set_adj_orf.add(l_by_end[index-i])
                elif orf.first_start - l_by_end[index-i].end > ListGene.dist_max :
                    break #on peut faire ça car les liste sont sorted :-) 
                i += 1
            
def translate(nucl_seq,codonTable=11):
    """Returns the translated proteic sequence based on the standard genetic code"""

    proteic_seq={}
    codons=getWords(nucl_seq,3)
    
   
    code=getGeneticCode(codonTable)

    nucl_seq["data"]=nucl_seq["data"].upper()
    #proteic_seq["description"]="translated proteic sequence from "+nucl_seq["description"]
    proteic_seq["data"]=""
    for codon in codons:
        if codon not in code:
           proteic_seq["data"]=proteic_seq["data"]+"?"
        else:
            proteic_seq["data"]=proteic_seq["data"]+code[codon]
    
    return proteic_seq

def getWords(seq,wlength):

    """Returns a list of non-overlapping words of a given length in a seq"""

    seq["data"]=seq["data"].upper()
    words=[]
    for i in range(0,len(seq["data"])-wlength+1,wlength):
        words.append(seq["data"][i:i+wlength])

    return words


def getGeneticCode(NCBI_ID):
  """Returns the genetic code corresponding to a given NCBI_ID

  This function is written by Theo Falgarone.

  Args:
      NCBI_ID: NCBI identifier for the genetic code
      
  Returns:
      table: genetic code table
  
  """
  if NCBI_ID==1:
    base1="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    base2="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    base3="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
    start="---M---------------M---------------M----------------------------"
    aas  ="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    table = {}
    codonstart=[]
    codonstop=[]
    for i in range(0,len(aas)):
      codon =(base1[i]+base2[i]+base3[i])
      aa =aas[i]
      table[codon] = aa
      if aa=='*':
        codonstop.append(codon)
      elif start[i]=='M':
        codonstart.append(codon)
    table['stop']=codonstop
    table['start']=codonstartAll
    return table
  elif NCBI_ID==11:
    base1="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    base2="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    base3="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
    start="---M---------------M------------MMMM---------------M------------"
    aas  ="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    table = {}
    codonstart=[]
    codonstop=[]
    for i in range(0,len(aas)):
      codon = (base1[i]+base2[i]+base3[i])
      aa =aas[i]
      table[codon] = aa
      if aa=='*':
        codonstop.append(codon)
      elif start[i]=='M':
        codonstart.append(codon)
    table['stop']=codonstop
    table['start']=codonstart
    return table
  else:
    print('ERROR : votre NCBI_ID ne correspond pas, veuillez entrer 1 ou 11')

#3_2_4 : Afficher le brin complementaire d une sequence adn (3'->5')
def complement(seq):
    result =''
    complement={}
    for i in seq['data']:
        if i=='A':
            result=result+'T'
        elif i=='T':
            result=result+'A'
        elif i=='C':
            result=result+'G'
        elif i=='G':
            result=result+'C'
    complement['data'] =result
    return complement

#3_2_5 : Afficher le brin d une sequence nucleique dans l autre sens (si 3'->5' alors 5'->3')
def reverse(seq):
    rev_seq={}
    rev_seq['data']=seq['data'][::-1]
    return rev_seq

#3_2_6 : Afficher le brin complementaire d une sequence nucleique dans l autre sens (si 3'->5' alors 5'->3')
def complement_reverse(seq):
    compl_rev_seq={}
    compl_rev_seq=reverse(seq)
    compl_rev_seq=complement(compl_rev_seq)
    return compl_rev_seq
