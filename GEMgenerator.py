#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 19:28:33 2022

@author: omidard
"""
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import os
from os.path import join
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
from os.path import join
from os.path import isfile, join
from os import listdir



#generate data directories
directories = ['gap_inf_dir','reference_genome_dir', 'target_genome_dir', 'prots_dir', 'nucls_dir', 'bbh_dir' ,'present_absence_dir','initial_models_dir', 'output_models_dir', 'gapfilled_models_dir', 'blast_exe_dir','temp_files','ref_model_dir','backupgenome']
for i in directories:
    cmd_line = 'mkdir '+i
    os.system(cmd_line)


def get_file_names (directory1,directory2):
    onlyfiles = [f for f in listdir(directory1) if isfile(join(directory1, f))]
    gl8 = []
    for i in onlyfiles:
        gl1 = i.replace('.gbk','')
        gl8.append(gl1)
    gl9=pd.DataFrame({'strain':gl8,'NCBI ID':gl8, 'Pathotype':range(len(onlyfiles))})
    path = directory2+'/StrainInformation.xlsx'
    StrainInformation = gl9.to_excel(path)


def dl_genome(id, folder='genomes'): # be sure get CORRECT ID
    files=glob('%s/*.gbk'%folder)
    out_file = '%s/%s.gbk'%(folder, id)

    if out_file in files:
        print (out_file, 'already downloaded')
        return
    else:
        print ('downloading %s from NCBI'%id)
        
    from Bio import Entrez
    Entrez.email = "omidard@biosustain.dtu.dk"     #Insert email here for NCBI
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    fout = open(out_file,'w')
    fout.write(handle.read())
    fout.close()


def get_strain_info(directory1):
    files = glob('%s/*.gbk'%directory1)
    strain_info = []
    
    for file in files:
        handle = open(file)
        record = SeqIO.read(handle, "genbank")
        for f in record.features:
            if f.type=='source':
                info = {}
                info['file'] = file
                info['id'] = file.split('\\')[-1].split('.')[0]
                for q in f.qualifiers.keys():
                    info[q] = '|'.join(f.qualifiers[q])
                strain_info.append(info)
    return pd.DataFrame(strain_info)


def parse_genome(id, type='prot', in_folder='genomes', out_folder='prots', overwrite=1):

    in_file = '%s/%s.gbk'%(in_folder, id)
    out_file='%s/%s.fa'%(out_folder, id)
    files =glob('%s/*.fa'%out_folder)
    
    if out_file in files and overwrite==0:
        print (out_file, 'already parsed')
        return
    else:
        print ('parsing %s'%id)
    
    handle = open(in_file)
    
    fout = open(out_file,'w')
    x = 0
    
    records = SeqIO.parse(handle, "genbank")
    for record in records:
        for f in record.features:
            if f.type=='CDS':
                seq=f.extract(record.seq)
                
                if type=='nucl':
                    seq=str(seq)
                else:
                    seq=str(seq.translate())
                    
                if 'locus_tag' in f.qualifiers.keys():
                    locus = f.qualifiers['locus_tag'][0]
                elif 'gene' in f.qualifiers.keys():
                    locus = f.qualifiers['gene'][0]
                else:
                    locus = 'gene_%i'%x
                    x+=1
                fout.write('>%s\n%s\n'%(locus, seq))
    fout.close()


def make_blast_db(id,folder='prots',db_type='prot'):
    import os
    
    out_file ='%s/%s.fa.pin'%(folder, id)
    files =glob('%s/*.fa.pin'%folder)
    
    if out_file in files:
        print (id, 'already has a blast db')
        return
    if db_type=='nucl':
        ext='fna'
    else:
        ext='fa'

    cmd_line="makeblastdb -in %s/%s.%s -dbtype %s" %(folder, id, ext, db_type)
    
    print (' making blast db with following command line...')
    print (cmd_line)
    os.system(cmd_line)

def run_blastp(seq,db,in_folder='prots_dir', out_folder='bbh_dir', out=None,outfmt=6,evalue=0.001,threads=64):
    import os
    if out==None:
        out='%s/%s_vs_%s.txt'%(out_folder, seq, db)
        print(out)
    
    files =glob('%s/*.txt'%out_folder)
    if out in files:
        print (seq, 'already blasted')
        return
    
    print ('blasting %s vs %s'%(seq, db))
    
    db = '%s/%s.fa'%(in_folder, db)
    seq = '%s/%s.fa'%(in_folder, seq)
    cmd_line='blastp -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastp with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out


def get_gene_lens(query, in_folder='prots'):

    file = '%s/%s.fa'%(in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out


def get_bbh(query, subject, in_folder='bbh'):    
    
    #Utilize the defined protein BLAST function
    run_blastp(query, subject)
    run_blastp(subject, query)
    
    query_lengths = get_gene_lens(query, in_folder='prots_dir')
    subject_lengths = get_gene_lens(subject, in_folder='prots_dir')
    
    #Define the output file of this BLAST
    out_file = '%s/%s_vs_%s_parsed.csv'%(in_folder,query, subject)
    files=glob('%s/*_parsed.csv'%in_folder)
    
    #Combine the results of the protein BLAST into a dataframe
    print ('parsing BBHs for', query, subject)
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    bbh=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,query, subject), sep='\t', names=cols)
    bbh = pd.merge(bbh, query_lengths) 
    bbh['COV'] = bbh['alnLength']/bbh['gene_length']
    
    bbh2=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,subject, query), sep='\t', names=cols)
    bbh2 = pd.merge(bbh2, subject_lengths) 
    bbh2['COV'] = bbh2['alnLength']/bbh2['gene_length']
    out = pd.DataFrame()
    
    # Filter the genes based on coverage
    bbh = bbh[bbh.COV>=0.30]
    bbh2 = bbh2[bbh2.COV>=0.30]
    
    #Delineate the best hits from the BLAST
    for g in bbh.gene.unique():
        res = bbh[bbh.gene==g]
        if len(res)==0:
            continue
        best_hit = res.loc[res.PID.idxmax()]
        best_gene = best_hit.subject
        res2 = bbh2[bbh2.gene==best_gene]
        if len(res2)==0:
            continue
        best_hit2 = res2.loc[res2.PID.idxmax()]
        best_gene2 = best_hit2.subject
        if g==best_gene2:
            best_hit['BBH'] = '<=>'
        else:
            best_hit['BBH'] = '->'
        out=pd.concat([out, pd.DataFrame(best_hit).transpose()])
    
    #Save the final file to a designated CSV file
    out.to_csv(out_file)


def gbk2fasta(gbk_filename):
    faa_filename = '.'.join(gbk_filename.split('.')[:-1])+'.fna'
    input_handle  = open(gbk_filename, "r")
    output_handle = open(faa_filename, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank") :
        print ("Converting GenBank record %s" % seq_record.id)
        output_handle.write(">%s %s\n%s\n" % (
               seq_record.id,
               seq_record.description,
               seq_record.seq))

    output_handle.close()
    input_handle.close()


def run_blastn(seq, db,outfmt=6,evalue=0.001,threads=64):
    import os
    out = 'nucls_dir/'+seq+'_vs_'+db+'.txt'
    seq = 'nucls_dir/'+seq+'.fa'
    db = 'target_genome_dir/'+db+'.fna'
    
    cmd_line='blastn -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastn with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out



def parse_nucl_blast(infile):
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    data = pd.read_csv(infile, sep='\t', names=cols)
    data = data[(data['PID']>50) & (data['alnLength']>0.5*data['queryEnd'])]
    data2=data.groupby('gene').first()
    return data2.reset_index()


def extract_seq(g, contig, start, end):
    from Bio import SeqIO
    handle = open(g)
    records = SeqIO.parse(handle, "fasta")
    
    for record in records:
        if record.name==contig:
            if end>start:
                section = record[start:end]
            else:
                section = record[end-1:start+1].reverse_complement()
                
            seq = str(section.seq)
    return seq





get_file_names('backupgenome','temp_files')
StrainsOfInterest=pd.read_excel('temp_files'+'/StrainInformation.xlsx')
print(StrainsOfInterest)
referenceStrainID='reactome'
targetStrainIDs=list(StrainsOfInterest['NCBI ID'])



directory1 ='target_genome_dir'
files = glob('%s/*.gbk'%directory1)
strain_info = []
for file in files:
        handle = open(file)
        record = list(SeqIO.parse(handle, "genbank"))
        for i in record:
            for f in i.features:
                if f.type=='source':
                    info = {}
                    info['file'] = file.replace('target_genome_dir','')
                    info['id'] = file.replace('target_genome_dir','')
                    info['genome_size'] = len(i.seq)/1000000
                    for q in f.qualifiers.keys():
                        info[q] = '|'.join(f.qualifiers[q])
                        strain_info.append(info)
sinf = pd.DataFrame(strain_info)
                        
                        
sinf2 = sinf.drop_duplicates(subset='id', keep='first', inplace=False, ignore_index=False)
print(sinf2)




for strain in targetStrainIDs:
    parse_genome(strain, type='prot', in_folder='target_genome_dir', out_folder='prots_dir')
    parse_genome(strain, type='nucl', in_folder='target_genome_dir', out_folder='nucls_dir')
    


for strain in targetStrainIDs:
    make_blast_db(strain,folder='prots_dir',db_type='prot')
make_blast_db(referenceStrainID,folder='prots_dir',db_type='prot')


for strain in targetStrainIDs:
    get_bbh(referenceStrainID,strain, in_folder='bbh_dir')


blast_files=glob('%s/*_parsed.csv'%'bbh_dir')
for blast in blast_files:
    bbh=pd.read_csv(blast)
    print (blast,bbh.shape)


data_dir = 'ref_model_dir'
model = cobra.io.read_sbml_model(join(data_dir, "marlbr2.xml"))
listGeneIDs=[]
for gene in model.genes:
    listGeneIDs.append(gene.id)


ortho_matrix=pd.DataFrame(index=listGeneIDs,columns=targetStrainIDs)
geneIDs_matrix=pd.DataFrame(index=listGeneIDs,columns=targetStrainIDs)
print(len(listGeneIDs))


for blast in blast_files:
    bbh=pd.read_csv(blast)
    listIDs=[]
    listPID=[]
    for r,row in ortho_matrix.iterrows():
        try:
            currentOrtholog=bbh[bbh['gene']==r].reset_index()
            listIDs.append(currentOrtholog.iloc[0]['subject'])
            listPID.append(currentOrtholog.iloc[0]['PID'])
        except:
            listIDs.append('None')
            listPID.append(0)
    for col in ortho_matrix.columns:
        if col in blast:
            ortho_matrix[col]=listPID
            geneIDs_matrix[col]=listIDs
print(sorted(ortho_matrix))



for column in ortho_matrix:
    ortho_matrix.loc[ortho_matrix[column]<=50.0,column]=0
    ortho_matrix.loc[ortho_matrix[column]>50.0,column]=1


for strain in targetStrainIDs:
    gbk2fasta('target_genome_dir/'+strain+'.gbk')



for strain in targetStrainIDs:
    make_blast_db(strain,folder='target_genome_dir',db_type='nucl')


genome_blast_res=[]
for strain in targetStrainIDs:
    res = run_blastn(referenceStrainID,strain)
    genome_blast_res.append(res)



na_matrix=pd.DataFrame()
for file in genome_blast_res:
    genes =parse_nucl_blast(file)
    name ='.'.join(file.split('_')[-1].split('.')[:-1])
    na_matrix = na_matrix.append(genes[['gene','subject','PID']])
na_matrix = pd.pivot_table(na_matrix, index='gene', columns='subject',values='PID')



ortho_matrix_w_unannotated = ortho_matrix.copy()
geneIDs_matrix_w_unannotated = geneIDs_matrix.copy()


nonModelGenes=[]
for g in na_matrix.index:
    if g not in listGeneIDs:
        nonModelGenes.append(g)

na_model_genes=na_matrix.drop(nonModelGenes)



pseudogenes = {}

for c in ortho_matrix.columns:
    
    orfs = ortho_matrix[c]
    genes = na_model_genes
    # All the Model Genes that met the BLASTp Requirements
    orfs2 = orfs[orfs==1].index.tolist()
    # All the Model Genes based off of BLASTn similarity above threshold of 80
    genes2 = genes[genes>=50].index.tolist()
    # By Definition find the genes that pass sequence threshold but were NOT in annotated ORFs:
    unannotated = set(genes2) -set(orfs2)
    
    # Obtain sequences of this list to check for premature stop codons:
    data = 'nucls_dir/reactome_vs_%s.txt'%c
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    data = pd.read_csv(data, sep='\t', names=cols)
    #
    pseudogenes[c] = {}
    unannotated_data = data[data['gene'].isin(list(unannotated))]
    for i in unannotated_data.index:
        gene = data.loc[i,'gene']
        contig = data.loc[i,'subject']
        start = data.loc[i,'subjectStart']
        end = data.loc[i,'subjectEnd']
        seq = extract_seq('target_genome_dir/%s.fna'%c,contig, start, end)
        # check for early stop codons - these are likely nonfunctional and shouldn't be included
        if '*' in seq:
            print (seq)
            pseudogenes[c][gene]=seq
            # Remove the gene from list of unannotated genes
            unannotated-set([gene])
            
    
    print (c, unannotated)
    
    # For pertinent genes, retain those based off of nucleotide similarity within the orthology matrix and geneIDs matrix
    ortho_matrix_w_unannotated.loc[unannotated,c]=1
    for g in unannotated:
        geneIDs_matrix_w_unannotated.loc[g,c] = '%s_ortholog'%g


#Save the Presence/Absence Matrix and geneIDs Matrix for future use
ortho_matrix_w_unannotated.to_csv('present_absence_dir/ortho_matrix.csv')
geneIDs_matrix_w_unannotated.to_csv('present_absence_dir/geneIDs_matrix.csv')

data_dir = 'ref_model_dir'
model = cobra.io.read_sbml_model(join(data_dir, "marlbr2.xml"))


hom_matrix=pd.read_csv('present_absence_dir/ortho_matrix.csv')
hom_matrix=hom_matrix.set_index('Unnamed: 0')



for strain in hom_matrix.columns:
    
    #Get the list of Gene IDs from the homology matrix dataframe for the current strain without a homolog
    currentStrain=hom_matrix[strain]
    nonHomologous=currentStrain[currentStrain==0.0]
    nonHomologous=nonHomologous.index.tolist()
    nonHomologous.remove('spontaneous')
    nonHomologous.remove('EXCHANGE')
    nonHomologous.remove('BIOMASS')
    nonHomologous.remove('Diffusion')
    nonHomologous.remove('GAP')
    nonHomologous.remove('DEMAND')
    nonHomologous.remove('SINK')
    nonHomologous.remove('ORPHAN')
    print(len(nonHomologous))
    #above genes are artificial genes used in reactome for spontaneous/exchange/gap/biomass/demand/diffusion reactions and as such has no homologs,
    #However, it is retained for these spontaneous reactions to function
    #nonHomologous.remove('s0001')
    #Define a list of Gene objects from the base reconstruction to be deleted from the current strain
    model = cobra.io.load_matlab_model(os.path.join('ref_model_dir/marlbr2.mat'))
    toDelete=[]
    for gene in nonHomologous:
        toDelete.append(model.genes.get_by_id(gene))

    #Establish a model copy and use the COBRApy function to remove the appropriate content and save this model
    modelCopy=model.copy()
    remove_genes(modelCopy, toDelete, remove_reactions=True)
    modelCopy.id=str(strain)
    cobra.io.save_json_model(modelCopy, str('initial_models_dir/'+strain+'.json'), pretty=False)
    
    
models=glob('%s/*.json'%'initial_models_dir')
geneIDs_matrix=pd.read_csv('present_absence_dir/geneIDs_matrix.csv')
geneIDs_matrix=geneIDs_matrix.set_index('Unnamed: 0')


from cobra.manipulation.modify import rename_genes
for mod in models:
    model=cobra.io.load_json_model(mod)
    for column in geneIDs_matrix.columns:
        if column in mod:
            currentStrain=column
    
    IDMapping=geneIDs_matrix[currentStrain].to_dict()
    IDMappingParsed = {k:v for k,v in IDMapping.items() if v != 'None'}
    #added
    """
    temp = []
    res = dict()
    for key, val in IDMappingParsed.items():
        if val not in temp:
            temp.append(val)
            res[key] = val
    res2 = dict()
    for key, val in res.items():
        if 'ortholog' not in val:
            res2[key] = val
    ###end --> res2 to IDMappingParsed
    """
    rename_genes(model,IDMappingParsed)
    cobra.io.save_json_model(model,mod, pretty=False)


modelgene = []
modelid = []
genumb = []
reanumb = []
for strain in hom_matrix.columns:
    model=cobra.io.load_json_model(str('initial_models_dir/'+strain+'.json'))
    print (model.id,'Number of Model Genes:',len(model.genes),'Number of Model Reactions:',len(model.reactions))
    modelgene.append(model.genes)
    modelid.append(model.id)
    genumb.append(len(model.genes))
    reanumb.append(len(model.reactions))
    
modelsinf = pd.DataFrame({'model_id':modelid,'reactions':reanumb,'genes':genumb})



models=glob('%s/*.json'%'initial_models_dir')
for i in range(len(models)):
    model=cobra.io.load_json_model(models[i])
    strain = models[i].replace('initial_models_dir/','')
    ort = []
    for ge in model.genes:
        if 'ortholog'in str(ge):
            ort.append(ge)
    modelCopy=model.copy()
    remove_genes(modelCopy, ort, remove_reactions=True)
    modelCopy.id=str(strain)
    cobra.io.json.save_json_model(modelCopy, str('output_models_dir/'+strain+'.json'),pretty=False)



print('\n','------------------------------------------------')
models2=glob('%s/*.json'%'output_models_dir/')
for i in range(len(models2)):
    model=cobra.io.load_json_model(models2[i])
    print (model.id,'Number of Model Genes:',len(model.genes),'Number of Model Reactions:',len(model.reactions))
    







