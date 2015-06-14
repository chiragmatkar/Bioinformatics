#
# Copyright (c) 2015 Chirag Matkar
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Chirag Matkar <chirag.matkar@gmail.com>
#
# Python wrapper to create modeller Homology modelling files 
# and evaluate the model at one go instead of several intermediate steps
# Customization of Homology modelling can be done as required

from modeller import *
from modeller.scripts import complete_pdb
from modeller.automodel import *
import urllib
#from Bio import SeqIO
#from modeller import soap_protein_od
log.verbose()
env = environ()


def sequence_alignment():
	pass

def profile_build(filename,db_file):
	sdb = sequence_db(env)
	sdb.read(seq_database_file= db_file+'.pir', seq_database_format='PIR',
         chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

	#-- Write the sequence database in binary form
	sdb.write(seq_database_file=db_file+'.bin', seq_database_format='BINARY',
          chains_list='ALL')

	#-- Now, read in the binary database
	sdb.read(seq_database_file= db_file+'.bin', seq_database_format='BINARY',
         chains_list='ALL')

	#-- Read in the target sequence/alignment
	aln = alignment(env)
	aln.append(file=filename, alignment_format='PIR', align_codes='ALL')

	#-- Convert the input sequence/alignment into
	#   profile format
	prf = aln.to_profile()

	#-- Scan sequence database to pick up homologous sequences
	prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

	#-- Write out the profile in text format
	prf.write(file='build_profile.prf', profile_format='TEXT')

	#-- Convert the profile back to alignment format
	aln = prf.to_alignment()

	#-- Write out the alignment file
	aln.write(file='build_profile.ali', alignment_format='PIR')



			  
def compare():
	aln = alignment(env)
	for (pdb, chain) in (('1b8p', 'A'), ('1bdm', 'A'), ('1civ', 'A'),
                     ('5mdh', 'A'), ('7mdh', 'A'), ('1smk', 'A')):
            m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
            aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
	    aln.malign()
	    aln.malign3d()
	    aln.compare_structures()
	    aln.id_table(matrix_file='family.mat')
	    env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)

	
def align_2d():
	aln = alignment(env)
	mdl = model(env, file='1bdm', model_segment=('FIRST:A','LAST:A'))
	aln.append_model(mdl, align_codes='1bdmA', atom_files='1bdm.pdb')
	aln.append(file='TvLDH.ali', align_codes='TvLDH')
	aln.align2d()
	aln.write(file='TvLDH-1bdmA.ali', alignment_format='PIR')
	aln.write(file='TvLDH-1bdmA.pap', alignment_format='PAP')

	
	
def build_single_model(aln,known,seq,total_models):
	a = automodel(env, alnfile=aln,
				knowns=known, sequence=seq,
				assess_methods=(assess.DOPE,assess.GA341))
	a.starting_model = 1
	a.ending_model = int(total_models)
	a.make()
			  
			  
			  
def r_enumerate(seq):
    """Enumerate a sequence in reverse order"""
    # Note that we don't use reversed() since Python 2.3 doesn't have it
    num = len(seq) - 1
    while num >= 0:
        yield num, seq[num]
        num -= 1

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    f = file(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in r_enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals			  
			  
	
def evaluate_single_model(pdb):
	env.libs.topology.read(file='$(LIB)/top_heav.lib') 
	env.libs.parameters.read(file='$(LIB)/par.lib') 
	mdl = complete_pdb(env, pdb)
	s = selection(mdl)   
	s.assess_dope(output='ENERGY_PROFILE NO_REPORT', 
					file='TvLDH.profile',
              normalize_profile=True, smoothing_window=15)


def get_PDB(file):
	url = 'http://www.rcsb.org/pdb/files/'+file
	print 'Downloading..... '+ url
        handle = urllib.urlopen(url)
        with open(file, 'wb') as out:
            while True:
                data = handle.read(1024)
                if len(data) == 0: break
                out.write(data)

				
				
				
#Modeller Steps


#Step 1 : create target sequence file in PIR format.If necessary convert FASTA to PIR	input.txt 	  


#Step 2 : create a sequence alignment file .PIR containing non-redundant PDB sequences at 95% sequence identity into the sdb database
#sequence_alignment()

#Step 3 : Search for similar structures in sequence_alignment file using profile.build 
#profile_build('Human.ali','pdb_human')
profile_build('TvLDH.ali','pdb_95')

#Step 4 : Search most  similar template from set of structures using compare
#compare()


#Step 5 : Download PDB files structures
#get_PDB('1H1B.pdb')


#Step 6: Align sequence with Structure 


#Step 7:Build Model
#build_single_model('TvLDH-1bdmA.ali','1bdmA','TvLDH',5)


#Step 8:Validate Model
#evaluate_single_model'TvLDH.B99990002.pdb')
