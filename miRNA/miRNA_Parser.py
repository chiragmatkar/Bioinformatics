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
# Construct miRNA Data pool from  miRNA.dat,mature.fa and miFam.dat extrracted from miRBase for  python bindings
# miRNA ids and sequences are now Biopython Sequence objects and can be used in various other programs to find complements,ORFS etc 
#
#
from Bio import SeqIO
import urllib
import gzip


def extract(urla,zipa,filea):
    print 'Downloading.. '+ urla
    handle = urllib.urlopen(urla)
    with open(zipa, 'wb') as out:
        while True:
            data = handle.read(1024)
            if len(data) == 0: break
            out.write(data)
	handle = gzip.open(zipa)
    with open(filea, 'w') as out:
        for line in handle:
            out.write(line)


def miRNA():
	dict={}
	for seq_record in SeqIO.parse('mature.fa','fasta'):
		dict[seq_record.id]=seq_record.seq
	return dict

def organism_miRNA(name,mirna={}):
	dict={}
	for key, value in mirna.items():
		if name in key: 
			dict[key]=value
	return dict


#Extract miRNA databases from mirbase.org
#extract('ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz','mature.fa.gz','mature.fa')
#extract('ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz','miRNA.dat.gz','miRNA.dat')
#extract('ftp://mirbase.org/pub/mirbase/CURRENT/miFam.dat.gz','miFam.dat.gz','miFam.dat')
	
#Whole miRNA dict
mirna=miRNA()
print len(mirna)
print mirna['cel-let-7-5p']
print mirna['rno-miR-1839-3p']

#Specified organism miRNA dict example Homo Sapiens
hsa=organism_miRNA('hsa',mirna)
print len(hsa)
print hsa['hsa-miR-153-5p']

#Specified organism miRNA dict example Drosophila melanogaster
dme=organism_miRNA('dme',mirna)
print len(dme)
print dme['dme-miR-4983-3p']
print len(dme['dme-miR-4983-3p'])

#find complement from seq object
print dme['dme-miR-4983-3p'].complement()
