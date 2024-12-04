import sys
import pysam
import gzip
import io


varlist=sys.argv[1]
bamList=sys.argv[2]
outf=sys.argv[3]


#Parse VCF lines
def getPositionInfo(line):
	global chr, pos, ref, alt, vartype, indelLength, mnvLength
	if line.split('\t')[0][0:3]=='chr':
		chr=line.split('\t')[0]
	else:
		chr=str(line.split('\t')[0])
	pos=int(line.split('\t')[1])
	ref=line.split('\t')[2]
	alt=line.rstrip('\n').split('\t')[3]
	if len(ref)==1 and len(alt)==1:
		vartype="SNV"
	elif len(ref)==1 and len(alt)>1:
		vartype="INS"
	elif len(ref)>1 and len(alt)==1:
		vartype="DEL"
	elif len(ref)>1 and len(alt)==len(ref):
		vartype="MNV"
		mnvLength=len(ref)
	if vartype=="INS" or vartype=="DEL":
		indelLength=len(alt)-len(ref)



### SNVs 

#Filtering functions for count_coverage

def f1r2(read):
	if read.is_read1 and not read.is_reverse:
		return True
	elif read.is_read2 and read.is_reverse:
		return True
	else:
		return False

def f2r1(read):
	if read.is_read2 and not read.is_reverse:
		return True
	elif read.is_read1 and read.is_reverse:
		return True
	else:
		return False

def revread(read):
	if read.is_reverse:
		return True
	else:
		return False

def forwread(read):
	if not read.is_reverse:
		return True
	else:
		return False



#Functions to get counts per allele and diff filters

nucleotides=["A","C","G","T"]

#Depth per allele - min BQ 20
def get_AD_counts(chr,pos,bamfile):
	aCounts=[]
	for i in bamfile.count_coverage(chr,pos-1,pos,quality_threshold=20):
		aCounts.append("".join(str(x) for x in i))
	dictaCounts=dict(zip(nucleotides,aCounts))
	return dictaCounts

#Depth per allele - forward read1 or reverse read2
def get_F1R2_counts(chr,pos,bamfile):
	f1r2_counts=[]
	for i in bamfile.count_coverage(chr,pos-1,pos,quality_threshold=20,read_callback=f1r2):
		f1r2_counts.append("".join(str(x) for x in i))
	dict_f1r2=dict(zip(nucleotides,f1r2_counts))
	return dict_f1r2

#Depth per allele - forward read2 or reverse read1
def get_F2R1_counts(chr,pos,bamfile):
	f2r1_counts=[]
	for i in bamfile.count_coverage(chr,pos-1,pos,quality_threshold=20,read_callback=f2r1):
		f2r1_counts.append("".join(str(x) for x in i))
	dict_f2r1=dict(zip(nucleotides,f2r1_counts))
	return dict_f2r1

#Depth per allele - reverse strand
def get_rev_counts(chr,pos,bamfile):
	revCounts=[]
	for i in bamfile.count_coverage(chr,pos-1,pos,quality_threshold=20,read_callback=revread):
		revCounts.append("".join(str(x) for x in i))
	dictrevCounts=dict(zip(nucleotides,revCounts))
	return dictrevCounts

#Depth per allele - forward strand
def get_forward_counts(chr,pos,bamfile):
	forwCounts=[]
	for i in bamfile.count_coverage(chr,pos-1,pos,quality_threshold=20,read_callback=forwread):
		forwCounts.append("".join(str(x) for x in i))
	dictforwCounts=dict(zip(nucleotides,forwCounts))
	return dictforwCounts

#Generate sample annotation FOR SNVs
def getAnnotation(chr,pos,bamfile):
	#Get counts
	dictaCounts=get_AD_counts(chr,pos,bamfile)
	dict_f1r2=get_F1R2_counts(chr,pos,bamfile)
	dict_f2r1=get_F2R1_counts(chr,pos,bamfile)
	dictrevCounts=get_rev_counts(chr,pos,bamfile)
	dictforwCounts=get_forward_counts(chr,pos,bamfile)
	dictforwCounts=get_forward_counts(chr,pos,bamfile)
	#Calculate VCF annotations
	AD=",".join([dictaCounts[ref],dictaCounts[alt]])
	DP=sum([int(x) for x in list(dictaCounts.values())])
	if DP==0:
		AF=0
	else:
		AF=int(dictaCounts[alt])/DP
	F1R2=",".join([dict_f1r2[ref],dict_f1r2[alt]])
	F2R1=",".join([dict_f2r1[ref],dict_f2r1[alt]])
	SB=",".join([dictforwCounts[ref],dictrevCounts[ref],dictforwCounts[alt],dictrevCounts[alt]])
	#Put together in one string
	info=":".join([str(x) for x in ["0/0",AD,round(AF,3),DP,F1R2,F2R1,SB]])
	return info



## INDELS 

def get_indel_info(chr,pos,bamfile,indelLength):
	for pileupcolumn in bamfile.pileup(chr, pos, pos+1, ignore_overlaps=False):
		depth=0
		depthPlus=0
		depthMinus=0
		nindelPlus=0
		nindelMinus=0
		depthF1R2=0
		depthF2R1=0
		nindelF1R2=0
		nindelF2R1=0
		#Get all reads aligned to position
		if pileupcolumn.pos==pos-1:
			#Iterate over the reads
			for pileupread in pileupcolumn.pileups:
				depth+=1
				#Reverse strand
				if pileupread.alignment.is_reverse:
					depthMinus+=1
					if pileupread.indel==indelLength:
						nindelMinus+=1
					#R1
					if pileupread.alignment.is_read1:
						depthF2R1+=1
						if pileupread.indel==indelLength:
							nindelF2R1+=1
					#R2
					elif pileupread.alignment.is_read2:
						depthF1R2+=1
						if pileupread.indel==indelLength:
							nindelF1R2+=1
				#Forward strand (there's no function to test forwardness)
				else:
					depthPlus+=1
					if pileupread.indel==indelLength:
						nindelPlus+=1
					#F1
					if pileupread.alignment.is_read1:
						depthF1R2+=1
						if pileupread.indel==indelLength:
							nindelF1R2+=1
					#F2
					elif pileupread.alignment.is_read2:
						depthF2R1+=1
						if pileupread.indel==indelLength:
							nindelF2R1+=1
			DP=depth
			altAD=nindelPlus+nindelMinus
			refAD=DP-altAD
			AD=",".join([str(x) for x in [refAD,altAD]])
			if DP==0:
				AF=0
			else:
				AF=round(altAD/DP,3)
			F1R2=",".join([str(x) for x in [depthF1R2-nindelF1R2,nindelF1R2]])
			F2R1=",".join([str(x) for x in [depthF2R1-nindelF2R1,nindelF2R1]])
			SB=",".join([str(x) for x in [depthPlus-nindelPlus,depthMinus-nindelMinus,nindelPlus,nindelMinus]])
			info=":".join([str(x) for x in ["0/0",AD,AF,DP,F1R2,F2R1,SB]])
			return info
	



### PROCESS FILE

#Read BAM list
bams=[line.rstrip('\n') for line in open(bamList)]
	
out=open(outf, 'a')

for line in open(varlist):
	if line.split("\t")[3].count(",")>0:
		sys.exit("The input VCF file contains multiallelic entries. Please normalise it.")
	else:
		#Makes chr, pos, ref, alt, vartype, indelLength, mnvLength - global variables
		getPositionInfo(line)
		lbroken=line.strip().split("\t")
		#Loop over each sample to add
		for bam in bams:
			#Import bam
			bamfile=pysam.AlignmentFile(bam,"rb")
			if vartype=="SNV":			
				info=getAnnotation(chr,pos,bamfile)
			elif vartype=="INS" or vartype=="DEL":
				info=get_indel_info(chr,pos,bamfile,indelLength)
			elif vartype=="MNV":
				info_MNV=[]
				for i in range(mnvLength):
					ref=line.split('\t')[2][i]
					alt=line.split('\t')[3][i]
					info=getAnnotation(chr,pos+i,bamfile)
					info_MNV.append(info)
				AD_MNV=[]
				for info in info_MNV:
					AD_MNV.append(info.split(":")[1].split(",")[1])
				info=info_MNV[AD_MNV.index(min(AD_MNV))]
			#If coverage==0, fill in info which will be empty
			if not(info):
				info="0/0:0,0:0.0:0:0,0:0,0:0,0,0,0"
			#Append to lbroken
			lbroken.append(info)
		#Write new line
		out.write('\t'.join(lbroken)+'\n')


out.close()






