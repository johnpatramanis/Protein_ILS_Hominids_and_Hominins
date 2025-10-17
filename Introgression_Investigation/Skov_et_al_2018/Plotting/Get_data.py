from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Find archaic segments overlapping regions')
parser.add_argument("-chrom",help="chromosome position", type=str)
parser.add_argument("-start",help="Start position", type=int)
parser.add_argument("-end",help="End position", type=int)
parser.add_argument("-outfile",help="Outfile", type=str)
parser.add_argument("-windowsize",help="Windowsize", type=int, default = 50000)

args = parser.parse_args()

GENE_START = args.start - args.windowsize
GENE_END = args.end + args.windowsize
CHROM = args.chrom


row_number_dict = defaultdict(int)
counter = 0

with open(args.outfile + '.segments', 'w') as out, open(args.outfile + '.snps', 'w') as out_snps:

	# write headers
	print('rownumber', 'individual_name', 'chrom', 'start', 'end', 'Archaic', 'region','population', sep = '\t', file = out)
	print('rownumber', 'chrom', 'pos', 'snptype', 'individual_name', 'region', 'population', 'Archaic', sep = '\t', file = out_snps)

	for dataset in ['HGDP', '1000g']:

		# Overlap segments with snps
		snptypes = defaultdict(lambda: 'unlinked')    
		with open(f'hg38_{dataset}_SNPS.txt') as data:
			for line in data:
				if not line.startswith('chrom'):
					chrom, pos, snptype, ancestral, derived, freq, NDtype, matches = line.strip().split()[0:9]
					ID = f'{chrom}_{pos}'
					snptypes[ID] = snptype

					if snptype == 'DAV':    
						snptypes[ID] = NDtype
					

	
			
		with open(f'../hg38_{dataset}_segments.txt') as data:
			for line in data:
				if not line.startswith('name'):

					name, haplotype, population, region, chrom, start, end, mean_prob, ND_type, snps, admixpopvariants, Altai, Vindija, Denisova, Chagyrskaya, variants = line.strip().split('\t')



					start, end = int(start), int(end)
					Altai, Vindija, Denisova, Chagyrskaya = int(Altai), int(Vindija), int(Denisova), int(Chagyrskaya)


					if chrom == CHROM and start < GENE_END and end > GENE_START:

						if float(mean_prob) > 0.9:

							ID = f'{name}_{haplotype}'
							if ID not in row_number_dict:
								counter += 1
								row_number_dict[ID] = counter 
								
							rownumber = row_number_dict[ID]


							print(rownumber, name, chrom, start, end, ND_type, region, population, sep = '\t', file = out)
							print(rownumber, name, chrom, start, end, ND_type, region, population, sep = '\t')

							for snp in variants.split(','):
								print(rownumber, chrom, snp, snptypes[f'{chrom}_{snp}'], name, region, population, ND_type, sep = '\t', file = out_snps)

