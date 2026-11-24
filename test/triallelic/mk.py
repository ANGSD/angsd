#makes triallelic/sites.glf: chr pos and 10 log likelihoods (AA AC AG AT CC CG CT GG GT TT) for 10 individuals.
#The true genotype has likelihood 0, genotypes that share an allele with it 'a', the others 'b'
G=['AA','AC','AG','AT','CC','CG','CT','GG','GT','TT']
def ind(g,a,b):
    return ['0' if x==g else (str(a) if set(x)&set(g) else str(b)) for x in G]
sites=[(1,['AA']*6+['AC']*4,-0.5,-1),#biallelic A/C
       (2,['GG']*6+['GT']*4,-0.5,-1),#biallelic G/T
       (3,['CC']*4+['CG']*3+['CT']*3,-5,-10)]#triallelic C/G/T
for p,gs,a,b in sites:
    print('\t'.join(['1',str(p)]+sum([ind(g,a,b) for g in gs],[])))
