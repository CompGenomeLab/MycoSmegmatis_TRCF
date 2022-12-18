import sys
 
distanceFromThreePrime = 6

for line in sys.stdin:
    ll = line.rstrip().split('\t')
    strand = ll[5]
    assert strand == '+' or strand == '-', "given strand " + strand + " is not expected."
    if strand == '+':
        ll[1] = int(ll[2]) - distanceFromThreePrime
        ll[2] = ll[1] + 1
    elif strand == '-':
        ll[1] = int(ll[1]) + distanceFromThreePrime - 1
        ll[2] = ll[1] + 1

    ll[1] = str(ll[1])
    ll[2] = str(ll[2])
    print('\t'.join(ll))


