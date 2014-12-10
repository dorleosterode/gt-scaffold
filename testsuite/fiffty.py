import sys

if len(sys.argv) < 4:
    sys.stderr.write("Usage: fiffty.py <fasta> <libPE.de> <libPe.astat> [<n> <m>]\n")
    sys.stderr.write("keeps the first <n> contigs in <fasta> that are longer than <m>\n")
    sys.stderr.write("filters the files <libPE.de> and <libPE.asta> accordingly\n")
    sys.stderr.write("n: default 50\n")
    sys.stderr.write("m: default 200\n")
    sys.stderr.write("the new files are named new.fa, new.de and new.astat\n")
    exit(1)

fasta = sys.argv[1]
de = sys.argv[2]
astat = sys.argv[3]

n = 50
m = 200

if len(sys.argv) >= 5:
    n = int(sys.argv[4])
if len(sys.argv) == 6:
    m = int(sys.argv[5])

new_fasta = "new.fa"
new_de = "new.de"
new_astat = "new.astat"

contigs = {}
count = 0

with open(fasta) as f, open(new_fasta, 'w') as out:
    get_entry = False
    for line in f:
        if line.startswith('>'):
            get_entry = False
            if int(line.split(' ')[1]) > m:
                get_entry = True
                count += 1
            else:
                continue
        if count == n + 1:
            break
        if get_entry:
            name = line.split(' ')[0][1:]
            contigs[name] = True
            out.write(line)

with open(de) as d, open(new_de, 'w') as out:
    for line in d:
        new_list = []
        lst = line.split(' ')
        try:
            b = contigs[lst[0]]
            new_list.append(lst[0])
        except KeyError:
            continue

        for elem in lst[1:]:
            try:
                b = contigs[elem.split(',')[0][:-1]]
                new_list.append(elem)
            except KeyError:
                if elem == ';':
                    new_list.append(elem)

        out.write(' '.join(new_list))
        if not new_list[-1].endswith('\n'):
            out.write('\n')

with open(astat) as a, open(new_astat, 'w') as out:
    for line in a:
        contig = line.split('\t')[0]
        try:
            b = contigs[contig]
            out.write(line)
        except KeyError:
            pass
