import sys
import os

def create_contigs_hash(d, n):
    contigs = {}
    count = 0

    # construct the hash
    for line in d:
        lst = line.split(' ')
        root = lst[0]
        try:
            b = contigs[root]
        except KeyError:
            if count <= n:
                contigs[root] = False
                count += 1
            else:
                break

        for record in lst[1:]:
            if record == ';':
                continue
            name = record.split(',')[0][:-1]
            try:
                b = contigs[name]
            except KeyError:
                if count <= n:
                    contigs[name] = False
                    count += 1
                else:
                    break
    return contigs

def filter_fasta(fasta, new_fasta, contigs):
    with open(fasta) as f, open(new_fasta, 'w') as out:
        get_entry = False
        for line in f:
            if line.startswith('>'):
                name = line.split(' ')[0][1:]
                try:
                    b = contigs[name]
                    get_entry = True
                except KeyError:
                    get_entry = False
            if get_entry:
                out.write(line)

def filter_astat(astat, new_astat, contigs):
    with open(astat) as a, open(new_astat, 'w') as out:
        for line in a:
            contig = line.split('\t')[0]
            try:
                b = contigs[contig]
                out.write(line)
            except KeyError:
                pass

def filter_de(d, new_de, contigs):
    with open(new_de, 'w') as out:
        # construct the new file
        for line in d:
            lst = line.split(' ')
            root = lst[0]
            new_line = [root]
            try:
                b = contigs[root]
                if not b:
                    contigs[root] = True
                else:
                    print("format is corrupted")
            except KeyError:
                # we havn't found a wanted contig as root
                continue
            for record in lst[1:]:
                if record == ';':
                    new_line.append(record)
                name = record.split(',')[0][:-1]
                try:
                    b = contigs[name]
                    new_line.append(record)
                except KeyError:
                    continue
                if not new_line[-1].endswith('\n'):
                    new_line[-1] = new_line[-1] + '\n'
                out.write(' '.join(new_line))


def main():
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        sys.stderr.write("Usage: traverse_de.py <fasta> <libPE.de> <libPe.astat> [<n>]\n")
        sys.stderr.write("keeps the first <n> contigs that occur in <libPE.de> and all connections between them\n")
        sys.stderr.write("filters the files <fasta> and <libPE.asta> accordingly\n")
        sys.stderr.write("n: default 50\n")
        sys.stderr.write("the new files are named new_<fasta>.fa, new_<libPE>.de and new_<libPE>.astat\n")
        exit(1)

    # get the filenames
    fasta = sys.argv[1]
    de = sys.argv[2]
    astat = sys.argv[3]

    # get the number of contigs
    if len(sys.argv) == 5:
        n = int(sys.argv[4])
    else:
        n = 50

    # create new filenames
    new_fasta = "new_{0}.fa".format(os.path.basename(fasta).split('.')[0])
    new_astat = "new_{0}.astat".format(os.path.basename(astat).split('.')[0])
    new_de = "new_{0}.de".format(os.path.basename(de).split('.')[0])

    # read the whole distance-file into d
    d = []
    with open(de) as dfile:
        d = dfile.readlines()

    contigs = create_contigs_hash(d, n)
    filter_fasta(fasta, new_fasta, contigs)
    filter_astat(astat, new_astat, contigs)
    filter_de(d, new_de, contigs)

if __name__ == "__main__":
    main()
