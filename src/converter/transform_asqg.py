import sys

if len(sys.argv) != 3:
    sys.stderr.write("Usage: {0} <asqg> <ovlfile out>\n".format(sys.argv[0]))
    sys.exit(1)

asqg = sys.argv[1]
ovlfile = sys.argv[2]

seq_nums = {}
seq_num = 0

with open(asqg) as a, open(ovlfile, 'w') as out:
    for line in a:
        if line.startswith("VT"):
            cid = line.split("\t")[1]
            seq_nums[cid] = seq_num
            seq_num += 1
        elif line.startswith("ED"):
            ed_list = line.split("\t")[1].split(" ")
            contig1 = ed_list[0]
            contig2 = ed_list[1]
            ovs1 = (int) (ed_list[2])
            ove1 = (int) (ed_list[3])
            len1 = (int) (ed_list[4])
            ovs2 = (int) (ed_list[5])
            ove2 = (int) (ed_list[6])
            len2 = (int) (ed_list[7])
            direction = (int) (ed_list[8])
            diff = (int) (ed_list[9])

            # print the ovlfile
            c1 = seq_nums[contig1]
            c2 = seq_nums[contig2]
            length = ove1 - ovs1 + 1
            # suffix of contig1
            if ove1 == (len1 - 1):
                # prefix of contig2
                if ovs2 == 0:
                    if direction != 0:
                        print("format is broken: {0}\n".format(line))
                        exit -1
                    # found suffix-prefix-match
                    out.write("{0} + {1} + {2}\n".format(c1, c2, length))
                else:
                    # suffix of contig2
                    if ove2 == (len2 - 1):
                        if direction != 1:
                            print("format is broken: %s\n".format(line))
                            exit -1
                        # found suffix-suffix-match
                        out.write("{0} + {1} - {2}\n".format(c1, c2, length))
                    else:
                        print("format is broken: {0}\n".format(line))
                        exit -1
            else:
                #prefix of contig1
                if ovs1 == 0:
                    # prefix of contig2
                    if ovs2 == 0:
                        if direction != 1:
                            print("format is broken: {0}\n".format(line))
                            exit -1
                        # found prefix-prefix-match
                        out.write("{0} - {1} + {2}\n".format(c2, c1, length))
                    else:
                        # suffix of contig2
                        if ove2 == (len2 - 1):
                            if direction != 0:
                                print("format is broken: {0}\n".format(line))
                                exit -1
                            # found prefix-suffix-match
                            out.write("{0} + {1} + {2}\n".format(c2, c1, length))
                        else:
                            print("format is broken: {0}\n".format(line))
                            exit -1
                else:
                    print("format is broken: {0}\n".format(line))
                    exit -1
