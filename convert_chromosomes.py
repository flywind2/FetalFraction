import sys

# Argv 1 - from file
# Argv 2 - to file

file = sys.argv[1]
to_write = sys.argv[2]

known = {"CM000663.2": "chr 1", "CM000664.2": "chr 2", "CM000665.2": "chr 3", "CM000666.2": "chr 4",
         "CM000667.2": "chr 5", "CM000668.2": "chr 6", "CM000669.2": "chr 7", "CM000670.2": "chr 8",
         "CM000671.2": "chr 9", "CM000672.2": "chr 10", "CM000673.2": "chr 11", "CM000674.2": "chr 12",
         "CM000675.2": "chr 13", "CM000676.2": "chr 14", "CM000677.2": "chr 15", "CM000678.2": "chr 16",
         "CM000679.2": "chr 17", "CM000680.2": "chr 18", "CM000681.2": "chr 19", "CM000682.2": "chr 20",
         "CM000683.2": "chr 21", "CM000684.2": "chr 22", "CM000685.2": "chr X", "CM000686.2": "chr Y"}

with open(file) as f, open(to_write, "w") as out:
    for line in f:
        try:
            splitted = line.split("\t")
            chromosome, pos = splitted[2], splitted[3]
            out.write(file.split("/")[-2].strip() + ";" + known[chromosome] + ";" + pos + "\n")
        except:
            pass
