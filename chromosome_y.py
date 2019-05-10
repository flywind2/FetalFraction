from collections import defaultdict
import sys

# Argv 1 - from file
# Argv 2 - to file
# Argv 3 - CHR Y PERFECT REGIONS (BED FILE)
# Argv 4 - BIN SIZE

bin_size = int(sys.argv[4])
CHR_TOTAL_LENGTH = 2633553926  # SUM of all chromosomes/{13,18,21,X,Y}


def convert(good_region):
    good_region = good_region.split("\t")
    start, end = int(good_region[1]), int(good_region[2])
    return start, end, abs(end - start)


# Read good Y regions from BED file -> (start, end, length)
with open(sys.argv[3], "r") as f:
    Y_GOOD_REGION = list(map(convert, f.readlines()))

Y_GOOD_REGION_LENGTH = sum([length for (_, _, length) in Y_GOOD_REGION])  # - Y good region length # 57227415 - Y length

bins_per_chromosome = list(
    map(lambda x: (x[0], int(x[1] / bin_size)),
        [["chr 1", 248956422],  # chr_name <-> chr length ----> chr_name <-> bin count
         ["chr 2", 242193529],
         ["chr 3", 198295559],
         ["chr 4", 190214555],
         ["chr 5", 181538259],
         ["chr 6", 170805979],
         ["chr 7", 159345973],
         ["chr 8", 145138636],
         ["chr 9", 138394717],
         ["chr 10", 133797422],
         ["chr 11", 135086622],
         ["chr 12", 133275309],
         ["chr 13", 114364328],
         ["chr 14", 107043718],
         ["chr 15", 101991189],
         ["chr 16", 90338345],
         ["chr 17", 83257441],
         ["chr 18", 80373285],
         ["chr 19", 58617616],
         ["chr 20", 64444167],
         ["chr 21", 46709983],
         ["chr 22", 50818468],
         ["chr X", 156040895],
         ["chr Y", 57227415]]))


# Find normalized value
def is_trisomy_chr(chr_):
    return chr_ == "chr 13" or chr_ == "chr 18" or chr_ == "chr 21" or chr_ == "chr X" or chr_ == "chr Y"


def get_writable_line(patient_name, read_counts, Y_count):
    # Patient total read count of all chromosomes/{13,18,21,X,Y}
    read_count = sum([sum(read_counts[chr_].values()) for (chr_, _) in bins_per_chromosome if not is_trisomy_chr(chr_)])

    # Patient total read count normalized. Y chr is haploid, autosomals are diploids.
    normalized_chr_count = (read_count / 2) / CHR_TOTAL_LENGTH
    norm_Y = Y_count / Y_GOOD_REGION_LENGTH

    # Calculate fetal fraction
    ff = norm_Y / normalized_chr_count * 100 if norm_Y != 0 else 0

    # Generate line for file
    to_file = patient_name + ";" + str(Y_count) + ";" + str(normalized_chr_count) + ";" + str(norm_Y) + ";" + str(
        ff) + ";" + str(read_count) + ";" + str(Y_GOOD_REGION_LENGTH)

    for chr_, bins in bins_per_chromosome:
        if not is_trisomy_chr(chr_):
            for i in range(bins):
                to_file += ";" + ("0" if read_counts[chr_][i] is None else str(read_counts[chr_][i]))
    return to_file


with open(sys.argv[1]) as f, open(sys.argv[2], "a", encoding="UTF-8") as out:
    old_patient = ""
    read_count_dict = defaultdict(lambda: defaultdict(int))
    Y_count = 0

    for line_number, line in enumerate(f):
        try:
            line = line.strip().split(";")
            patient, read_chromosome, read_pos = line[0], line[1], int(line[2])
        except IndexError:
            print(line)
            continue

        if patient == old_patient:
            # Add 1 to correct bin in correct chromosome in read counts
            read_count_dict[read_chromosome][(int(read_pos / bin_size))] += 1

            # Increase Y_count if read in suitable position
            for start, end, length in Y_GOOD_REGION:
                if start <= read_pos <= end and read_chromosome == "chr Y":
                    Y_count += 1
                    break  # Actually not needed, should not have overlapping regions

        else:
            if line_number != 0:
                out.write(get_writable_line(old_patient, read_count_dict, Y_count) + "\n")

            old_patient = patient
            read_count_dict = defaultdict(lambda: defaultdict(int))
            Y_count = 0

    out.write(get_writable_line(old_patient, read_count_dict, Y_count) + "\n")


def get_writable_header():
    to_file = "patient;Y total count;norm all ch;norm good Y;FF;all;Y good region length"
    for chromosome, bins in bins_per_chromosome:
        if not is_trisomy_chr(chromosome):
            for i in range(bins):
                to_file += ";" + (chromosome.replace("chr ", "c") + "b" + str(i))
    return to_file


# Write header to separate file
with open("header.txt", "w", encoding="UTF-8") as header:
    header.write(get_writable_header() + "\n")
