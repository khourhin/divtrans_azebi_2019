import logging
import tempfile
import subprocess
import os
import pandas as pd
import shutil


def simplify_outpath(filename, prefix="", suffix="", keep_path=False, check_exist=True):
    """
    Produce a name for an output file
    """

    output_name = os.path.basename(filename)
    path = os.path.dirname(filename)

    # Remove all extensions
    output_name = os.path.splitext(output_name)[0]

    output_name = prefix + output_name + suffix

    if keep_path:
        output_name = os.path.join(path, output_name)

    if check_exist and os.path.exists(output_name):
        logging.error("File already exists: {}".format(output_name))
        raise IOError("File already exists: {}".format(output_name))

    return output_name


def exec_command(cmd, silent=False):
    """ Wrapper for proper execution of subprocesses
    """

    try:
        proc = subprocess.run(
            cmd,
            shell=True,
            check=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True,
        )
        if not silent and proc.stdout:
            logging.debug(proc.stdout.strip())

    except subprocess.CalledProcessError as exc:
        raise ChildProcessError(exc.stderr)


def df_to_bed(df, path):
    """ Create a bed object from a pandas dataframe
    """

    df.to_csv(path, sep="\t", header=None, index=False)
    return Bed(path)


class Bed(object):
    """Bed object
    """

    def __init__(self, path):

        if not os.path.exists(path):
            raise IOError("The bed file {} does not exist".format(path))

        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]

    def __len__(self):
        return sum(1 for _ in open(self.path))

    def __repr__(self):
        return f"<Bed object: {self.name}>\n"

    def head(self, nlines=10):
        return self.to_df().head(nlines)

    def move_to(self, newpath):
        """Move the bed file associated
        """
        shutil.move(self.path, newpath)
        self.path = newpath
        self.name = os.path.basename(self.path)

        return self

    def to_df(self):
        """ Bed file to pandas dataframe
        """
        return pd.read_csv(self.path, sep="\t", header=None, low_memory=False)

    def sort(self, outfolder="bed_outfolder"):
        """ Sort bed file
        """
        os.makedirs(outfolder, exist_ok=True)
        sort_path = os.path.join(outfolder, self.name + "_S.bed")

        cmd = "sort -k1,1 -k2,2n {0} > {1}".format(self.path, sort_path)

        exec_command(cmd)

        return Bed(sort_path)

    def merge(self, outfolder="bed_outfolder", supp_args=""):
        """ Merge bed
        """
        os.makedirs(outfolder, exist_ok=True)
        merged_path = os.path.join(outfolder, self.name + "_M.bed")

        cmd = "bedtools merge -i {0} {1} > {2}".format(
            self.path, supp_args, merged_path
        )
        exec_command(cmd)

        return Bed(merged_path)

    def concat(self, bedobj, outfolder="bed_outfolder"):
        """Concatenate 2 bed files
        """

        os.makedirs(outfolder, exist_ok=True)
        concat_path = os.path.join(
            outfolder, self.name + "_concat_" + bedobj.name + ".bed"
        )

        cmd = "cat {0} {1} | sort -k1,1 -k2,2n > {2}".format(
            self.path, bedobj.path, concat_path
        )
        subprocess.check_output(cmd, shell=True)

        exec_command(cmd)

        return Bed(concat_path)

    def closest(self, bed_obj, outfolder="bed_outfolder", supp_args=""):
        """ Make default bedtools closest
        """

        os.makedirs(outfolder, exist_ok=True)
        closest_path = os.path.join(
            outfolder, self.name + "-closest-" + bed_obj.name + ".bed"
        )

        cmd = "bedtools closest -a {0} -b {1} {2} > {3}".format(
            self.path, bed_obj.path, supp_args, closest_path
        )
        exec_command(cmd)

        return Bed(closest_path)


class DivTransFromBed(object):
    """Object to store methods for divergent transcription detection
    readng with paired reads with bedtools fragments """

    def __init__(self, bam, name_sorted=False):
        """
        :param bam: A name sorted bam file
        :param name_sorted: If the bam files are already name sorted.
        """
        self.bam = bam
        self.name_sorted = name_sorted
        self.tmp_dir = tempfile.mkdtemp(prefix="divtrans_")
        logging.info(f"Intermediate files will be stored in {self.tmp_dir}")

    def __repr__(self):
        return f"DivTrans by beds Object from: {simplify_outpath(self.bam)}"

    def sort_by_name(self):
        """ Sort the bam file by read name 
        """

        logging.debug("START: Sorting by names")
        sorted_bam = simplify_outpath(self.bam, suffix="_nameSort.bam", keep_path=True)
        cmd = f"samtools sort -n -o {sorted_bam} {self.bam}"

        exec_command(cmd)

        self.bam = sorted_bam

    def filter(self):
        """ Filters bam file based on mapping, mapping of mate, splicing
        """

        logging.debug("START: Filtering bams")
        self.filt_bam = os.path.splitext(self.bam)[0] + "_filt.bam"

        # Filtering bam, removing unmapped, mate unmapped, spliced, secondary mapping
        cmd = "samtools view -h -F 0x4 -F 0x8 -F 0x100 {0} | awk '{{if ($1 ~ /^@/) {{print}} else if ($6 !~ /N/) {{print}}}}' | samtools view -bh > {1}".format(
            self.bam, self.filt_bam
        )

        subprocess.check_output(cmd, shell=True)
        logging.debug("DONE: {}".format(cmd))

        return self

    def to_bedpe(self):
        """From mapped reads in a bam, get bed intervals of the pairs joined
        into fragments
        """

        logging.debug("START: Generate bedpe")
        self.bedpe = os.path.join(
            self.tmp_dir, (simplify_outpath(self.bam, suffix="_frag.bed"))
        )
        # Converting to fragments and bed format
        cmd = "bedtools bamtobed -bedpe -mate1 -i {0} > {1} 2> bedpe.log" "".format(
            self.filt_bam, self.bedpe
        )

        subprocess.check_output(cmd, shell=True)
        logging.debug("DONE: {}".format(cmd))

        return self

    def to_fragments(self, max_insert_size=500):

        logging.debug("START: Generate fragments bed")

        diff_chroms_count = 0
        frag_lengths = []

        bed_out = os.path.splitext(self.bedpe)[0] + "_clean.bed"

        frag_filt_n = 0

        with open(bed_out, "w") as f:
            with open(self.bedpe, "r") as bedin:
                for l in bedin:
                    # Check if different chromosomes
                    if l.split()[0] != l.split()[3]:
                        diff_chroms_count += 1
                        continue

                    coords = [int(l.split()[i]) for i in [1, 2, 4, 5]]
                    start = min(coords)
                    stop = max(coords)

                    frag_l = stop - start
                    frag_lengths.append(frag_l)

                    # This is setting the correct strand (strand of the second of pair)
                    if frag_l <= max_insert_size:
                        f.write(
                            "{0}\t{1}\t{2}\t{3}\t.\t{4}\n".format(
                                l.split()[0], start, stop, l.split()[6], l.split()[9]
                            )
                        )
                    else:
                        frag_filt_n += 1

        logging.info(
            "Number of pairs mapping on different chromosomes: {}".format(
                diff_chroms_count
            )
        )
        logging.info(
            "Min frag_length: {0}; Max frag_length: {1};".format(
                min(frag_lengths), max(frag_lengths)
            )
        )
        logging.info(
            "Number of fragments filtered out because insert size > {0}b: {1}".format(
                max_insert_size, frag_filt_n
            )
        )

        self.fragbed = Bed(bed_out)
        return self

    def merge(self):
        """ Merge bed fragments keeping the strand info
        """

        logging.debug("START: Merging bed")

        # -c 6 -o distinct is to keep the strand info
        self.mergedbed = self.fragbed.sort(outfolder=self.tmp_dir).merge(
            supp_args="-s -c 4,5,6 -o min,distinct,distinct", outfolder=self.tmp_dir
        )

        return self

    def get_div_reads(self):
        """ Get couples of closest reads on opposite strands"""

        df = self.mergedbed.to_df()

        # Separate + and - strands
        pos_strand_df = df[df[5] == "+"]
        neg_strand_df = df[df[5] == "-"]

        pos_bed = df_to_bed(
            pos_strand_df,
            simplify_outpath(self.mergedbed.path, suffix="_plus.bed", keep_path=True),
        )
        neg_bed = df_to_bed(
            neg_strand_df,
            simplify_outpath(self.mergedbed.path, suffix="_minus.bed", keep_path=True),
        )

        # The divergent read should always be looked upstream considering
        # the orientation of bedA (see bedtools closest manual)
        bed1 = pos_bed.closest(neg_bed, supp_args="-D a -id", outfolder=self.tmp_dir)
        bed2 = neg_bed.closest(pos_bed, supp_args="-D a -id", outfolder=self.tmp_dir)
        closest_bed = bed1.concat(bed2, outfolder=self.tmp_dir).move_to(
            simplify_outpath(
                self.mergedbed.path, suffix="_closest_concat.bed", keep_path=True
            )
        )

        # Ignore overlaping reads to look for reads further away
        # Treat special cases where overlapping hide divergent transcription (test with test/frag.bed)
        # This will overwrite the previous beds so need to be concat sequentially
        bed3 = pos_bed.closest(
            neg_bed, supp_args="-D a -id -io", outfolder=self.tmp_dir
        )
        closest_bed = closest_bed.concat(bed3, outfolder=self.tmp_dir).move_to(
            simplify_outpath(
                self.mergedbed.path,
                suffix="_closest_concat.bed",
                keep_path=True,
                check_exist=False,
            )
        )
        bed4 = neg_bed.closest(
            pos_bed, supp_args="-D a -id -io", outfolder=self.tmp_dir
        )
        closest_bed = closest_bed.concat(bed4, outfolder=self.tmp_dir).move_to(
            simplify_outpath(
                self.mergedbed.path,
                suffix="_closest_concat.bed",
                keep_path=True,
                check_exist=False,
            )
        )

        self.closest_bed = closest_bed

        return self

    def get_div_read_intervals(self, row, internal=True):
        """ Defines the bed intervals of divergent transcription event
        Either defined as the interval separating divergetn reads (internal=True)
        Or as the min and max of divergent reads coordinates (internal=False)"""

        chro = row.iloc[0]
        coords = sorted(row.iloc[[1, 2, 7, 8]])
        if internal:
            return [chro, coords[1], coords[2], ".", ".", "."]
        else:
            return [chro, min(coords), max(coords), ".", ".", "."]

    def cleanup_intervals(self):
        """ Process fragments beds to extract divergent transcription intervals"""

        logging.debug("START: Divergent intervals detection and cleaning")

        def filter_convergence(row):
            """ Keep only divergent transcription (no convergent)
            """
            # 2 cases for divergent transcription (hence removing convergente transcription)
            if row.iloc[1] <= row.iloc[7] and row.iloc[5] == "-":
                return True
            elif row.iloc[1] >= row.iloc[7] and row.iloc[5] == "+":
                return True
            else:
                return False

        df = self.closest_bed.to_df()

        # Filtering convergence (still detected when reads are overlapping)
        df = df.loc[df.apply(filter_convergence, axis=1)]

        # Filter by event size
        df = df[(abs(df[12]) <= 500)]

        # remove GL contigs
        df = df[[not str(x).startswith("GL") for x in df.iloc[:, 0]]]

        # For cases where no close reads can be found, bedtools add -1 as coordinates
        # This filters this case out
        sel = (
            (df.iloc[:, 1] != -1)
            & (df.iloc[:, 2] != -1)
            & (df.iloc[:, 7] != -1)
            & (df.iloc[:, 8] != -1)
        )
        df = df.loc[sel, :]

        df = pd.DataFrame(list(df.apply(self.get_div_read_intervals, axis=1)))

        div_read_intervals = df_to_bed(
            df,
            simplify_outpath(
                self.mergedbed.path, suffix="_div_read.bed", keep_path=True
            ),
        )

        self.divtrans_bed = (
            div_read_intervals.sort(outfolder=self.tmp_dir)
            .merge(outfolder=self.tmp_dir)
            .move_to(
                simplify_outpath(self.bam, suffix="_bed_divtrans.bed", keep_path=True)
            )
        )

        return self

    def clean(self):
        """
        Remove intermediate files
        """
        shutil.rmtree(self.tmp_dir)

    def run(self):
        """ Filter the bam and produce a bed of fragments"""
        if not self.name_sorted:
            self.sort_by_name()
        # Remove read unmapped, mate unmapped, only primary, no spliced
        self.filter()
        # Generate fragments from read pairs and orient them properly
        self.to_bedpe()
        self.to_fragments()
        self.merge()
        self.get_div_reads()
        self.cleanup_intervals()

        return self
