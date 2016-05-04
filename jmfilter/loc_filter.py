#!/usr/bin/env python
import argparse
import logging
import sys
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s')

class LocFile(object):
    """
        LocFile 
    """

    def load_file(self):
        """
            Read loc file into python.
        """
        self._nloc = 0   
        self._header = []
        self._lines = []
        with open(self._filename) as loc_input:
            for i, line in enumerate(loc_input):
                if i < 4:
                    self._header.append(line.strip())
                else:
                    if (i % 2) == 0:
                        stripped = line.strip()
                        stripped += "\t"
                        self._nloc += 1
                    else:
                        stripped += line.strip()
                        self._lines.append((stripped))
        logging.info("Loaded linkage map")
        logging.info("Number of loci = %s" % str(self._nloc))

    def __init__(self, loc_file):
        self._filename = loc_file
        self.load_file() 

    def get_header(self):
        """
            Get header for the output locus file. 
        """
        self._header[2] = "nloc\t=\t" + str(self._nloc)
        return "\n".join(self._header) + "\n"

    def update_lines(self, lines):
        self._lines = lines
        nloc = 0
        logging.info("NLOC before update: %s" % self._nloc)
        for line in self._lines:
            nloc += 1
        self._nloc = nloc
        logging.info("NLOC after update: %s" % self._nloc)

    @property 
    def lines(self):
        return self._lines

__VALID__ = ["<nnxnp>","<lmxll>","<hkxhk>"]

def remove_missing_parent_markers(loc_file):
    """
        Identify lines with input that makes sense.
    """
    logging.info("Remove data with missing parent information")
    lines =  []
    total = 0
    for i, line in enumerate(loc_file.lines):
        split_line = line.split("\t")
        id = split_line[0]
        p_genotype = split_line[1]
        if any(p_genotype in s for s in __VALID__):
            lines.append(line) 
            total += 1
        #elif "SNP" not in id: 
        #    lines.append(line)
        #    total += 1
    loc_file.update_lines(lines)


def remove_bad_call_rate(loc_file, missing):
    """
        Identify and remove sites with a low call rate
    """

    logging.info("Remove sites that have a bad call rate")
    lines = []
    call_rate_thresh = missing
    for i, line in enumerate(loc_file.lines):
        split_line = line.split("\t")
        id = split_line[0]
        genotypes = split_line[2:]
        total_missing = sum([x == "--" for x in genotypes])
        if total_missing < (call_rate_thresh * len(genotypes)):
            lines.append(line)

    loc_file.update_lines(lines)

from scipy.stats import chisquare
import math
import numpy as np

def check_mendel(genotypes):
   
    groups = {}
    total_gt = 0
    for gt in genotypes:
        if gt != "--":
            try:
                groups[gt] += 1
            except KeyError:
                groups[gt] = 1
            total_gt += 1
    if (groups.keys()[0] in ["hk","hh","kk"]):
        keys = groups.keys()
        try: 
            idx = keys.index("hk")
        except:
            groups['hk'] = 0
            pass
        try: 
            idx = keys.index("kk")
        except:
            groups['kk'] = 0
            pass
        try: 
            idx = keys.index("hh")
        except:
            groups['hh'] = 0
            pass


    if len(groups.keys()) == 3:
        expected = np.array([0.25,0.5,0.25]) * total_gt
        observed = np.array([0,0,0])
        i= 0 
        for gt, counts in groups.items():
            partial_chi = counts 
            observed[i] = partial_chi 
            i += 1
         
    elif len(groups.keys()) == 2:
        expected = np.array([0.5,.5]) * total_gt
        observed = np.array([0,0]) 
        i= 0 
        for gt, counts in groups.items():
            partial_chi = counts 
            observed[i] = partial_chi 
            i += 1
    else:
        return 0.5
    p_value = chisquare(observed, expected)[1]
    if len(groups.keys()) == 0:
        return float("nan")

    #print id, p_value , groups, observed, expected 
    #if p_value <= 1:
    #    keys = groups.keys()
    #    oe = observed - expected
    #    if len(keys) == 3:
    #        print oe[1]
    #    elif len(keys) == 2:
    #        print oe[0]
    #        try:
    #            idx =keys.index("lm")
    #        except:
    #            pass
    #        try:
    #            idx =keys.index("np")
    #        except:
    #            pass
    #        try: 
    #            idx = keys.index("hk")
    #            print oe[idx]
    #        except:
    #            pass

                
    return(p_value)

def check_impossible_genotypes(loc_file,f2_pvalue, backcross_pvalue):

    logging.info("Look for impossible genotypes")
    lines = []
    for i, line in enumerate(loc_file.lines):
        split_line = line.split("\t")
        p_genotype = split_line[1] 
        genotypes = split_line[2:]
        id = split_line[0]
        if p_genotype == "<nnxnp>":
            p_value= check_mendel(genotypes) 
            if p_value > backcross_pvalue:
                lines.append(line)
        elif p_genotype == "<lmxll>":
            p_value=  check_mendel(genotypes)
            if p_value > backcross_pvalue:
                lines.append(line)
        else:
            p_value = check_mendel(genotypes)
            if p_value > f2_pvalue: 
                lines.append(line)
    loc_file.update_lines(lines)

def output_loc_map_loc_files(loc_file, output_file):
    """
        Output locus map file
    """
    output_loc = output_file 
    output_loc_f = open(output_loc, 'w')
    output_loc_f.write(loc_file.get_header())
    for line in loc_file.lines:
        output_loc_f.write(line + "\n")
        split_line = line.split('\t')

def main():
    """
        Main function for loc filter.
    """
    parser = argparse.ArgumentParser(description="Loc_file_filter")
    parser.add_argument("loc_file", help="Loc file to format")
    parser.add_argument("-m", "--missing", dest="missing", help="Maximum proportion of missing data allowed (default = 0.5)", default=0.5)
    parser.add_argument("-b", "--backcross-pvalue", dest="backcross_pvalue",
                        help="Segregant Chi-square P-value filter for maternally and paternally"
                        "informative back crosses: 1:1 expected ratio (default = 0.0001)", default=0.0001)
    parser.add_argument("-f", "--f2-pvalue", dest="f2_pvalue", help="Segregant Chi-square P-value"
                        "filter for`f2 crosses: 1:2:1 expected ratio (default = 0.001)", default=0.001)
    parser.add_argument("-o", "--output_file", dest="output_file", help="Output locus file" ,default="output_loc.loc")  
    logging.info("Started LOC filter for linkage map creation")
    args = parser.parse_args()
    try: 
        args.backcross_pvalue = float(args.backcross_pvalue)
        args.missing = float(args.missing)
        args.f2_pvalue = float(args.f2_pvalue)
    except:
        logging.error("Could not convert number arguments please check you command-line arguments")
        sys.exit(1)
    loc_file = args.loc_file
    loc_file = LocFile(loc_file)
    remove_missing_parent_markers(loc_file)
    remove_bad_call_rate(loc_file, args.missing)
    check_impossible_genotypes(loc_file, args.f2_pvalue,args.backcross_pvalue)
    output_loc_map_loc_files(loc_file, args.output_file)
if __name__ == "__main__":
    main()
