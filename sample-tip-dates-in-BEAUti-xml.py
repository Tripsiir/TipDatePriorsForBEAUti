#! /usr/bin/env python

'''
Author: Pieter Moris 27/06/2016

This script takes a.BEAUti .xml file as input and adds
prior distribution id's, sample operators and logger entries 
for each of the samples/taxa in the tree.
The input xml file should be generated in the following manner:
    - Import the alignment
    - Use tip dates -> since before the present -> guess from file lsd_input_dates.txt
        file can be created using the lsdDataInputFromSeqs.py
    - Set site/clock model, priors and MCMC options. 
    - Make sure that any other taxon sets do not follow the naming convention 
      that is used by this script: "tip.taxonid".

The output is saved in the working directory as update.BEAUti.xml.

describe input seq list format

POSSIBLE IMPROVEMENTS:
    - Retrieve the sample id's from the lsd-input-dates.txt file instead of the xml file directly.
      This would prevent entries being created for sequences that don't have a known date.
        + check if all provides id's are actually in the xml file.
        OR extract from    <trait id="dateTrait.t:renamed-SNP-L4" spec="beast.evolution.tree.TraitSet" traitname="date"> section?
    - Instead of opening the .xml file three times, extract all the required information in one go.
    - Error catching: what if .xml file is not complete? E.g. missing the prior or logger sections?
    - Add option to disable logger section.
    - check if dates are provided as time in past
     <trait id="dateTrait.t:renamed-SNP-L4" spec="beast.evolution.tree.TraitSet" traitname="date-forward">
                 <trait id="dateTrait.t:renamed-SNP-L4" spec="beast.evolution.tree.TraitSet" traitname="date-backward"> this is what we need
    - add optional prior distribution option

Help: python.BEAUti-xml-sample-from-tip-dates.py -h

Usage: python.BEAUti-xml-sample-from-tip-dates.py.BEAUti.xml
'''

from __future__ import print_function, division, with_statement
import sys, os, re, time, argparse

# Check provided arguments
parser = argparse.ArgumentParser(description="Script to add tip date prior "
                                             "sampling to beati .xml files.")
parser.add_argument("inputFile", type=str, help="The path or name "
    # https://stackoverflow.com/questions/18862836/how-to-open-file-using-argparse
    # type=argparse.FileType('r'),
                    "(if located in working directory) of the input.BEAUti .xml file.",
                    metavar="inputFile")
parser.add_argument("-o", "--output", type=argparse.FileType('w'), dest="outFile", 
                    default='updated-beati.xml',
                    help="Specifies the name or path of the updated.BEAUti .xml file.",
                    metavar='')
parser.add_argument("-s", "--sequences", type=str, dest="sequences", 
                    default='',
                    help="A file that specifies the sequences to which "
                    "date priors should be assigned. Each line should "
                    "contain exactly one sequence id (default = all dated "
                    "sequences found in the .BEAUti .xml file).",
                    # metavar="Sequences to assign date priors to")
                    metavar='')
parser.add_argument("-d", "--priorDist", type=str.lower, dest="priorDistribution",
                    default="exponential",
                    choices=["exponential", "poisson", "log-normal", "beta",
                             "gamma", "inverse-gamma", "laplace", "uniform", 
                             "normal", "1/x", ],
                    help="The prior distribution to be used (default = exponential). "
                    "Possibilities are: 'exponential', 'poisson', 'log-normal', "
                    "'gamma', 'inverse-gamma', beta', 'laplace', 'uniform', "
                    "'normal' and '1/x'."
                    "\n Please refer to BEAUti for more details on these distributions.",
                    # metavar="Prior distribution")
                    metavar='')
parser.add_argument("-p1", "--parameter1", type=float, dest="parameter1",
                    default=1.0,
                    help="The first parameter (mean, alpha, lambda, mu "
                    "lower or M) the prior distribution (default = 1).",
                    # metavar="First parameter of prior distribution")
                    metavar='')
parser.add_argument("-p2", "--parameter2", type=float, dest="parameter2",
                    default=2.0,
                    help="The second parameter (sigma, S, beta, scale or "
                    "upper) of the prior distribution (default = 2).",
                    # metavar="Second parameter of prior distribution")
                    metavar='')
parser.add_argument("-po", "--parametero", type=float, dest="parametero",
                    default=0.0,
                    help="The off-set parameter of the prior distribution (default = 0).",
                    # metavar="Offset parameter of prior distribution")
                    metavar='')
parser.add_argument("-r", "--realspace", dest="meanInRealSpace",
                    action="store_true",
                    help="If flag is provided, the mean of the log normal "
                    "distribution is treated as being in real space, "
                    "rather than log-transformed space")
parser.add_argument("-e", "--estimate", dest="estimateParameters",
                    action="store_true",
                    help="If flag is provided, all the parameters of the "
                    "chosen distribution will be estimated.")
args = parser.parse_args()

# Map input distributions to beast distribution names
beastDistributionsDict = {"exponential": "Exponential",
                          "log-normal": "LogNormal",
                          "gamma": "Gamma",
                          "beta": "Beta",
                          "inverse-gamma": "InverseGamma",
                          "laplace": "LaplaceDistribution",
                          "normal": "Normal",
                          "1/x": "OneOnX"}
                          # poisson and uniform have different section:
                          # grep 'distribution id="tip.AUSA2s61.prior"' * -A 7

# Define parameters for different distributions
parameterNamesDict = {"exponential": ["mean", "offset"],
                      "log-normal": ["M", "S", "offset"],
                      "gamma": ["alpha", "beta", "offset"],
                      "beta": ["alpha", "beta", "offset"],
                      "inverse-gamma": ["alpha", "beta", "offset"],
                      "poisson": ["lambda", "offset"],
                      "laplace": ["mu", "scale", "offset"],
                      "1/x": ["offset"],
                      "normal": ["mean", "sigma", "offset"],
                      "uniform": ["lower", "upper", "offset"]}

# Store parameters in list
parameterList = [args.parameter1, args.parameter2, args.parametero]

# Extract path from input argument
path = os.path.abspath(args.inputFile)

# Check file extension
if not os.path.splitext(path)[1].lower() == '.xml':
    raise Exception('Not a .xml file.')
'''
if not os.path.splitext(path)[1].lower() == '.xml':
    print("ERROR!\nProvided file does not have the correct .xml extension.")
    sys.exit()
'''

# Read file into memory as a continuous string
try:
    with open(path, "r") as inputXML:
        xmlContents = inputXML.read()
except EnvironmentError:
    print("Input file not found. Please check the path/filename and try again.")
    sys.exit()

# Extract tree ID from file
try: 
    treeID = re.search(r'<tree id="(?P<treeID>[\w:.-]+)"', xmlContents).group("treeID")
except AttributeError:
    print("ERROR!\nFailed to find tree id in BEAUti .xml file.\n"
          "Please check if .xml file was generated correctly.")
    sys.exit()

# Find section containing dated sequences
try:
    dateSection = re.search(r'traitname="date-backward">\s*([\w=.\-\s,]*)<',xmlContents).group(1)
except AttributeError:
    print("ERROR!\nCould not find any dated sequences in BEAUti .xml file.\n"
          "Please check if.BEAUti .xml file was generated correctly.\n"
          "Possible causes:\n - Tip dates were not added in BEAUti.\n"
          " - Dates were not specified as \"before the present\" ")
    sys.exit()

# Extract sequence id's from this section
taxonList = re.findall(r'([\w\-]*)=', dateSection)

# Only keep specified sequences if a sequence input file was provided
if args.sequences: 
    try: 
        with open(os.path.abspath(args.sequences)) as sequenceFile:
            sequenceList = sequenceFile.read().splitlines()

            if set(sequenceList).issubset(taxonList):
                taxonList = sequenceList
            else:
                print("Provided sequence file contained sequence id's which "
                      "are not present in the BEAUti .xml file. Please try again.")
                sys.exit()

    except EnvironmentError:
        print("Failed to open input sequence file. "
              "Please check the path/filename and try again.")
        sys.exit()

print("\nPPoTD - Putting Priors on Tip Dates")
print("Preparing to add required sections to the beati .xml file.\n")
print("The following options will be used:")
print("- Prior distribution:", args.priorDistribution)
for i, param in enumerate(parameterNamesDict[args.priorDistribution]):
    print("\t -", param, "=", parameterList[i])
if args.priorDistribution == 'log-normal':
    if args.meanInRealSpace:
        print("\t - Mean is specified in real space.")
    else:
        print("\t - Mean is specified in log-transformed space.")
print("\nInput BEAti .xml file:",path)
print("Output BEAti .xml file:",os.path.normcase(os.getcwd() + '/' + args.outFile.name))
print("\nThe input BEAUti .xml file seems to be in order...\n")
if args.sequences:
    print("Adding prior distributions, sample operators and logger entries for the tip dates",
          str(len(taxonList)), "sequences provided in the input file.")
else:
    print("Adding prior distributions, sample operators and logger entries for the tip dates of all",
          str(len(taxonList)), "dated sequences found in the BEAUti .xml file.")