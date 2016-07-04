#! /usr/bin/env python

'''
Author: Pieter Moris 27/06/2016

This script takes a .BEAUti .xml file as input and adds prior distribution id's, 
sample operators and logger entries for each of the samples/taxa in the tree.

The input xml file should be generated in the following manner:
    - Import the alignment in .BEAUti
    - Load in tip dates -> set time since before the present or time since the past
    - Set site/clock model, priors and MCMC options. 
    - Do not manually create any priors for the tip dates! 
    - Make sure that any other taxon sets do not follow the naming convention 
      that is used by this script: "tip.taxonid".

By default the output is saved to "updated-beauti.xml" in the current working directory.

A list of sequence id's for which to add tip priors can be provided as an input argument.
This should be a text file with each line containing the sequence id as specified in the alignment.

Help: python.BEAUti-xml-sample-from-tip-dates.py -h

Usage: python.BEAUti-xml-sample-from-tip-dates.py. BEAUti.xml

POSSIBLE IMPROVEMENTS:
    - Add option to disable logger section.
    - check if dates are provided as time in past
     <trait id="dateTrait.t:renamed-SNP-L4" spec="beast.evolution.tree.TraitSet" traitname="date-forward">
                 <trait id="dateTrait.t:renamed-SNP-L4" spec="beast.evolution.tree.TraitSet" traitname="date-backward"> this is what we need
    - add option for estimating prior distribution parameters
    - change beastDistributionsDict[args.priorDistribution] to stored variable name
    - remove argparse writer, manually open file and change default name path?
    - fix output names!
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
                    "rather than log-transformed space. ")
parser.add_argument("-t", "--time", type=str, dest="timeDirection",
                    default="present",
                    choices=["present","past"],
                    help="Specifies whether dates are interpreted as 'time before the present' "
                    "or 'time since the past'.",
                    metavar='')
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
                          "1/x": "OneOnX",
                          "uniform": "Uniform"}
                          # poisson has different section:
                          # grep 'distribution id="tip.AUSA2s61.prior"' * -A 7

# Define parameters for different distributions
parameterDict = {"exponential": {"mean": args.parameter1},
                      "log-normal": {"M": args.parameter1,
                                     "S": args.parameter2},
                      "gamma": {"alpha": args.parameter1,
                                "beta": args.parameter2},
                      "beta": {"alpha": args.parameter1,
                               "beta": args.parameter2},
                      "inverse-gamma": {"alpha": args.parameter1,
                                        "beta": args.parameter2},
                      "poisson": {"lambda": args.parameter1},
                      "laplace": {"mu": args.parameter1,
                                  "scale": args.parameter2},
                      "1/x": {"offset": args.parametero},
                      "normal": {"mean": args.parameter1,
                                 "sigma": args.parameter2},
                      "uniform": {"lower": args.parameter1,
                                  "upper": args.parameter2}}

# parameterNamesDict = {"exponential": ["mean", "offset"],
#                       "log-normal": ["M", "S", "offset"],
#                       "gamma": ["alpha", "beta", "offset"],
#                       "beta": ["alpha", "beta", "offset"],
#                       "inverse-gamma": ["alpha", "beta", "offset"],
#                       "poisson": ["lambda", "offset"],
#                       "laplace": ["mu", "scale", "offset"],
#                       "1/x": ["offset"],
#                       "normal": ["mean", "sigma", "offset"],
#                       "uniform": ["lower", "upper", "offset"]}
# Store parameters in list
# parameterList = [args.parameter1, args.parameter2, args.parametero]

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
    if args.timeDirection == 'present':
        dateSection = re.search(r'traitname="date-backward">\s*([\w=.\-\s,]*)<',xmlContents).group(1)
    else:
        dateSection = re.search(r'traitname="date-forward">\s*([\w=.\-\s,]*)<',xmlContents).group(1)
except AttributeError:
    print("ERROR!\nCould not find any dated sequences in BEAUti .xml file.\n"
          "Please check if.BEAUti .xml file was generated correctly.\n"
          "Possible causes:\n - Tip dates were not added in BEAUti.\n"
          " - Dates were specified as \"before the present\" or \"since some time in the past\" "
          "in the .XML file, but a different direction argument was provided as an input argument to the script.")
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
# for i, param in enumerate(parameterNamesDict[args.priorDistribution]):
#     print("\t -", param, "=", parameterList[i])
for parameter in parameterDict[args.priorDistribution]:
    print("\t -", parameter, "=", parameterDict[args.priorDistribution][parameter])
print("\t - offset =", args.parametero)

if args.priorDistribution == 'log-normal':
    if args.meanInRealSpace:
        print("\t - Mean is specified in real space.")
    else:
        print("\t - Mean is specified in log-transformed space.")
if args.timeDirection == "present":
    print("\nTip dates are interpreted as time before the present.")
else:
    print("\nTip dates are interpreted as since some time in the past.")
print("\nInput BEAti .xml file:",path)
print("\nOutput BEAti .xml file:",os.path.normcase(os.getcwd() + '/' + args.outFile.name))
print("\nThe input BEAUti .xml file seems to be in order...\n")
if args.sequences:
    print("Adding prior distributions, sample operators and logger entries for the tip dates",
          str(len(taxonList)), "sequences provided in the input file.")
else:
    print("Adding prior distributions, sample operators and logger entries for the tip dates of all",
          str(len(taxonList)), "dated sequences found in the BEAUti .xml file.")

# Split input file line by line into a list to allow for easy output writing
xmlContentsList = xmlContents.split('\n')

for i, line in enumerate(xmlContentsList):
    # Check if current section needs to be expanded with an operator,
    # logger or prior distribution for the sample tips.

    # Print progress bar to screen
    # time.sleep(0.001)
    # sys.stdout.write('\rProgress: {0:.2%}'.format(i / len(xmlContentsList)))
    # sys.stdout.flush()

    newLine = ''

    # Add prior distributions for each chosen tip date in the prior section of the BEAUti .xml file
    if '<distribution id="prior" spec="util.CompoundDistribution">' in line:
        for taxon in taxonList:

            # Create parameter id's for specified distribution
            ### NOTE: ignoring lower and upper in parameter id for S for log-normal!

            parameterLines = ''
            for parameter in parameterDict[args.priorDistribution]:
                parameterLines += '\t\t    <parameter id="RealParameter.' + parameter + '.' + taxon + \
                                  '" estimate="false" name="' + parameter + '">' + \
                                  str(parameterDict[args.priorDistribution][parameter]) + \
                                  '</parameter>\n'

            # Create distribution section to surround parameters
            if args.priorDistribution != "poisson": 
                distributionSection = '\t\t<' + beastDistributionsDict[args.priorDistribution] + \
                                      ' id="' + beastDistributionsDict[args.priorDistribution] + \
                                      '.' + taxon + '" name="distr" offset="' + str(args.parametero) + \
                                      '">\n' + parameterLines + '\t\t</' + \
                                      beastDistributionsDict[args.priorDistribution] + '>\n'
            # Requires slightly different syntax for Poisson distribution!
            else: 
                distributionSection = '\t\t<distr id="Poisson.' + taxon + '" spec="beast.math.distributions.Poisson" ' \
                                      'offset="1.0">\n' + parameterLines + '\t\t</distr>\n'

            # Create remainder of the prior distribution id body and insert the distributionSection
            newLine += '\t    <distribution id="tip.' + taxon + '.prior"' \
                       'spec="beast.math.distributions.MRCAPrior" ' \
                       'tipsonly="true" tree="@' + treeID + '">\n' \
                       '\t\t<taxonset id=" tip.' + taxon + '" spec="TaxonSet">\n' \
                       '\t\t    <taxon id="' + taxon + '" spec="Taxon"/>\n\t\t</taxonset>\n' + \
                       distributionSection + '\t    </distribution>\n'

        # Append newly created sections to the current line, i.e. insert them at
        # the top of the <distribution id="prior" spec="util.CompoundDistribution">' section
        line = line + "\n" + newLine

    # Add loggers for tip priors at the top of the logger section
    elif '<logger id="tracelog"' in line:
        for taxon in taxonList:
            newLine += '\t<log idref="@tip.' + taxon + '.prior"/>\n'
        line = line + "\n" + newLine

    # Add sample operator for tip priors at the end of the file
    elif '</run>' in line:
        newLine = ''
        for taxon in taxonList:
            newLine += '<operator id="TipDatesRandomWalker.' + taxon + \
                       '"\nwindowSize="1"\n' \
                       'spec="TipDatesRandomWalker"\n' \
                       'taxonset="@tip.' + taxon + \
                       '"\ntree="@' + treeID + \
                       '"\nweight="1.0"/>\n'
        line = newLine + line

    # write edited or original line to the new file
    args.outFile.write(line+"\n")

'''
add option for estimating dist parameters!

fix formatting with tabs

add loggers  before </logger>? no because screenlog has same name!
'''

# Set progess bar to done.
# sys.stdout.write('\r{0:.2%} done.'.format(1))

# Print name and location of new.BEAUti .xml file.
print('\nUpdated .BEAUti .xml file was saved to',
     os.path.normcase(os.getcwd() + "/updated-" + args.outFile.name))

#print('\nUpdated .BEAUti .xml file was saved to',
#      os.path.normcase(os.getcwd() + "/SampledTips-" + os.path.basename(path)))

sys.exit()