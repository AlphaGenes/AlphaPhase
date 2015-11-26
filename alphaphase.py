#!/usr/bin/python2.7
# encoding: utf-8
'''
alphaphase -- Generate the Spec file for AlphaPhase
Different
@author:     Roberto Antolín
@copyright:  2015 Roberto Antolín. All rights reserved.
@license:    license
@contact:    roberto dot antolin at roslin dot ed dot ac dot uk
@deffield    updated: Updated
'''

import sys, os
import subprocess

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.2
__date__ = '2015-06-03'
__updated__ = '2015-06-03'

DEBUG = 0
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def main(argv=None):    # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
  Created by rantolin on %s.
  Copyright 2015 Roslin Institute. All rights reserved.
  Licensed under the GPL 3.0v
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument(dest="output", help="Output file", metavar="file", type=str)
        parser.add_argument("-P", "--pedigree", dest="pedigree",
            help="File containing the pedigree information", metavar="file", type=file, required=False)
        parser.add_argument("-G", "--genotype", dest="genotype",
            help="File containing the genotypes", metavar="file", type=file)
        parser.add_argument("-f", "--format", help="Format of the genotypes [Default: %(default)s]", metavar=str,
            required=False, dest="genoFormat", choices=["GenotypeFormat", "UnorderedFormat"], default="GenotypeFormat")
        parser.add_argument("-S", "--snp", dest='nSNPs',
            help="Number of SNP in the genotype file", type=int, metavar="nSNP", required=True)
        parser.add_argument("-t", "--tails", help="Core and Tail lengths [Default: %(default)d]",
            default=100, type=int, metavar="int", dest='tails')
        parser.add_argument("-c", "--cores", dest='cores',
            help="Core lengths [Default: %(default)d]", type=int, default=100, metavar="int")
        parser.add_argument("--offset", help='Create offset between cores',
            action="store_true", dest='offset')
        parser.add_argument("-n", "--surrogates", help='Number of surrogates [Default: %(default)d]',
            type=int, default=10, metavar=int, dest='surrogates')
        parser.add_argument("-d", "--PercentageSurrDisagree",
            help="Percentage of surrogates in disagreement [Default: %(default)3.1f]",
            type=float, default=10.0, dest="PercentageSurrDisagree", metavar="float") 
        parser.add_argument("-g", "--PercentageGenoHaploDisagree",
            help="Percentage of genos in disagreement across SNPs [Default: %(default)3.1f]",
            type=float, default=00.0, dest="PercentageGenoHaploDisagree", metavar="float")
        parser.add_argument("-e", "--GenotypeMissingErrorPercentage",
            help="Percentage of SNPs missingdisagreement across surrogates [Default: %(default)3.1f]",
            type=float, default=00.0, dest="GenotypeMissingErrorPercentage", metavar="float")
        parser.add_argument("-r", "--NrmThresh",
            help="Max. coeff. of relationship for surrogates parents [Default: %(default)3.1f]",
            type=float, default=0.0, dest="NrmThresh", metavar="float")
        parser.add_argument("--fullOutput", help='Extra output files are suppressed',
            action="store_true", dest='fullOutput')
        parser.add_argument("--graphics", help='Graphical output',
            action="store_true", dest='graphics')
        parser.add_argument("--simulation",
            help='Analysis involves simulated data where the true phase is known',
            action="store_true", dest='simulation')
        # parser.add_argument("-M", "--hmm",
        #     help="Use Hidden Markov Model [Default: %(default)s]",
        #     dest="hmm", choices=["No", "Yes"], metavar="str")
        # parser.add_argument("-m", "--hmm_param",
        #     help="Hidden Markov Model parameters[Default: %(default)s]",
        #     nargs=7, default=[5,100,-123432345,10.5,10.0,'Centre',2], metavar="str", dest="hmmParam")
        parser.add_argument("-T", "--truephase",
            help="True Phase File", metavar="file", type=file, dest='truePhase')
        parser.add_argument("-V", "--version", action="version", version=program_version_message)


        # Process arguments
        args = parser.parse_args()

        outputFile = args.output
        pedigreeFile = args.pedigree
        genotypeFile = args.genotype
        genoFormat = args.genoFormat
        nSNPs = args.nSNPs
        tailsLength = args.tails
        coresLength = args.cores
        useOffset = args.offset
        nSurrogates = args.surrogates
        PercentageSurrDisagree = args.PercentageSurrDisagree
        PercentageGenoHaploDisagree = args.PercentageGenoHaploDisagree
        GenotypeMissingErrorPercentage = args.GenotypeMissingErrorPercentage
        NrmThresh = args.NrmThresh
        useFullOutput = args.fullOutput
        useGraphics = args.graphics
        useSimulation = args.simulation
        # hmm = args.hmm
        # hmmParameters = args.hmmParam
        truePhaseFile = args.truePhase
        

        # Construct file
        if pedigreeFile:
            spec= 'PedigreeFile\t\t\t\t,"{0}"\n'.format(pedigreeFile.name)
        else:
            spec= 'PedigreeFile\t\t\t\t,NoPedigree\n'
        spec+= 'GenotypeFile\t\t\t\t,"{0}",{1}\n'.format(genotypeFile.name,genoFormat)
        spec+= 'NumberOfSnp\t\t\t\t,{0}\n'.format(nSNPs)
        spec+= 'GeneralCoreAndTailLength\t\t,{0}\n'.format(tailsLength)
        spec+= 'GeneralCoreLength\t\t\t,{0}'.format(coresLength)
        if useOffset:
            spec+= ',{0}\n'.format('Offset')
        else:
            spec+= ',{0}\n'.format('NotOffset')
        spec+= 'UseThisNumberOfSurrogates\t\t,{0}\n'.format(nSurrogates)
        spec+= 'PercentageSurrDisagree\t\t\t,{0}\n'.format(PercentageSurrDisagree)
        spec+= 'PercentageGenoHaploDisagree\t\t,{0}\n'.format(PercentageGenoHaploDisagree)
        spec+= 'GenotypeMissingErrorPercentage\t\t,{0}\n'.format(GenotypeMissingErrorPercentage)
        spec+= 'NrmThresh\t\t\t\t,{0}\n'.format(NrmThresh)
        if useFullOutput:
            spec+= 'FullOutput\t\t\t\t,1\n'
        else:
            spec+= 'FullOutput\t\t\t\t,0\n'
        if useGraphics:
            spec+= 'Graphics\t\t\t\t,1\n'
        else:
            spec+= 'Graphics\t\t\t\t,0\n'
        if useSimulation:
            spec+= 'Simulation\t\t\t\t,1\n'
        else:
            spec+= 'Simulation\t\t\t\t,0\n'
        # if hmm=="Yes":
        #     spec+= 'ManageHMM\t\t\t\t,Yes'
        # else:
        #     spec+= 'ManageHMM\t\t\t\t,No'

        # for param in hmmParameters:
        #     spec+= ',' + str(param)
        # spec+= '\n'
        if truePhaseFile is not None and os.path.isfile(truePhaseFile.name):
            spec+= 'TruePhaseFile\t\t\t\t,{0}\n'.format(truePhaseFile.name)
        else:
            spec+= 'TruePhaseFile\t\t\t\t,None\n'


        # Write and close file
        # specFile = open(outputFile, 'w')
        with open(outputFile, 'w') as specFile:
            specFile.write(spec)

        specFile.close()

        # Bye bye
        print outputFile + " has been created successfully.\n"

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help\n")
        return 2

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'argparse_module_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
