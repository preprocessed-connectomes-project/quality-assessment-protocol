#!/usr/bin/env python
import argparse
import os
#import nibabel
#import numpy
#from glob import glob
from subprocess import Popen, PIPE
from shutil import rmtree
import subprocess
import yaml
import CPAC.utils as cpac_utils

print("hello!")


def run(command, env={}):
    process = Popen(command, stdout=PIPE, stderr=subprocess.STDOUT,
        shell=True, env=env)
    while True:
        line = process.stdout.readline()
        line = str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break

parser = argparse.ArgumentParser(description='C-PAC Pipeline Runner')
parser.add_argument('bids_dir', help='The directory with the input dataset '
    'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the output files '
    'should be stored. If you are running group level analysis '
    'this folder should be prepopulated with the results of the'
    'participant level analysis.')
parser.add_argument('analysis_level', help='Level of the analysis that will '
    ' be performed. Multiple participant level analyses can be run '
    ' independently (in parallel) using the same output_dir.',
    choices=['participant', 'group'])
parser.add_argument('--participant_label', help='The label of the participant'
    ' that should be analyzed. The label '
    'corresponds to sub-<participant_label> from the BIDS spec '
    '(so it does not include "sub-"). If this parameter is not '
    'provided all subjects should be analyzed. Multiple '
    'participants can be specified with a space separated list.', nargs="+")
parser.add_argument('--pipeline_file', help='Name for the pipeline '
    ' configuration file to use',
    default="/cpac_resources/default_pipeline.yaml")
parser.add_argument('--n_cpus', help='Number of execution '
    ' resources available for the pipeline', default="1")
parser.add_argument('--mem', help='Amount of RAM available to the pipeline'
    '(GB).', default="6")
parser.add_argument('--save_working_dir', action='store_true',
    help='Save the contents of the working directory.', default=False)

# get the command line arguments
args = parser.parse_args()

print(args)
# get and set configuration
c = cpac_utils.Configuration(yaml.load(open(\
    os.path.realpath(args.pipeline_file), 'r')))

# set the parameters using the command line arguements

# we will need to check that the directories exist, and
# make them if they do not
c.outputDirectory = os.path.join(args.output_dir, "output")
c.workingDirectory = os.path.join(args.output_dir, "working")
c.crashLogDirectory = os.path.join(args.output_dir, "crash")
c.logDirectory = os.path.join(args.output_dir, "log")

c.memoryAllocatedPerSubject = int(args.mem)
c.numCoresPerSubject = int(args.n_cpus)
c.numSubjectsAtOnce = 1
c.num_ants_threads = min(args.n_cpus,4)
if( args.save_working_dir == True ):
    c.removeWorkingDir = False
else:
    c.removeWorkingDir = True

print ("#### Running C-PAC on %s"%(args.participant_label))
print ("Number of subjects to run in parallel: %d"%(c.numSubjectsAtOnce))
print ("Output directory: %s"%(c.outputDirectory))
print ("Working directory: %s"%(c.workingDirectory))
print ("Crash directory: %s"%(c.crashLogDirectory))
print ("Log directory: %s"%(c.logDirectory))
print ("Remove working directory: %s"%(c.removeWorkingDir))
print ("Available memory: %d (GB)"%(c.memoryAllocatedPerSubject))
print ("Available threads: %d"%(c.numCoresPerSubject))
print ("Number of threads for ANTs: %d"%(c.num_ants_threads))

# read in the directory to find the input files



#
#
#subjects_to_analyze = []
## only for a subset of subjects
#if args.participant_label:
    #subjects_to_analyze = args.participant_label
## for all subjects
#else:
    #subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
    #subjects_to_analyze = \
        #[subject_dir.split("-")[-1] for subject_dir in subject_dirs]
#
## running participant level
#if args.analysis_level == "participant":
    ## find all T1s and skullstrip them
    #for subject_label in subjects_to_analyze:
        ## grab all T1s from all sessions
        #input_args = " ".join(["-i %s"%f for f in \
            #glob(os.path.join(args.bids_dir,"sub-%s"%subject_label,"anat",
                #"*_T1w.nii*")) + \
            #glob(os.path.join(args.bids_dir,"sub-%s"%subject_label,"ses-*",
                #"anat", "*_T1w.nii*"))])
        #cmd = "echo 'CPAC participant analysis: %s %s %s'"%(subject_label,
            #args.output_dir, input_args)
        #print(cmd)
        #if os.path.exists(os.path.join(args.output_dir, subject_label)):
            #rmtree(os.path.join(args.output_dir, subject_label))
        #run(cmd)
#
   ## Import packages
    #import commands
    #commands.getoutput('source ~/.bashrc')
    #import yaml
    #
    #
    ## Try and load in the subject list
    #try:
        #sublist = yaml.load(open(os.path.realpath(subject_list_file), 'r'))
    #except:
        #raise Exception ("Subject list is not in proper YAML format. Please check your file")
    #
    ## Grab the subject of interest
    #sub_dict = sublist[int(indx)-1]
    #sub_id = sub_dict['subject_id']
#
    #try:
        ## Build and run the pipeline
        #prep_workflow(sub_dict, c, pickle.load(open(strategies, 'r')), 1, p_name, plugin=plugin, plugin_args=plugin_args)
    #except Exception as e:
        #print 'Could not complete cpac run for subject: %s!' % sub_id
        #print 'Error: %s' % e
#
#
#
#elif args.analysis_level == "group":
    ## running group level
    ## generate study specific template
    #cmd = "echo 'CPAC group analysis " +  " ".join(subjects_to_analyze) + "'"
    #print(cmd)
    #run(cmd, env={"SUBJECTS_DIR": args.output_dir})
