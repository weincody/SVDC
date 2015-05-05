#~#~ Script for running cancer genome SV detection python

#!/usr/bin/env python
import sys, getopt
import subprocess

#http://stackoverflow.com/questions/4256107/running-bash-commands-in-python

def bash_command(cmd):
    return subprocess.Popen(['/bin/bash', '-c', cmd])

def main(argv):
    som_ref = ""
    canc_ref = ""
    som_reads = ""
    canc_reads = ""
    if(len(argv) != 8):
        print "Need exactly 4 arguments: \npython normIdealSim.py --s_ref <somatic bowtie reference root> --c_ref <cancer bowtie reference root> --s_reads <somatic reads bam> --c_reads <cancer reads bam>"
        sys.exit()
    try:
        opts, args = getopt.getopt(argv,"h",["s_ref=","c_ref=", "s_reads=", "c_reads="])
    except getopt.GetoptError:
        print "python nestedSimDriver.py --s_ref <somatic bowtie reference root> --c_ref <cancer bowtie reference root> --s_reads <somatic reads bam> --c_reads <cancer reads bam>"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print "python nestedSimDriver.py -s_ref <somatic bowtie reference root> -c_ref <cancer bowtie reference root> -s_reads <somatic reads bam> -c_reads <cancer reads bam>"
            sys.exit()
        elif opt in ("--s_ref"):
            som_ref = arg
        elif opt in ("--c_ref"):
            canc_ref = arg
        elif opt in ("--s_reads"):
            som_reads = arg
        elif opt in ("--c_reads"):
            canc_reads = arg
    print "Somatic reference is ", som_ref
    print "Cancer reference is ", canc_ref
    print "Somatic reads are ", som_reads
    print "Cancer reads are ", canc_reads
    return som_ref, canc_ref, som_reads, canc_reads
        

if __name__ == "__main__":

    som_ref, canc_ref, som_reads, canc_reads = main(sys.argv[1:])

    run_socrates_S = "../Tools/Socrates/Socrates all ../Data/nestedSim/" + canc_reads + " --bowtie2_db ../Data/normIdealSim/" + som_ref
    run_socrates_C = "../Tools/Socrates/Socrates all ../Data/normIdealSim/" + som_reads + " --bowtie2_db ../Data/normIdealSim/" + canc_ref
    
    process0 = bash_command(run_socrates_S)
    communicate0 = process0.communicate()
    process1 = bash_command(run_socrates_C)
    communicate1 = process1.communicate()
    
    compile_exportFlanks = "javac BP.java BPs.java exportFlanks.java"
    process2 = bash_command(compile_exportFlanks)
    communicate2 = process2.communicate()

    canc_reads_root = canc_reads[0:len(canc_reads) - 4]
    run_exportFlanks_S = "java exportFlanks ../Tools/Socrates/results_Socrates_paired_" + canc_reads_root + "* S > ../Data/nestedSim/SomFlanks.fa"

    som_reads_root = som_reads[0:len(som_reads) - 4]
    run_exportFlanks_C = "java exportFlanks ../Tools/Socrates/results_Socrates_paired_" + som_reads_root + "* C > ../Data/nestedSim/CancFlanks.fa"

    process3 = bash_command(run_exportFlanks_S)
    communicate3 = process3.communicate()
    process4 = bash_command(run_exportFlanks_C)
    communicate4 = process4.communicate()

    index_S = "../Tools/bowtie2-2.2.3/bowtie2-build ../Data/nestedSim/SomFlanks.fa ../Data/nestedSim/REF_SomFlanks"
    index_C = "../Tools/bowtie2-2.2.3/bowtie2-build ../Data/nestedSim/CancFlanks.fa ../Data/nestedSim/REF_CancFlanks"

    process5 = bash_command(index_S)
    communicate5 = process5.communicate()
    process6 = bash_command(index_C)
    communicate6 = process6.communicate()
    
    map_CtoS = "../Tools/bowtie2-2.2.3/bowtie2 -x ../Data/nestedSim/REF_SomFlanks -c -f ../Data/nestedSim/CancFlanks.fa > ../Data/nestedSim/MAP_C_to_S.txt"
    map_StoC = "../Tools/bowtie2-2.2.3/bowtie2 -x ../Data/nestedSim/REF_CancFlanks -c -f ../Data/nestedSim/SomFlanks.fa > ../Data/nestedSim/MAP_S_to_C.txt"
    
    process7 = bash_command(map_CtoS)
    communicate7 = process7.communicate()
    process8 = bash_command(map_StoC)
    communicate8 = process8.communicate()

    compile_exportMatches = "javac SAM.java SAMs.java exportMatches.java"
    process9 = bash_command(compile_exportMatches)
    communicate9 = process9.communicate()

    run_exportMatches = 'java exportMatches ../Data/nestedSim/MAP_C_to_S.txt ../Data/nestedSim/MAP_S_to_C.txt'
    process10 = bash_command(run_exportMatches)
