from __future__ import print_function
import sys
import getopt
import mmap
import re
import numpy as np
import random
import pdb

def clean(ref_file):
    # flag = 1
    with open(ref_file, "r") as f:
        with open("clean_"+ref_file, "w") as g:
            header=f.readline()
            #print(header, file=g)
            for x in f:
                x = x.rstrip()
                if not x: continue

                # remove leading Ns
                # if((x[0] == "N") & flag):
                #    continue

                # remove other Ns
                if(np.where(x == 'N') == len(x)):
                    continue
                # flag = 0
                print(x, end='', file=g)

# for debugging
def checkOrder(data):
    last = -1
    misordered = 0
    #print(range(data[0].size - 1))
    for i in range(data[0].size - 1):
        if(type(data[0,i]) == int):
            #print(i)
            #print(last)
            #print(data[0,i])
            if data[0,i] != (last + 1):
                misordered = 1
                print("MISORDERED")
                print(i)
                print(data[:,i])
                return misordered
            last = last + 1
    return misordered

def SimSV(ref_file, out_file, num_SVs):
    with open(ref_file, "r") as f:

        # string of reference sequence for tumor genome
        ref = f.read()
        # somatic genome (only initial insertions)
        refS = ref

        # temporary, shorten if not real run
        # ref = ref[0:(len(ref)/1000)] # TEMP
        # refS = ref # TEMP

        L = len(ref)
        N = len(ref)
        L0 = len(ref)
        # length tumor genome has changed
        cumIns = 0
        # length difference between somatic and tumor genomes
        somMutL = 0

        # data type in data structure (object is flexible)
        # dt = np.dtype('>H')
        # dtc = np.dtype('a1')
        dt = np.dtype(object)

        # new index, old index, sequence
        data = np.matrix([range(L), range(L), list(ref)], dtype=dt)
        # SV EVENTS
        events = np.matrix([0]*L, dtype=np.int8)

        # 1: initial insertion
        initial_insertion = 1
        # 2: novel insertion
        insertion = 2
        # 3: deletion
        deletion = 3
        # 4: inversion
        inversion = 4

        # for summary SV counts
        summary = [0, 0, 0, 0]
        # buf = min size of novel insertion, 1/2 min distance between nested SV and novel insertion edge
        buf = 10
        # c = counter for # events
        c = 0
        #loop through each SV
        for i in np.arange(num_SVs):

            # determine where event occurs (deleted region for translocations)
            ev_loc0 = int((L0/num_SVs)*(i) + abs(random.randint(buf, L0/num_SVs - 3*buf)) + cumIns)
            ev_loc0S = ev_loc0 - somMutL

            ### first create novel INSERTION
            insert = ""
            # don't want insertion size very small (so make 3*buf)
            insert_size0 = abs(random.randint(4*buf, 10000))
            j = 0
            while(j < insert_size0):
                insert+=random.choice('ATGC')
                j+=1
            ref = ref[:ev_loc0] + insert + ref[ev_loc0:]
            refS = refS[:ev_loc0S] + insert + refS[ev_loc0S:]
                
            # find data index where new index = ev_loc
            ind = int(np.where(data[0,:] == ev_loc0)[1])
            # format data to insert
            to_insert = np.matrix([[int(x + ev_loc0) for x in range(insert_size0)], list('N'*insert_size0), list(insert)], dtype=dt)
            # concatenate data back together
            data = np.concatenate((data[:,0:ind], to_insert, data[:,ind:N]), axis=1)
            for k in range((ind + insert_size0),(N + insert_size0)):
                if(type(data[0,k]) == int):
                    data[0,k] = data[0,k] + insert_size0
                
            # create blank filler matrix to insert into event track list
            to_insert_event = np.matrix(np.zeros(shape=(c+1,insert_size0), dtype=np.int8))
            # combine filler matrix to event matrix so size correct
            events = np.concatenate((events[:,0:ind], to_insert_event, events[:,ind:N]), axis=1)

            # update length of matrices
            N = N + insert_size0
            L = L + insert_size0
            cumIns += insert_size0
            # add new event column to event matrix
            events = np.concatenate((events, np.matrix(np.zeros(shape=(1,N), dtype=np.int8))), axis=0)
            # update new event column with event type at event locations
            events[c+1,ind:(ind+insert_size0)] = initial_insertion
            # update number of events counter
            c = c + 1
            # update SV count for initial insertions
            summary[0] += 1

            ### Type of event within novel insertion
            event_type = abs(random.randint(2, 4))
            # update SV count for that type
            summary[event_type - 1] += 1

            # determine where nested event occurs (within novel inserted region), has buffer of at least buf to edge of novel insertion
            ev_loc = int(abs(random.randint((ev_loc0 + buf), (ev_loc0 + insert_size0 - 2*buf))))

            ### (nested) INSERTION
            if(event_type == insertion):
                insert = ""
                # insert size less than insert size it is contained in
                insert_size = abs(random.randint(2*buf, insert_size0))
                j = 0
                while(j < insert_size):
                    insert+=random.choice('ATGC')
                    j+=1
                ref = ref[:ev_loc] + insert + ref[ev_loc:]
                
                # find data index where new index = ev_loc
                ind = int(np.where(data[0,:] == ev_loc)[1])
                # format data to insert
                to_insert = np.matrix([[int(x + ev_loc) for x in range(insert_size)], list('N'*insert_size), list(insert)], dtype=dt)
                # concatenate data back together
                data = np.concatenate((data[:,0:ind], to_insert, data[:,ind:N]), axis=1)
                for k in range((ind + insert_size),(N + insert_size)):
                    if(type(data[0,k]) == int):
                        data[0,k] = data[0,k] + insert_size

                # create blank filler matrix to insert into event track list
                to_insert_event = np.matrix(np.zeros(shape=(c+1,insert_size), dtype=np.int8))
                # combine filler matrix to event matrix so size correct
                events = np.concatenate((events[:,0:ind], to_insert_event, events[:,ind:N]), axis=1)
                # update length of matrices
                N = N + insert_size
                L = L + insert_size
                cumIns += insert_size
                somMutL += insert_size
                # add new event column to event matrix
                events = np.concatenate((events, np.matrix(np.zeros(shape=(1,N), dtype=np.int8))), axis=0)
                # update new event column with event type at event locations
                events[c+1,ind:(ind+insert_size)] = event_type


            ### DELETION
            if(event_type == deletion):
                # limit deletion size so can't exceed length of initial insertion (completely contained)
                deletion_size = abs(random.randint(buf, min(10000, ev_loc0 + insert_size0 - ev_loc - buf)))
                ref = ref[:ev_loc] + ref[(ev_loc + deletion_size):]

                # find data index where new index = ev_loc
                ind = int(np.where(data[0,:] == ev_loc)[1])
                # declare ind_2 bc might be existing deletion within bounds of this deletion, so need to extend event more than absolute deletion_size
                ind_2 = int(np.where(data[0,:] == (ev_loc + deletion_size))[1])

                # alter indices of deleted region
                data[0,ind:ind_2] = 'N'

                for k in range(ind_2, N):
                    if(type(data[0,k]) == int):
                        data[0,k] = data[0,k] - deletion_size                

                # update length of matrices
                L = L - deletion_size
                cumIns -= deletion_size
                somMutL -= deletion_size

                # add new event column to event matrix
                events = np.concatenate((events, np.matrix(np.zeros(shape=(1,N), dtype=np.int8))), axis=0)
                # update new event column with event type at event locations
                events[c+1,ind:ind_2] = event_type


            ### INVERSION
            if(event_type == inversion):
                inversion_size = abs(random.randint(buf, min(10000, ev_loc0 + insert_size0 - ev_loc - buf)))
                inverted = ref[ev_loc:(ev_loc + inversion_size)]
                inverted = inverted[::-1]
                ref = ref[:ev_loc] + inverted + ref[(ev_loc + inversion_size):]
                
                # find data index where new index = ev_loc (deletion)
                ind1 = int(np.where(data[0,:] == ev_loc)[1])
                # declare ind_2 bc might be existing deletion within bounds of this deletion, so need to extend event more than absolute deletion_size
                ind1_2 = int(np.where(data[0,:] == (ev_loc + inversion_size))[1])

                # region in index matrix that will be inverted
                data_to_invert = data[:,ind1:ind1_2]
                # invert region
                data_to_invert[:,:] = data_to_invert[:,::-1]
                # note: sequence and old indices can totally invert, while new indices must be adjusted because of possible N's

                # adjust new indices
                counter = 0
                # minimum new index of inverted region
                m = np.amin(data_to_invert[0,:])
                for k in range(ind1_2 - ind1):
                    # if actual int there, replace with inverted int
                    if(type(data_to_invert[0,k]) == int):
                        data_to_invert[0,k] = m + counter
                        counter += 1
                    # (if not int there, keep as N)

                # alter index matrix to match rearrangement
                data = np.concatenate((data[:,0:ind1], data_to_invert, data[:,ind1_2:N]), axis=1)

                # region in events matrix that will be inverted
                events_to_invert = events[:,ind1:ind1_2]
                # rearrange event matrix so indices still align
                events = np.concatenate((events[:,0:ind1], events_to_invert[:,::-1], events[:,ind1_2:N]), axis=1)

                # add new event column to event matrix
                events = np.concatenate((events, np.matrix(np.zeros(shape=(1,N), dtype=np.int8))), axis=0)
                # update new event column with event type at event locations
                events[c+1, ind1:ind1_2] = event_type

            # update number of events counter
            c = c + 1

        with open(out_file+"_tumor.fa", "wb") as s:
            # write the simulated file to output
            s.write(ref)
            # close the simulated output
            s.close()

        with open(out_file+"_normal.fa", "wb") as s:
            # write the simulated file to output
            s.write(refS)
            # close the simulated output
            s.close()

        with open(out_file+"_INDS.txt", "wb") as s: 
            for j in range(data.shape[1]):
                # write the simulated file to output
                print(str(data[0,j]) + " " + str(data[1,j]) + " " + str(data[2,j]), file=s)
            # close the simulated output
            s.close()

        with open(out_file+"_EVENTS.txt", "wb") as s: 
            for j in range(events.shape[1]):
                line = ""
                for i in range(events.shape[0]):
                    line = line + str(events[i,j]) + " "

                # write the simulated file to output
                print(line, file=s)
            # close the simulated output
            s.close()

        with open(out_file+"_SUMMARY.txt", "wb") as s:
                line0 = "initial_insertion " + str(summary[0])
                line1 = "insertion " + str(summary[1])
                line2 = "deletion " + str(summary[2])
                line3 = "inversion " + str(summary[3])
                line4 = "total " + str(sum(summary))
                print(line0, file=s)
                print(line1, file=s)
                print(line2, file=s)
                print(line3, file=s)
                print(line4, file=s)
                s.close()
                

def main():
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.error, msg:
        print(msg)
        print("for help use --help")
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__)
            sys.exit(0)

    # process arguments
    ref_file = args[0]
    out_file = args[1]
    num_SVs = int(args[2])
    if len(args) > 3:
        is_clean = args[3]
        if(is_clean == 'F'):
            #clean reference file (remove newlines)
            clean(ref_file)
            # call simulator
            SimSV("clean_"+ref_file, out_file, num_SVs)
    else:
        # call simulator
        #SimSV(ref_file, out_file, num_SVs)
        
        # try reading with standard I/O
        SimSV(ref_file, out_file, num_SVs)

if __name__ == "__main__":
    main()


# ref_file = "Homo_sapiens.GRCh37.75.dna.chromosome.11.fa"
# out_file = "Chrom11_sim.fa"

### USAGE
#cmd [practice]: python simNestedSV.py ../Data/normIdealSim/Homo_sapiens.GRCh37.75.dna.chromosome.11.fa Chrom11_sim 10 F
#cmd (already clean): python simNestedSV.py clean_Homo_sapiens.GRCh37.75.dna.chromosome.11.fa Chrom11_sim 10
