
import argparse
import random
import sys
import queue
import multiprocessing

def shuffle_seqs_np(fast_f, seq_lines, out_f):
    seq_breaks = np.array([0], dtype=int)
    print("Reading in sequence byte positions...")
    with open(fast_f, 'r') as IN, open(out_f, 'w') as OUT:
        # read in the sequence start positions as a series of byte positions
        while True:
            try:
                line = IN.readline()
            except:
                break

            for i in range(seq_lines - 1):
                IN.readline()

            seq_breaks = np.append(seq_breaks, IN.tell())
            if not len(seq_breaks) % 1000:
                print(len(seq_breaks))

        # print randomly
        print("Printing sequences in random order...")
        length = len(seq_breaks)
        while length:
            indx = random.random_int(0, length - 1)
            
            start = seq_breaks[indx]

            to_write = ''
            for i in range(seq_lines):
                IN.seek(start)
                to_write += IN.readline()
            OUT.write(to_write)

            # remove the index from the array
            seq_breaks = np.delete(seq_breaks, indx)
            length -= 1

        sys.exit()
            
            
        # randomize the start byte positions
        randomize(seq_breaks)
        
        # print the seqs in the new order
        for start in seq_breaks:
            to_write = ''
            IN.seek(start)
            for i in range(seq_lines):
                to_write += IN.readline()
            OUT.write(to_write)

def shuffle_seqs2(fast_f, seq_lines, out_f):
    seq_breaks = []
    print("Reading in sequence byte positions...")
    with open(fast_f, 'r') as IN, open(out_f, 'w') as OUT:
        # read in the sequence start positions as a series of byte positions
        while True:
            try:
                line = IN.readline()
                assert line != ""
            except:
                break

            for i in range(seq_lines - 1):
                IN.readline()

            seq_breaks.append(IN.tell())
            if not len(seq_breaks) % 100000:
                print(len(seq_breaks))

        # print randomly
        print("Printing sequences in random order...")
        length = len(seq_breaks)
        while length:
            indx = random.randint(0, length - 1)
            
            start = seq_breaks.pop(indx)

            to_write = ''
            IN.seek(start)
            for i in range(seq_lines):
                to_write += IN.readline()
            OUT.write(to_write)

            length -= 1


def randomize(array, swaps=None):
    """ Randomize an array by swapping indices """
    length = len(array)
    if not swaps:
        swaps = length * 2

    for i in range(swaps):
        x1 = random.randint(0, length - 1)
        x2 = random.randint(0, length - 1)
        array[x1], array[x2] = array[x2], array[x1]

    return array


class FastxIter(object):

    def __init__(self, fastx_f, lines_per_seq, reverse=False, blocksize=4096):
        self.fh = open(fastx_f, 'r')
        self.lines_per_seq = lines_per_seq
        self.reverse = reverse
        self.blocksize = blocksize

        self.finished = False
        self.queue = queue.Queue(maxsize=1000)


    def get_seqs(self, num_seqs):
        seqs = []
        for i in range(num_seqs):
            try:
                seqs.append(self.queue.get(block=False, timeout=10))
            except queue.Empty:
                if self.finished:
                    if seqs:
                        return seqs
                    else:
                        raise GeneratorExit("All seqs processed")
        #print(seqs)
        return seqs


    def process_seqs(self):
        """ Process seqs in chunks and adds them to the queue. Should be threaded. """
        
        if self.reverse:
            #print("reverse")
            leftovers = []
            for chunk in self._iter_reverse():
                # split the chunk into lines and add the leftover lines from
                # from the last chunk
                lines = chunk.split("\n")
                lines += leftovers

                seqs = int(len(lines) / self.lines_per_seq)

                # calc lines in the begining of the chunk that are not a
                # complete sequence
                orphan_lines = len(lines) % self.lines_per_seq
                
                dup_check = []
                # add each seq to the queue as a single string (reverse order)
                for i in reversed(range(seqs)):
                    # get the base start and end of the sequence
                    start, end = i*self.lines_per_seq, (i+1)*self.lines_per_seq
                    
                    # if there are orphan lines, need to add them to the beginning
                    # and subtract 1 to offset the base 1 or orphan lines
                    if orphan_lines:
                        start, end = start + orphan_lines -1, end + orphan_lines -1
                    
                    self.queue.put("\n".join(lines[start:end]))
                #print("Added {} reverse seqs to queue.".format(str(seqs)))

                # add lines that weren't a full seq to leftovers
                leftovers = lines[0:orphan_lines]

        else:
            leftovers = []
            #print("forward")
            for chunk in self._iter_forward():
                # split the chunk into lines and add the leftover lines from
                # the last chunk
                lines = chunk.split("\n")
                lines += leftovers

                seqs = int(len(lines) / self.lines_per_seq)

                # add each seq to the queue as a single line
                for i in range(seqs):
                    start, end = i*self.lines_per_seq, (i+1)*self.lines_per_seq

                    self.queue.put("\n".join(lines[start:end]))
                #print("Added {} forward seqs to queue.".format(str(seqs)))
        
                # add lines that weren't a full seq to leftovers
                leftovers = lines[seqs*self.lines_per_seq:]

        self.finished = True

    def _iter_reverse(self):
        self.fh.seek(0, 2)
        position = self.fh.tell()

        prev_line = ''
        while position > 0:
            # get either the full block or whatever is left
            delta = min(self.blocksize, position)
            position -= delta
            self.fh.seek(position, 0)
            
            # first line may not be a full line so split it off
            chunk = self.fh.read(delta)
            try:
                first_line, rest = chunk.split("\n", 1)
                yield rest + prev_line
                prev_line = first_line

            except ValueError:      # entire chunk has no newline
                first_line = chunk.split("\n", 1)
                prev_line = first_line[0] + prev_line
        
        # yield the first line in the file
        yield prev_line

    def _iter_forward(self):
        position = 0

        self.fh.seek(0, 2)
        total_size = self.fh.tell()

        self.fh.seek(0, 0)
        prev_line = ''
        while position < total_size:
            delta = min(self.blocksize, total_size - position)
            position += delta
            # last line may not be a full line so split it off
            chunk = self.fh.read(delta)
            try:
                rest, last_line = chunk.rsplit("\n", 1)
                yield prev_line + rest
                prev_line = last_line

            except ValueError:      # entire chunk has no newline
                last_line = chunk.rsplit("\n", 1)
                prev_line += last_line[0]

        # yield the last line in the file
        yield prev_line



def shuffle_seqs(fast_f, lines_per_seq, num_subfiles=6):
    # split fast_f into subfiles
    out_handles = [open("shuf_subf{}.fastx".format(subf), 'w') for subf in range(num_subfiles)]
    
    with open(fast_f, 'r') as IN:
        seqs_remaining = 0
        lines = []
        while True:
            # pick number of seqs to write
            if seqs_remaining == 0:
                seqs_remaining = random.randint(1, 50)


            # collect the lines from one sequence
            try:
                line = IN.readline()
                assert line != ""
                lines.append(line)
            except:
                break

            for i in range(lines_per_seq - 1):
                lines.append(IN.readline())
            
            seqs_remaining -= 1

            # check if collected enough seqs
            if seqs_remaining == 0:
                # pick a random file handle and write the seqs to it
                fh = random.randint(0, num_subfiles - 1)
                out_handles[fh].write("".join(lines))
                lines = []

        # write the last chunk (if needed)
        if lines:
            fh = random.randint(0, num_subfiles - 1)
            out_handles[fh].write("".join(lines))

    # close the open handles
    for fh in out_handles:
        fh.close()

def rewrite_fasta(lines_per_seq, fasta):
    f_iter = FastxIter(fasta, lines_per_seq, reverse=False)
    f_iter_r = FastxIter(fasta, lines_per_seq, reverse=True)

    f_iter.process_seqs()
    f_iter_r.process_seqs()

    iters = [f_iter, f_iter_r]
    

    OUT = open("fwd.fasta", 'w')
    first = True
    while f_iter:
        num_seqs = 1
        try:
            to_write = f_iter.get_seqs(num_seqs)
        except GeneratorExit:   # remove the iter if it is done
            print("Queue finished, removing.")
            f_iter = False
     
        #print(to_write)
        if to_write:    # check if there is anything
            if len([line for line in "\n".join(to_write).split("\n")]) % 8:
                #print(to_write)
                print()
            if first:   # check if first, if so don't new line
                OUT.write("\n".join(to_write))
                first = False
            else:
                OUT.write("\n" + "\n".join(to_write))

    OUT.close()


    OUT = open("rev.fasta", 'w')
    first = True
    while f_iter_r:
        num_seqs = random.randint(0, 25)
        try:
            to_write = f_iter_r.get_seqs(num_seqs)
        except GeneratorExit:   # remove the iter if it is done
            print("Queue finished, removing.")
            f_iter_r = False
     
        #print(to_write)
        if to_write:    # check if there is anything
            if len([line for line in "\n".join(to_write).split("\n")]) % 8:
                #print(to_write)
                print()
            if first:   # check if first, if so don't new line
                OUT.write("\n".join(to_write))
                first = False
            else:
                OUT.write("\n" + "\n".join(to_write))



def random_merge(lines_per_seq, num_subfiles, outfile):
   
    iters = []
    for i in range(num_subfiles):
        if random.randint(0, 1):
            f_iter = FastxIter("shuf_subf{}.fastx".format(str(i)), lines_per_seq, reverse=True)
        else:
            #print(i)
            f_iter = FastxIter("shuf_subf{}.fastx".format(str(i)), lines_per_seq, reverse=False)


        f_iter.process_seqs()
       
        iters.append(f_iter)

        #p = multiprocessing.Process(target=f_iter.process_seqs())
        #p.start()
       
        #processes.append(p)

    with open(outfile, 'w') as OUT:
        first = True
        while iters:
            iter_indx = random.randint(0, len(iters) - 1)
            num_seqs = random.randint(0, 25)
            try:
                to_write = iters[iter_indx].get_seqs(num_seqs)
                print("Wrote {} seqs from iter {}".format(str(num_seqs), str(iter_indx)))
            except GeneratorExit:   # remove the iter if it is done
                print("Queue finished, removing.")
                
                iters.pop(iter_indx)
        
        #print(to_write)
            if to_write:    # check if there is anything
                if len([line for line in "\n".join(to_write).split("\n")]) % 8:
                    #print(to_write)
                    print()
                if first:   # check if first, if so don't new line
                    OUT.write("\n".join(to_write))
                    first = False
                else:
                    OUT.write("\n" + "\n".join(to_write))

    #for p in processes:
        #p.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="fastx file to shuffle", required=True)
    parser.add_argument("-n", help="number of lines per sequence; 2 for fasta, 4 for paired fasta, 4 for fastq, 8 for paired fastq", required=True, type=int)
    parser.add_argument("-s", help="number of subfiles", type=int, default=6)
    parser.add_argument("-o", help="name for out file", required=True)
    args = parser.parse_args()

    #shuffle_seqs(args.f, args.n, args.s)
    rewrite_fasta(args.n, args.f)

    #random_merge(args.n, args.s, args.o)
