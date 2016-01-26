
import inspect
import dill as pickle
import tempfile
import os
import sys
import subprocess
import time
from mypyli.samparser import SamRecord

import logging
logging.basicConfig()
LOG = logging.getLogger()
LOG.setLevel("INFO")


class NotFinishedError(ValueError):
    """ Custom error to raise when the job is not finished """


class LSFParallelizer(object):
    """ Parallelizes a command out over multiple nodes """

    def __init__(self):
        pass

    def run_command(self, command):
        pass

    @classmethod
    def wait_for_job(cls, job_name):
        """ waits for job to complete, checks every 10 seconds """
        while cls._job_running(job_name):
            time.sleep(10)

    @staticmethod
    def _job_running(job_name):

        output = subprocess.check_output([
                    "bjobs",
                    "-J", "{job_name}".format(job_name=job_name)
                ])  

        if output:
            return True
        else:
            return False

class Parallelizer(object):
    """ Spawns a function on another LSF node """

    devnull = open("/dev/null", 'w')

    def __init__(self, function, nodes, processors_per_node, imports=[], output_dir="paral_results", job_prefix="job"):

        # set function and optional imports
        self.function = function
        self.imports = imports

        # set some node power values
        self.nodes = nodes
        self.processors_per_node = processors_per_node
    

        self.output_dir = output_dir

        # make directory if doesn't exist -- subject to unlikely race condition, leaving this for simplicity
        if os.path.isdir(output_dir):
            raise ValueError("Output directory '{}' already exists. Please delete it or supply a new output directory.".format(output_dir))
        else:
            os.mkdir(output_dir)

       
        self.jobs = {}
        self.job_prefix = job_prefix
        self.job_id = 1

        
        # pickle the function 
        self.function_pkl = self.output_dir + "/" + "function.pkl"
        with open(self.function_pkl, 'wb') as OUT:
            pickle.dump(function, OUT)


    def run(self, job_args):
        """ 
        Submits a dict of arguments to the cluster to run through the function as soon as possible.
        
        Returns a batch id to use to check status of the job.
        """
        if type(job_args) != dict:
            raise ValueError("job_args argument must be a dict of names parameters to the function supplied on initialization.")

        job_id = self.job_id
        self.job_id += 1
        job_name = self.job_prefix + str(job_id)
            
        args_pkl = self.output_dir + "/" + job_name + "_arguments.pkl"
        results_pkl = self.output_dir + "/" + job_name + "_results.pkl"
        job_script = self.output_dir + "/" + job_name + "_script.py"
        log_base = self.output_dir + "/" + job_name
            
        # pickle the batch
        with open(args_pkl, 'wb') as OUT:
            pickle.dump(job_args, OUT)

        # add the job to the registry
        self.jobs[job_id] = LSFJob(job_id, job_name)
        self.jobs[job_id].results = results_pkl


        # write the script to execute on the compute node
        with open(job_script, 'w') as OUT:
            OUT.write(self.generate_script(args_pkl, results_pkl) + "\n")


        # wait until there is a spot available
        while self.count_jobs_in_queue() >= self.nodes:
            time.sleep(60)

        # run the job
        self.execute_command(job_id, "python {}".format(job_script), log_base)

        return job_id
 
    def get_all_results(self, only_new=False, wait=False):
        """ 
        Yields results of all jobs (or all new jobs)
        
        Optionally blocks until all submitted jobs either complete or fail
        """

        jobs_remaining = self.jobs.values()

        while jobs_remaining:
            jobs_in_next_loop = []
            for job in jobs_remaining:

                # check if job already has finished status and skip if so
                if job.status in ("completed", "exited"):
                    if only_new:
                        continue

                try:
                    yield self.get_results(job.id)
                except NotFinishedError:
                    jobs_in_next_loop.append(job)

            jobs_remaining = jobs_in_next_loop
            jobs_in_next_loop = []

            # block or return
            if wait:
                time.sleep(60)
            else:
                return 




    def get_results(self, job_id, wait=False):
        """ 
        Returns the results of a job or a ValueError if job has failed or a NotFinishedError is job is still in progress
        """

        LOG.debug("Getting result for job id '{}'".format(job_id))

        # make sure the job exists
        try:
            job = self.jobs[job_id]
        except KeyError:
            raise ValueError("No submitted job with id '{}'".format(job_id))


        while True:

            if job.status == "in progress":
                job.update_status()

            # update the status is the previous status was in progress
            if job.status == "in progress":
                if wait:
                    LOG.debug("Waiting for result...")
                    time.sleep(60)
                else:
                    raise NotFinishedError("Job with id '{}' is still in progress.".format(job_id))

            elif job.status == "exited":
                raise ValueError("Job with id '{}' exited with error.".format(job_id))

            elif job.status == "completed":
                with open(job.results, 'rb') as IN:
                    return pickle.load(IN)

            else:
                raise ValueError("Job status '{}' unknown. -- Module error.".format(job.status))


    #
    ## Submitting jobs
    #
    def submit_to_lsf(self, iterable):
        """ Submits jobs up until the maximum node count """
        while self.count_jobs_in_queue() < self.nodes:
            batch = next(iterable)
            batch_name = self.batch_prefix + str(self.batch_number)
            self.batch_number += 1
            
            

            batch_pkl = self.output_dir + "/" + batch_name + "_arguments.pkl"
            results_pkl = self.output_dir + "/" + batch_name + "_results.pkl"
            batch_script = self.output_dir + "/" + batch_name + "_script.py"
            log_base = self.output_dir + "/" + batch_name
            
            # pickle the batch
            with open(batch_pkl, 'wb') as OUT:
                pickle.dump(batch, OUT)

            with open(batch_scrupt, 'w') as OUT:
                OUT.write(self.generate_script(batch_pkl, results_pkl) + "\n")

            self.execute_command(self, batch_name, "python {}".format(batch_script))

            self.submitted_jobs[job_name] = {'results': results_pkl}
            
    def execute_command(self, job_id, python_command, log_base):
        """ Submits an LSF job """
        lsf_command = self.generate_lsf_command(job_id, o=log_base + ".out", e=log_base + ".err")

        full_command = lsf_command + " " + python_command

        subprocess.call(full_command.split(" "))

        self.jobs[job_id].command = full_command
        self.jobs[job_id].status = "in progress"

    def generate_lsf_command(self, job_id, lsf_program="bsub", q="week", R="span[hosts=1]", M="30", o="%J.out", e="%J.err"):
        """ Generates the lsf portion of the job submission """

        job = self.jobs[job_id]

        command = " ".join([
                lsf_program,
                "-q", str(q),
                "-n", str(self.processors_per_node),
                "-R", str(R),
                "-M", str(M),
                "-J", str(job.name),
                "-o", str(o),
                "-e", str(e)
                ])

        job.stdout = o
        job.stderr = e

        return command

    def generate_script(self, arguments_pkl, output_pkl):
        """ Generates a script to run on the other nodes """

        # set the shebang line
        shebang = "#!/usr/bin/env python"

        # generate the import statements
        imports = "\n".join(["import {}".format(i) for i in self.imports + ["dill"]])

        # load the function and arguments 
        # load function

        load_function = "\n".join([
                    "with open('{}', 'rb') as IN:".format(self.function_pkl),
                    "    function = dill.load(IN)"
                    ])

        load_args = "\n".join([
                    "with open('{}', 'rb') as IN:".format(arguments_pkl),
                    "    args = dill.load(IN)"
                    ])


        # run the function
        running = "results = function(**args)"

        # store the results
        storing = "\n".join([
                    "with open('{}', 'wb') as OUT:".format(output_pkl),
                    "    dill.dump(results, OUT)"
                    ])

        script_code = "\n\n".join([shebang, imports, load_function, load_args, running, storing])

        return script_code


    #
    ## Retrieving results
    #
    @classmethod
    def wait_for_job(cls, job_name):
        """ waits for job to complete, checks every 10 seconds """
        while cls._job_running(job_name):
            time.sleep(10)

    @classmethod
    def check_jobs(cls, job_name):

        output = subprocess.check_output([
                    "bjobs",
                    "-J", "{job_name}".format(job_name=job_name)
                ], stderr=cls.devnull).decode(sys.stdout.encoding) 

        return output
   
    def count_jobs_in_queue(self):
        """ Returns the number of jobs still listed in the queue """
        result = self.check_jobs(self.job_prefix + "*")
        
        if result:
            # have one line for each job (skip the header)
            jobs = result.split("\n")[1:]

            # - 1 because there is a header line
            return len(jobs) - 1

        else:
            return 0


class LSFJob(object):

    devnull = open("/dev/null", 'w')

    def __init__(self, id, name):
        self.id = id
        self.name = name
        
        self.status = None


        self.results = None
        self.command = None
        self.stdout = None
        self.stderr = None

    def update_status(self):
        """ Sets the job status by querying bjobs """


        # check if job is still pending/running
        if self._is_in_queue():
            status = "in progress"

        else:
            # wait until stdout is written
            while not os.path.isfile(self.stdout):
                time.sleep(1)
                
            # open stdout to get status
            with open(self.stdout, 'r') as IN:
                for line in IN:
                    if line == "Successfully completed.\n":
                        status = "completed"
                        break
                    elif line.startswith("Exited with exit code"):
                        status =  "exited"
                        break
                else:
                    raise ValueError("Could not find job status in stdout logs.")

        LOG.debug("Updating job '{}' status to '{}'".format(self.id, status))
        self.status = status


    def _is_in_queue(self):
        """ Returns True if the job is in the queue; False otherwise """

        output = subprocess.check_output([
                    "bjobs",
                    "-J", "{job_name}".format(job_name=self.name)
                ], stderr=self.devnull)  

        if output:
            return True
        else:
            return False


class ClassParallelizer(object):
    """ 
    Acts as an interface to break a class out of a program for parallel processing like a rough MPI.
    
    Classes are pickled, loaded for processing in a separate job, repickled, and then loaded back into the main program.
    """
    
    def __init__(self, cls_to_parallelize, processors=1):
        self.cls = cls_to_parallelize
        self.cls_path, self.cls_file = os.path.split(os.path.splitext(inspect.getsourcefile(cls_to_parallelize))[0])
        print((self.cls_path, self.cls_file))

        self.jobs = {}

    def run(self, obj, method):
        """ Method is a string representation of the method. Ex. to call the foo() method of obj bar, method="foo()" """
        # get a tmpfile for this run
        tmp_fh, tmpfile = tempfile.mkstemp()
        
        # pickle the object
        pickle.dump(obj, tmp_fh)

        # close fh
        tmp_fh.close()

        # write the script to execute the class
        script = tmpfile + ".py"
        with open(script, 'w') as OUT:
            OUT.write(self.gen_script_code(tmpfile, method))

        # run the new script

    def run_script(self, script_path, job_name):
        """ Runs a script using parameters stored in the Parellelizer """
        
        command = "bsub -o %J.out -e %J.err -J {}".format(str(job_name))

        if self.processors:
            command += " -n {} -R 'span[hosts=1]'".format(str(self.processors))

        command += " " + script_path

        subprocess.call(command.split(" "))

    def gen_script_codes(self, obj_tmpfile, method):

        script_code = """
import sys
import pickle

# prepend path
sys.path = {cls_path} + sys.path
from {cls_file} import {cls_name}

obj = pickle.load({obj_tmpfile})
obj.{method}
with open({obj_tmpfile}, 'wb') as OUT:
    pickle.dump(obj)

""".format(cls_path=path, cls_file=file, cls_name=name, obj_tmpfile=tmpfile)

        return script_code

    @staticmethod
    def _get_cls_code(cls):
        return inspect.getsource(cls)

     

if __name__ == "__main__":

    cp = ClassParallelizer(SamRecord)
