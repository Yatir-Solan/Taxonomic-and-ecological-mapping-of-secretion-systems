# modules import:
import io
import os
import re
import string
import random
import datetime
import argparse
from time import sleep
import subprocess as sp

def job_submit(script, job_name='job_runner', cpu=4, queue='duduheavy', run_after=None, wait=60):
    if run_after:
        wait_for_job_to_end(run_after, wait) # sending the submition to hold until the given job/s is ending.

    with open(f'{job_name}_{random.randint(0, 10000000)}', 'w') as template: 
        for line in get_template_file_lines(job_name, cpu, queue): 
            template.write(line)
            template.write('\n')
        template.write(script)
        template.write('\n')

    # This all section was produced due to a problem for sending jobs to the cluster I have encountered.
    return_code = True
    job_stdout = None
    while return_code != 0: 
        if job_stdout:
            print("bad_jobID")
            sleep(5)
        job_stdout = sp.run(f'qsub {template.name}', shell=True, capture_output=True, text=True)
        return_code = job_stdout.returncode
        return_code = False

    os.remove(template.name)
    job_id = job_stdout.stdout.split('.',1)[0] # job_id is returned, as 3567546.power9.tau.ac.il -> 3567546
    return job_id
            
def get_template_file_lines(job_name, cpu, queue):
    return ['#!/bin/bash',
           f'#PBS -q {queue}',
           f'#PBS -N {job_name}',
           f'#PBS -l ncpus={cpu}',
            '#PBS -V',
            '#PBS -j oe',
            'cd $PBS_O_WORKDIR',
            'module load python/python-anaconda3.2019.10']

def is_job_alive(job_id):
    job_id = job_id if type(job_id) is list else [job_id]
    return any(str(jb) in running_jobs() for jb in job_id)

def wait_for_job_to_end(job_id_to_run_after, wait):
    while True:
        if is_job_alive(job_id_to_run_after):
            sleep(wait)
            continue
        break

def running_jobs():
    qstat_output = sp.run('qstat', shell=True, capture_output=True, text=True).stdout
    running_jobs_ids = [job_line.split('.power', 1)[0] for job_line in filter(lambda x:'power' in x, io.StringIO(qstat_output))]
    return running_jobs_ids

def report_file_delete(id_name_dic):
    for job_id,job_name in id_name_dic.items():
        os.remove(f'{job_name}.o{job_id}')

def kill_user_jobs(username):
    print(username)
    qstat_output = sp.run(f'qstat -u {username}', shell=True, capture_output=True, text=True).stdout
    running_jobs_ids = [job_line.split('.power', 1)[0] for job_line in filter(lambda x:'power' in x, io.StringIO(qstat_output))]
    for job_id in running_jobs_ids:
        sp.run(f'qdel {job_id}', shell=True, capture_output=True, text=True)

     
if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-S', '--script', type=str, required=True, help='the script string')
    argparse.add_argument('-N', '--job_name', type=str, default='job_runner', required=False, help='job name')
    argparse.add_argument('-C', '--cpu', type=str, default=4, help='number of CPUs (default=4)')
    argparse.add_argument('-Q', '--queue', type=str, choices=['dudu24h', 'duduheavy', 'duduhimem', 'dudulight', 'duduguest'], default='dudu24h', help='choose the desired queue (default=dudu24h)')
    argparse.add_argument('-R', '--run_after', type=str, default=None, help='start the job run only after a given one id is finshed (single one - 123 / pythonic list - [123,456,789])')
    argparse.add_argument('-W', '--wait_time', type=int, default=60, help="time intervals checks for a 'run_after' to be finshed")
    input_details = argparse.parse_args()

    script = '\n'.join(input_details.script.split('\\n'))
    job_name = input_details.job_name
    cpu = input_details.cpu
    queue = input_details.queue
    run_after = re.sub(r'[\[\]\s]', '', input_details.run_after).split(',') if input_details.run_after else input_details.run_after #converting pythonic string that looks like list to be a real list : '[123,456,789]' -> ['123', '456', '789'] or '123' -> ['123']
    wait_time = input_details.wait_time

    job_submit(script, job_name, cpu, queue, run_after, wait_time)