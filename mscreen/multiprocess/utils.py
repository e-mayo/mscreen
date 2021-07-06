# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 23:43:31 2021

@author: o_o
"""

import multiprocessing



def count_processors(num_inputs, num_processors):
    """
    Checks processors available and returns a safe number of them to
    utilize.

    :param int num_inputs: The number of inputs.
    :param int num_processors: The number of desired processors.

    :returns: The number of processors to use.
    """
    # first, if num_processors<= 0, determine the number of processors to
    # use programatically
    if num_processors<= 0:
        num_processors = multiprocessing.cpu_count()

    # reduce the number of processors if too many have been specified
    if num_inputs < num_processors:
        num_processors = num_inputs

    return num_processors


def worker(input, output):
    for seq, job in iter(input.get, 'STOP'):
        func, args = job
        result = func(*args)
        ret_val = (seq, result)
        output.put(ret_val)


def start_processes(inputs, num_processors):
    """
    Creates a queue of inputs and outputs
    """

    # Create queues
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()

    # Submit tasks
    for item in inputs:
        task_queue.put(item)

    # Start worker processes
    for i in range(num_processors):
        multiprocessing.Process(target = worker, args = (task_queue, done_queue)).start()

    # Get and print results
    results = []
    for i in range(len(inputs)):
        results.append(done_queue.get())

    # Tell child processes to stop
    for i in range(num_processors):
        task_queue.put('STOP')

    results.sort(key = lambda tup: tup[0])

    return  [item[1] for item in map(list, results)]

def multi_threading(inputs, num_processors, task_name):
    """Initialize this object.

    Args:
        inputs ([data]): A list of data. Each datum contains the details to
            run a single job on a single processor.
        num_processors (int): The number of processors to use.
        task_class_name (class): The class that governs what to do for each
            job on each processor.
    """

    results = []

    # If there are no inputs, just return an empty list.
    if len(inputs) == 0:
        return results

    num_processors = count_processors(len(inputs), num_processors)

    tasks = []

    for index, item in enumerate(inputs):
        if not isinstance(item, tuple):
            item = (item,)
        task = (index, (task_name, item))
        tasks.append(task)

    if num_processors == 1:
        for item in tasks:
            job, args = item[1]
            output = job(*args)
            results.append(output)
    else:
        results = start_processes(tasks, num_processors)

    return results

#%%
from pathlib import Path
import os
import random
from multiprocessing import Pool    

def new_out():
    _id = random.randint(0,100)
    return f'./out-vina-multiprocessing_{_id:0>3d}.pdbqt'

def run_dock_multithread(docking_object, out):
    """
    Run the docking of a single molecule.

    Inputs:
    :param object docking_object: the class for running the chosen docking
        method
    :param str pdb: the path to the pdb of a molecule

    Returns:
    :returns: list failed_smiles_names: any smiles which were deleted (ie.
        docking failed)
    """

    print("Attempt to Dock complete: ", out)
    failed_smiles_names = docking_object.run_dock(out)
    return failed_smiles_names

class docking:
    def run_dock(out):        
        ligand = Path('../../data/prepare/prepared_ligands_vina/lig_3b2x.pdbqt')
        receptor = Path('../../data/prepare/prepared_receptors_vina/receptor_3b2x.pdbqt')
        conf = Path('../../data/config_vina_ex1.txt')
        
        
        torun = (f"vina --config {conf}\
                        --receptor {receptor}\
                        --ligand {ligand}\
                        --out {out}\
                        --cpu 1")
        return os.system(torun)




# if __name__ == '__main__':
        
    # outs = tuple([tuple([docking,new_out()]) for i in range(5)])
    # run_dock_multithread(*outs[0])
    # multi_threading(outs,4,run_dock_multithread)
    # with Pool(5) as p:
        # print(p.map(new_process, outs))



