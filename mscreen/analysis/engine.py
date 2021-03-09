# -*- coding: utf-8 -*-
"""
Created on Sat Jun 06 02:00:33 2020

@author: eduardo
"""

import numpy as np
from rdkit.Chem import AllChem


def calc_rms_matrix(mol_list):
    """
    calculate the rms matrix of a mol list
    intut:
    -----
    mol_list = [n_mols]

    return:
    ------
    np.array(n,n)
    """
    rmsMatrix = np.zeros((len(mol_list), len(mol_list)))
    # print((len(mol_list)))
    for i in range(len(mol_list)):
        for j in range(i):
            # print(i,j)
            rmsMatrix[i, j] = AllChem.GetBestRMS(mol_list[i], mol_list[j])
            rmsMatrix[j, i] = rmsMatrix[i, j]
    return rmsMatrix


def qt_clust(matrix, threshold):
    """
    A rough quality threashold algorith intended to cluster vina_results

    input:
    -----
    matrix: propierty matrix\n
    threashold: the max diference between element in a cluster

    return:
    ------
    qt_clust: list storing the clusters
    """
    # this set store the elements that have been clustered
    already_in_clust = set()
    # the number of elements
    poses = matrix.shape[0]
    # the clusters
    qt_clust = []
    poselist = set(range(poses))
    while poselist != already_in_clust:

        test_cluster = []  # the clusters test
        length_cluster = []  # this store the lengh of each cluster
        for i in range(poses):
            clust = []  # this list store all the molecules that satisfies the condition
            if i in already_in_clust:
                continue
        #    clust.append(i)
            go_out = False  # didn't figured out a better way to do this at 3:00 am lol
            for j in range(poses):
                if j in already_in_clust:
                    continue
                # print(i, j, matrix[i, j])

                # looking for molecules that satisfies the condition
                if matrix[i, j] <= threshold:
                    for a in clust:
                        # each molecule must satisfies the condition to all other molecules
                        if matrix[a, j] > threshold:
                            go_out = True
                    if go_out:
                        break
                    clust.append(j)
            test_cluster.append(clust)
            length_cluster.append(len(clust))
            # break #didn't remember this break where come from probaly from debugging
        max_clust = max(length_cluster)
        max_clust_index = length_cluster.index(max_clust)
        for i in test_cluster[max_clust_index]:
            already_in_clust.add(i)
        qt_clust.append(test_cluster[max_clust_index])
    return qt_clust


def qt_cluster_a_mol_list(mol_list, threashold, prop_matrix=None):
    """
    clusterize a mol list using the qt algorithm over the rms values
    store the cluster_id and the lengh of the cluster in each molecule
    input:
    -----
    mol_list --- molecules list
    threashold --- the quality threashold
    prop_matrix --- the propierties matrix

    return:
    ------
    mol_list --- a mol list with all the cluster stuff
    """
    # print('calculating rmsdMatrix')
    rmsMatrix = calc_rms_matrix(mol_list)
    # print('end of calculating rmsdMatrix')
    cluster = qt_clust(rmsMatrix, threashold)
    for i, j in enumerate(cluster):
        for m in j:
            mol_list[m].SetProp('cluster_id', str(i))
            mol_list[m].SetProp('clust_lenght', str(len(j)))
    return(mol_list)


def get_representative_clust(mol_list, cluster = None):
    """
    Figure out the representative pose. Attempting the criteria:
        max cluster population
        low cluster energy
    """
    cluster_id = [m.GetProp('cluster_id') for m in mol_list]
    clust_lenght = [m.GetProp('clust_lenght') for m in mol_list]
    enrgies = [m.GetProp('REMARK').split()[2] for m in mol_list]

    max_cluster_leght = clust_lenght.index(max(clust_lenght))
    min_energy = enrgies.index(min(enrgies))

    if cluster_id[min_energy] == cluster_id[max_cluster_leght]:
        mol_list[min_energy].SetProp('best_pose', 'is representative')
        return mol_list
    else:
        mol_list[min_energy].SetProp('best_pose', 'not representative best energy')
        mol_list[max_cluster_leght].SetProp('best_pose', 'not representative longest cluster')
        return mol_list

def find_mol_within_site_radial(mol,site,radius,min_atoms):
        c = mol.GetConformer()
        xyz_molecule = c.GetPositions()
        atoms_in_site = 0
        mol.SetProp('in_site', str('False'))
        in_site = False
        for atom_xyz in xyz_molecule:
            distance = ((atom_xyz - site)**2).sum(-1)
            if np.sqrt(distance) < radius:
                atoms_in_site += 1
            if atoms_in_site >= min_atoms:
                mol.SetProp('in_site', str('True'))
                in_site = True
                break
        return in_site


def __old__find_mol_within_site_radial__(mol,site,radius):
        c = mol.GetConformer()
        xyz_molecule = c.GetPositions()
        atoms_in_site = 0
        for atom_xyz in xyz_molecule:
            distance = ((atom_xyz - site)**2).sum(-1)
            if np.sqrt(distance) < radius:
                atoms_in_site += 1


def find_mol_within_site_square(mol,site,side,min_atoms):
            c = mol.GetConformer()
            xyz_molecule = c.GetPositions()
            mol.SetProp('in_site', str('False'))
            atoms_in_site = 0
            in_site = False
            for atom_xyz in xyz_molecule:
                if site[0] - side < atom_xyz[0] < site[0] + side:
                    if site[1] - side < atom_xyz[1] < site[1] + side:
                        if site[2] - side < atom_xyz[2] < site[2] + side:
                            atoms_in_site += 1
                            if atoms_in_site >= min_atoms:
                                mol.SetProp('in_site', str('True'))
                                in_site = True
                                break
            return in_site