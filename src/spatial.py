#! ./venv_ea python
# -*- coding: utf-8 -*-

__author__ = 'Subhodip Biswas'
__email__ = 'subhodip@cs.vt.edu'
__version__ = '0.1'

import math
import json
import copy
import time
from tqdm import tqdm
import multiprocessing as mp
from os.path import join
from random import randrange, random, shuffle
from utilities import make_move, local_move,\
    set_params, gen_solutions, get_solution, get_func_val,\
    update_property, get_neighbors, get_borderareas,\
    create_dir, get_results, get_params, parameters, check_solution
from get_inputs import GetInputs


def prob_selection(weights):
    """
    Probabilistic selection
    :param weights: List of values based on which the selection is made
    :return: The selected list element
    """
    total = 0
    winner = 0
    for i, w in enumerate(weights):
        total += w
        if random() * total < w:
            winner = i

    return winner


def make_swap(args, i, j, zi, zj, solutions):
    """
    Swap nodes between solutions i and j
    :param i: ID for solution i
    :param j: ID for solution j (more fit)
    :param zi: zone ID of a given spatial unit in solution i
    :param zj: zone ID of the same spatial unit in solution j
    :param solutions: Copy of the solutions to the spatial optimization problem
    :return: Update solution if successful swap has happened
    """
    zones1, zone_ids1 = copy.deepcopy(solutions[i]['zones']), copy.deepcopy(solutions[i]['zone_ids'])
    zones2, zone_ids2 = copy.deepcopy(solutions[j]['zones']), copy.deepcopy(solutions[j]['zone_ids'])
    swap, solution = False, solutions[i]
    spas_nbr = args[4]

    try:
        # Find mutually exclusive set of spatial units for swapping
        # Set operation: Zj - Zi
        areas_inzj_notzi = [
            m for m in zones2[zj]['members']
            if m not in zones1[zi]['members']
        ]
        # Set operation: Zi - Zj
        areas_inzi_notzj = [
            m for m in zones1[zi]['members']
            if m not in zones2[zj]['members']
        ]

        # The spatial units should not be located in the center of the respective zones, otherwise it will lead to holes
        neighbors_areas = get_neighbors(args, zi, zones1, zone_ids1)
        border_areas = get_borderareas(args, zi, zones1, zone_ids1)
        # Check for overlap
        incoming_areas = [m for m in areas_inzj_notzi if m in neighbors_areas]    # areas for moving into Zi
        outgoing_areas = [m for m in areas_inzi_notzj if m in border_areas]       # areas for moving out of Zi

        # Simultaneously move an (incoming) area into Zi and an (outgoing) area out of Zj
        if len(incoming_areas) > 0 and len(outgoing_areas):  # both sets should be non-empty
            in_area = incoming_areas[randrange(len(incoming_areas))]

            # Determine the incoming spatial unit
            in_id = zone_ids1[in_area]
            incoming = make_move(args, [in_id, zi], in_area, zones1, zone_ids1)

            # Determine the outgoing spatial unit
            outgoing = False
            c = 0
            while c < 10 and not outgoing:
                out_area = outgoing_areas[randrange(len(outgoing_areas))]
                c += 1
                if in_area == out_area:
                    continue
                else:
                    out_ids = [zone_ids1[a] for a in spas_nbr[out_area] if zone_ids1[a] != zi]
                    if len(out_ids) > 0:
                        out_id = out_ids[randrange(len(out_ids))]
                        outgoing = make_move(args, [zi, out_id], out_area, zones1, zone_ids1)

            swap = incoming and outgoing
    except Exception as e:
        pass

    if swap:
        # print('Swap happened!!')
        update_property(args, zones1)
        solution = get_solution(zones1, zone_ids1)

    return swap, solution


def spatial_recombination(args, solutions, i):
    """
    Perform spatially-aware recombination to modify solution 'i'
    :param args: Tuple for spatial computation (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schlattr)
    :param solutions: Copy of the solutions to the spatial optimization problem
    :param i: ID for solution i
    :return: Updated solution 'i' by recombination operation
    """
    pop_size = len(solutions)
    fitness = [1.0/(1 + abs(solutions[p]['func_val'])) for p in range(pop_size)]    # Fitness function for minimization
    weights = [fitness[p]/sum(fitness) for p in range(pop_size)]
    # Select another solution
    j = i
    while j == i:
        j = prob_selection(weights)

    # Recombination operator using the two solutions
    center_areas = [z for z in args[5].keys()]
    solution = copy.deepcopy(solutions[i])
    swap = False

    while not swap and len(center_areas) > 0:
        area = center_areas[randrange(len(center_areas))]
        # Find the respective zones containing 'area'
        zi, zj = solutions[i]['zone_ids'][area], solutions[j]['zone_ids'][area]

        if zi == zj:
            swap, solution = make_swap(args, i, j, zi, zj, solutions)

        if not swap:
            center_areas.remove(area)

    return solution


def local_improvement(args, solutions, i):
    """
    Perform local improvement to modify solution 'i'
    :param args: Tuple for spatial computation (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schlattr)
    :param solutions: Copy of the solutions to the spatial optimization problem
    :param i: ID for solution i
    :return: Updated solution 'i' by local search
    """
    solution = copy.deepcopy(solutions[i])
    zones, zone_ids = solution['zones'], solution['zone_ids']
    zone_list = [x for x in zones.keys()]
    shuffle(zone_list)
    moved = False

    while not moved and len(zone_list) > 0:
        recip_id = zone_list[randrange(len(zone_list))]
        neighbors = get_neighbors(args, recip_id, zones, zone_ids)

        while not moved and len(neighbors) > 0:
            # Randomly select areas until there is a local improvement
            area = neighbors[randrange(len(neighbors))]
            donor_id = zone_ids[area]  # Define the donor zone
            moved = make_move(args, [donor_id, recip_id], area, zones, zone_ids)
            neighbors.remove(area)  # remove it

        zone_list.remove(recip_id)  # remove it so that new zones can be picked

    if moved:
        # print('Updating fitness')
        update_property(args, zones)
        solution = get_solution(zones, zone_ids)

    return solution


def run_module(args, phase, solutions):
    """
    Run the improvement modules
    :param args: Tuple for spatial computation (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schlattr)
    :param phase: Binary variable - 0 for local improvement, 1 for spatial recombination
    :param solutions: Copy of the solutions to the spatial optimization problem
    :return: Updated solutions
    """
    # Run parallel processing to update the solutions
    pop_size = len(solutions)
    num_process = min(pop_size, mp.cpu_count())
    pool = mp.Pool(processes=num_process)
    modules = [local_improvement, spatial_recombination]
    try:
        output = [
            (p, pool.apply_async(modules[phase],
                                 args=(args, solutions, p)
                                 )
             )
            for p in range(pop_size)
        ]
        # Fitness based selection
        for p, o in output:
            solution = o.get()
            if solution['func_val'] < solutions[p]['func_val']:
                solutions[p] = solution

    except Exception as e:
        print("Error: {} in run_module()!!".format(e))

    pool.close()


def find_best_sol(solutions, best_f, stagnate):
    """
    Determine the best solution in the population.
    :param solutions: Copy of the solutions to the spatial optimization problem
    :param best_f: The functional value of the present best solution (to be updated)
    :param stagnate: Stagnation counter
    :return: New best solution and its fitness value, updated value of stagnation counter
    """
    new_best_f = math.inf
    new_best_sol = None
    for p in range(len(solutions)):
        if solutions[p]['func_val'] < new_best_f:
            new_best_f, new_best_sol = solutions[p]['func_val'], p

    # Check if solutions have converged to a local optimum
    if best_f - new_best_f < pow(10, -4):
        stagnate += 1
    else:
        stagnate = 0

    return new_best_f, new_best_sol, stagnate


def run_spatial(options):
    """
    Simulate runs for the SPATIAL algorithm
    """
    sch = options.school
    init = options.initialization
    runs = options.runs
    algo = options.algo

    print('\n... Starting {} algorithm ...\n'.format(algo))
    # SPATIAL hyper-parameters
    pop_size = 10
    iter_max = 1000

    # Read data files
    inputs = GetInputs(sch)
    args = inputs.get_inputs()

    # Values for the data set
    weight = 7    # Set weight in the range [0, 1] for calculating F = w * F1 + (1 - w) * F2
    set_params(weight, iter_max)
    # init_type = {1: 'seeded', 2: 'infeasible', 3: 'existing'}
    existing = gen_solutions(args, 3)

    print(' Present school boundary configuration has functional value : {:.3f}'.format(existing['func_val']))

    for r in range(runs):
        print('\n Run {} of {}\n'.format(r+1, algo))
        try:
            # Generate random seeds and use them for instantiating initial trial solutions
            seeds = [s + randrange(1000000) for s in range(pop_size)]
            solutions = copy.deepcopy(gen_solutions(args, init, pop_size, seeds))
            # Find the best solution and its functional value
            best_func_val, best_sol, stagnation = find_best_sol(solutions, math.inf, 0)
            print('Iter: 0 \t Best func_val: {:.3f}'.format(best_func_val))

            iteration = 0
            t_start = time.time()
            fval_iter = [(0, best_func_val)]    # List to save best solutions
            time_iter = [(0, time.time() - t_start)]

            # Iteratively improve the solutions using the two search operators
            for it in tqdm(range(iter_max)):
                run_module(args, 0, solutions)    # Local improvement
                run_module(args, 1, solutions)    # Spatially-aware recombination
                iteration = it + 1
                best_func_val, best_sol, stagnation = find_best_sol(solutions,
                                                                    best_func_val,
                                                                    stagnation)
                # Print the best result
                if iteration % 20 == 0:
                    print('Iter: {} \t Best func_val: {:.3f}'.format(iteration, best_func_val))
                    fval_iter.append((iteration, best_func_val))
                    time_iter.append((iteration, time.time() - t_start))

            # Printing the results
            t_elapsed = (time.time() - t_start) / 60.0  # measures in minutes
            best_solution = solutions[best_sol]

            print("Run: {} took {:.2f} min to execute {} iterations...\n"
                  " Obtained FuncVal: {:.3f} ".format(r + 1, t_elapsed, iteration, best_func_val))

            # Save the results
            print(' Checking the correctness of the solution...')
            correct = check_solution(args, best_solution['zones'])
            if correct:
                print('\n Correct solution.. Saving results .... ')
                w1, w2, epsilon, _ = parameters()
                params = get_params(AlgParams=get_params(w1=w1, w2=w2,
                                                         epsilon=epsilon,
                                                         pop_size=pop_size),
                                    Iteration=iteration,
                                    TimeElapsed=t_elapsed,
                                    School=sch)
                solu_info = {'Existing': existing,
                             'Final': best_solution,
                             'fval_vs_iter': fval_iter}
                run_results = {'properties': params,
                               'info': solu_info}
                write_path = create_dir('results', algo, sch)
                with open(join(write_path,
                               "run{}_{}_{}.json".format(r + 1, algo, sch)),
                          'w') as outfile:
                    json.dump(run_results, outfile)
            else:
                print('\n Incorrect solution. Has disconnected zones... \n')

        except Exception as e:
            print('Run {} incomplete due to error: {}'.format(r + 1, e))
