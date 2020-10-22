import os
import copy
import time
import math
import random
import multiprocessing as mp
from shapely.ops import unary_union
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

w1 = 0.0
w2 = 0.0
max_iter = 0
epsilon = pow(10,0)


def print_params():
    print(' w1 (Capacity) : {:.2}\n w2 (Compactness) : {:.2}'.format(w1, w2))
    print(' Epsilon : {}\n MaxIter : {}'.format(epsilon,
                                                max_iter
                                                )
          )


def set_params(w, iter_max = 1000):
    global w1, w2, epsilon, max_iter
    w1 = 0.1*w
    w2 = 1-w1
    epsilon = pow(10, -5)  # global constant
    max_iter = iter_max

    print_params()


def parameters():
    return w1, w2, epsilon, max_iter


def gen_solutions(args, init=1, num_sol=1, seeds=None):
    """Generate solutions using the initialize()"""

    solutions = dict()
    try:
        if num_sol > 0:
            if 0 < init < 3:
                print('\n....Generating {} solutions....\n'.format(num_sol))
                t = time.time()
                # Parallel (asynchronous) solution initialization
                pool = mp.Pool(processes=min(mp.cpu_count() - 1, num_sol))
                output = [(i,
                           pool.apply_async(initialize,
                                            args=(args, init,
                                                  (None, seeds[i])[seeds is not None]
                                                  )
                                            )
                           )
                          for i in range(num_sol)
                          ]
                pool.close()
                for i, p in output:
                    zones, zone_ids = p.get()
                    update_property(args, zones)
                    solution = get_solution(zones, zone_ids)
                    solutions[i] = solution
                '''
                # Serial execution
                for i in tqdm(range(num_sol)):
                    solutions[i] = initialize(init, seeds[i])
                '''
                print('\n Done.. ')
                t_elapsed = time.time() - t
                print('\n Time taken: {:.2} min\n'.format(t_elapsed/60.0))

            if init == 3:
                zones, zone_ids = initialize(args, init)
                update_property(args, zones)
                solutions = get_solution(zones, zone_ids)

    except Exception as e:
        print(e)
        print("Couldn\'t generate solution(s)!!")

    return solutions


def initialize(args, init, seed=None):
    """Initializes zones with different schemes"""

    zones, zone_ids = None, None
    try:
        # 'seeded' initialization
        if init == 1:
            zones, zone_ids = seeded_init(args, seed)
        # 'infeasible' initialization that doesn't satisfy contiguity of zones
        elif init == 2:
            zones, zone_ids = infeasible(args, seed)
        # 'existing' initialization
        elif init == 3:
            zones, zone_ids = exist_init(args)

    except Exception as e:
        print("Error: {} in initialize()".format(e))

    return zones, zone_ids


def seeded_init(args, seed=None):
    """
    Seeded initialization starting with school containing polygons
    """

    zones, zone_ids = dict(), dict()
    try:
        random.seed(seed)
        args = args
        # args = (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schl_attr)
        spas, spas_nbr, sch_spas_nbr, schl_attr = args[3], args[4], args[5], args[6]

        # Enumerate zones and set status
        spas_list = [s for s in spas_nbr.keys()]
        for area in spas_list:
            zone_ids[area] = -1  # -1 means not assigned to a zones

        # Initialize the zones with school-containing polygons
        sch_list = []
        for area in sch_spas_nbr.keys():
            sch_code = spas[schl_attr][area]
            zone_ids[area] = sch_code  # Assign the SCH_CODE to the area
            spas_list.remove(area)      # Remove areas that have already been assigned to a zone
            sch_list.append(sch_code)
            zones[sch_code] = get_params(members=[area], SCH=sch_code)

        num_zones = len(sch_list)     # No. of schools

        while len(spas_list) > 0:
            # Pick a random zones
            zoneid = sch_list[random.randrange(num_zones)]
            members = [x for x in zones[zoneid]['members']]
            neighbors = get_adj_areas(members, spas_nbr, zone_ids)   # Get list of free areas around it

            if len(neighbors) > 0:
                area = neighbors[random.randrange(len(neighbors))]
                zone_ids[area] = zoneid
                zones[zoneid]['members'].append(area)
                spas_list.remove(area)

        if list(zone_ids.values()).count(-1) > 0:
            print('There are unassigned polygons present. Error!!')

    except Exception as e:
        print("{} error in seeding initialization!!".format(e))

    return zones, zone_ids


def exist_init(args):
    """
    Extracts the existing partition for evaluation
    """
    zones, zone_ids = dict(), dict()
    # args = (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schl_attr, zones, zone_ids)

    # Enumerate zones and set status
    spas, spas_nbr, schl_attr = args[3], args[4], args[6]
    spas_list = [s for s in spas_nbr.keys()]

    for area in spas_list:
        zone_ids[area] = -1  # -1 means not assigned to a zones

    # Get the existing partition from the data
    for index, spa in spas.iterrows():
        area = spa['SPA']
        zoneid = spa[schl_attr]
        zone_ids[area] = zoneid

        if zoneid not in zones.keys():
            zones[zoneid] = get_params(members=[], SCH=zoneid)

        zones[zoneid]['members'].append(area)

    return zones, zone_ids


def infeasible(args, seed=None):
    """
    The resultant zones maintain contiguity constraint but might not contain one school per partition.
    """

    zones, zone_ids = dict(), dict()
    random.seed(seed)

    # args = (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schl_attr, zones, zone_ids)
    # Enumerate zones and set status
    spas_nbr = args[4]
    spas_list = [x for x in spas_nbr.keys()]

    for area in spas_list:
        zone_ids[area] = -1  # -1 means not assigned to a zones

    exist_zones, _ = exist_init(args)  # get existing partition

    # Initialize the M zones with M areas (keeping some familiarity with existing partition)
    for zoneid in exist_zones.keys():
        exist_zones = [m for m in exist_zones[zoneid]['members']]
        area = exist_zones[random.randrange(len(exist_zones))]  # pick random area from zones

        zones[zoneid] = dict()
        zones[zoneid]['members'] = [area]
        zone_ids[area] = zoneid  # key is the 'school_name'
        spas_list.remove(area)

    sch_list = [s for s in exist_zones.keys()]
    num_zones = len(sch_list)

    # Put the unassigned polygons in the zones
    while len(spas_list) > 0:
        zoneid = sch_list[random.randrange(num_zones)]  # pick randomly a zones to grow
        members = [m for m in zones[zoneid]['members']]
        neighbors = get_adj_areas(members, spas_nbr, zone_ids)

        if len(neighbors) > 0:
            index = random.randrange(len(neighbors))
            area = neighbors[index]
            zones[zoneid]['members'].append(area)
            zone_ids[area] = zoneid
            spas_list.remove(area)

    return zones, zone_ids


def get_borderareas(args, zid, zones, zone_ids):
    """
    Get the list of areas on the border of the base zone 'zid'
    """

    # args = (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schl_attr)
    spas_nbr, sch_spas_nbr = args[4], args[5]
    centers = [x for x in sch_spas_nbr.keys()]
    border_areas = []

    for area in zones[zid]['members']:
        for x in spas_nbr[area]:
            if zone_ids[x] != zid and area not in centers:
                border_areas += [area]
                break

    border_areas = list(set(border_areas))  # get unique list
    return border_areas


def get_neighbors(args, zid, zones, zone_ids):
    """
    Get the list of areas adjacent to the base zones
    """

    spas_nbr, sch_spas_nbr = args[4], args[5]
    centers = [m for m in sch_spas_nbr.keys()]
    neighbors = []

    for area in zones[zid]['members']:
        neighbors = neighbors + [x for x in spas_nbr[area] if zone_ids[x] != zid and x not in centers]

    neighbors = list(set(neighbors))  # get unique list

    return neighbors


def get_adj_areas(areas, spas_nbr, zone_ids):
    """
    Returns adjacent unassigned area polygons to a cluster
    """
    adj_areas = []
    if len(areas) > 0:
        for area in areas:
            adj_areas = adj_areas + [a for a in spas_nbr[area]
                                     if zone_ids[a] == -1]

        adj_areas = list(set(adj_areas))

    return adj_areas


def get_solution(zs, zids, func_val=None):
    """
    Generate a solution of the problem
    """
    zones = copy.deepcopy(zs)
    zone_ids = copy.deepcopy(zids)

    if func_val is None:
        func_val = get_func_val(zones)

    solution = get_params(zones=zones, zone_ids=zone_ids, func_val=func_val)

    return solution


def check_solution(args, zones, zids=None):
    """
    Checks if all the individual regions are geographically connected.
    """
    if zids is None:
        zids = zones.keys()

    connected = True
    adjacency = args[2]
    for zid in zids:
        zone = [m for m in zones[zid]['members']]
        connected = if_connected(adjacency, zone)
        if not connected:
            print('Solution error since zone {} is not connected!!'.format(zid))
            break

    return connected


def get_results(r, args, solution):
    """
    Get the output results in a dictionary form
    """

    existing = gen_solutions(args, 3)
    print("Run: {}\t  Present FuncVal {:.4f}   Obtained FuncVal: {:.4f}".format(r,
                                                                                existing['func_val'],
                                                                                solution['func_val']
                                                                                )
          )
    solu_info = {'Existing': existing, 'Final': solution}
    return solu_info


def create_dir(*args):
    """
    Function to create directories and generate paths.
    """
    write_path = "../"
    for name in args:
        try:
            write_path =  write_path + "{}/".format(name)
            os.mkdir(write_path)
        except Exception as e:
            print(".")

    return write_path


def get_params(**argsv):
    """
    Returns a dictionary
    """

    return argsv


def get_args(*args):
    """
    Returns a tuple
    """
    return args


def if_connected(adjacency, zone):
    """
    Determines if the spatial units forming a zone are connected or not.

    :param adjacency: Adjacency matrix of the spatial units
    :param zone: Set of spatial units forming the zone
    :return: boolean value 'True' if the zone is connected else 'False'
    """
    connected = False

    if len(zone) > 0:  # If the cluster is not a singleton
        zone_adj = adjacency.loc[zone, zone].values
        adj_mat = csr_matrix(zone_adj)
        num_connect_comp = connected_components(adj_mat,
                                                directed=False,
                                                return_labels=False)

        connected = num_connect_comp == 1

    return connected


def check_move(args, cids, area, zones):
    """
    Check if moving a unit  between the spatial units preserves contiguity
    """
    adjacency = args[2]
    donor_id, recip_id = cids[0], cids[1]
    donor_zone = [m for m in zones[donor_id]['members']]
    recip_zone = [m for m in zones[recip_id]['members']]
    donor_connected, recip_connected = False, False

    try:
        # Move 'area' from donor_zone to recipient_zone
        donor_zone.remove(area)
        recip_zone.append(area)

        donor_connected = if_connected(adjacency, donor_zone)
        recip_connected = if_connected(adjacency, recip_zone)

    except Exception as e:
        pass

    return donor_connected, recip_connected


def repair(args, zids, zones, zone_ids, connecteds):
    """
    Apply the repair operation if the spatial contiguity of a zone if broken
    """
    adjacency, spas_nbr, sch_spas_nbr = args[2], args[4], args[5]
    centers = [x for x in sch_spas_nbr.keys()]
    involved_ids = []

    for z in range(len(zids)):
        zid = zids[z]
        involved_ids.append(zid)

        if not connecteds[z]:
            zone = [m for m in zones[zid]['members']]
            zone_adj = adjacency.loc[zone, zone].values
            adj_mat = csr_matrix(zone_adj)
            num_connect_comp, labels = connected_components(adj_mat, directed=False, return_labels=True)
            assert num_connect_comp > 1    # if it is not connected

            # Perform repair process by scanning each connected component
            for c in range(num_connect_comp):
                connect_sub_zone = [zone[i] for i in range(len(zone)) if labels[i] == c]
                good = False

                for area in connect_sub_zone:
                    if area in centers:
                        good = True    # A 'good' connected component should be kept as it is
                        break

                if not good:
                    # Perform reassignment to neighboring zones
                    for area in connect_sub_zone:
                        nbr_ids = [zone_ids[nbr] for nbr in spas_nbr[area] if zone_ids[nbr] != zid]
                        if len(nbr_ids) > 0:
                            nbr_id = nbr_ids[random.randrange(len(nbr_ids))]
                            # Make the move
                            zones[zid]['members'].remove(area)
                            zones[nbr_id]['members'].append(area)
                            zone_ids[area] = nbr_id
                            involved_ids.append(nbr_id)

    # Make sure that the 'involved ids' are connected
    involved_ids = list(set(involved_ids))
    for zid in involved_ids:
        zone = [m for m in zones[zid]['members']]
        zone_adj = adjacency.loc[zone, zone].values
        adj_mat = csr_matrix(zone_adj)
        num_connect_comp = connected_components(adj_mat, directed=False, return_labels=False)
        assert num_connect_comp == 1, "Disconnected components"    # it is not connected


def make_move(args, zids, area, zones, zone_ids, do_repair=True):
    """
    Moves a spatial unit between two zones
    """
    moved = False
    try:
        donor_id, recip_id = zids[0], zids[1]
        assert donor_id != recip_id, "Donor and recipient zones are same!!"
        assert len(zones[donor_id]['members']) > 1, "Will result in empty \'donor\' zone!!"    # to prevent empty zones
        donor_connected, recipient_connected = check_move(args, zids, area, zones)
        zones[donor_id]['members'].remove(area)
        zones[recip_id]['members'].append(area)
        zone_ids[area] = recip_id
        # Perform repair if either of zones are not connected as a result of the move
        if do_repair and (not donor_connected or not recipient_connected):
            repair(args, zids, zones, zone_ids, [donor_connected, recipient_connected])

        moved = True
    except Exception as e:
        pass    # Exception won't effect final solution

    return moved


def local_move(args, zids, area, zs, z_ids):
    """
    Computes the cost function of a set of zones

    :param zones: Dictionary containing the properties of the zones
    :param zids: List of zones
    :return: Aggregated functional value of the list of zones
    """
    zones, zone_ids = copy.deepcopy(zs), copy.deepcopy(z_ids)
    possible, change_func_val = False, 0
    try:
        donor_id, recip_id = zids[0], zids[1]
        assert donor_id != recip_id, "Donor and recipient zones are same!!"
        assert len(zones[donor_id]['members']) > 1, "Will result in empty \'donor\' zone!!"  # prevent empty zones
        donor_connected, recipient_connected = check_move(args, zids, area, zones)
        possible = donor_connected and recipient_connected

        if possible:
            change_func_val -= sum([zones[zid]['func_val'] for zid in zids])
            zones[donor_id]['members'].remove(area)
            zones[recip_id]['members'].append(area)
            zone_ids[area] = recip_id
            update_property(args, zones, zids)
            change_func_val += sum([zones[zid]['func_val'] for zid in zids])
    except Exception as e:
        print('{} exception in local move'.format(e))  # Exception won't effect final solution

    return possible, change_func_val


# Methods related to cost function
def get_func_val(zones, zids=None):
    """
    Computes the cost function of a set of zones

    :param zones: Dictionary containing the properties of the zones
    :param zids: List of zones
    :return: Aggregated functional value of the list of zones
    """
    if zids is None:
        zids = [h for h in zones.keys()]

    f = sum(zones[zid]['func_val'] for zid in zids)
    return f


def target_balance(args, members):
    """
    Compute the target balance (capacity utilization) score of a zone (a set of spatial units)

    :param args: Tuple for spatial computation (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schlattr)
    :param members: List of spatial units forming the zone
    :return: Population, Capacity and capacity utilization of the zone
    """
    pop, cap = args[0], args[1]
    p, c, score = 0, 0, 0

    try:
        p = sum([pop[m] for m in members])
        c = sum([cap[m] for m in members])
        score = (p + 0.001) / (c + epsilon)

    except Exception as e:
        print(e)

    f_pop = abs(1 - score)
    return p, c, f_pop


def target_compactness(args, members):
    """
    Compute the target compactness of a zone (a set of spatial units)

    :param args: Tuple for spatial computation (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schlattr)
    :param members: List of spatial units forming the zone
    :return: Area, Perimeter and Compactness score of the zone
    """

    shapes = args[3]
    shape_list = [shapes['geometry'][m] for m in members]
    area, peri, score = 0, 0, 0

    try:
        total = unary_union(shape_list)    # Aggregate the geometry by combining the zones
        area = total.area
        peri = total.length
        score = (4 * math.pi * area) / (peri ** 2)    # IPQ score or Polsby Popper score
    except Exception as e:
        print(e)

    f_compact = 1 - score
    return area, peri, f_compact


def update_property(args, zones, zids=None):
    """
    Update the properties of zones

    :param args: Tuple for spatial computation (population, capacity, adjacency, spas, spas_nbr, sch_spas_nbr, schlattr)
    :param zones: Dictionary containing the properties of the zones
    :param zids: The zone ids whose properties need to be updated
    :return: None
    """
    global w1, w2
    if zids is None:
        zids = [h for h in zones.keys()]

    for zid in zids:
        try:
            members = [m for m in zones[zid]['members']]
            # Get population and capacity statistics
            pop, cap, f1 = target_balance(args, members)
            zones[zid]['Capacity'] = cap
            zones[zid]['Population'] = pop
            zones[zid]['F1'] = f1
            # Get area and perimeter statistics
            area, peri, f2 = target_compactness(args, members)
            zones[zid]['Area'] = area
            zones[zid]['Perimeter'] = peri
            zones[zid]['F2'] = f2

            zones[zid]['func_val'] = w1 * f1 + w2 * f2
        except Exception as e:
            print('{} in computation()'.format(e))
