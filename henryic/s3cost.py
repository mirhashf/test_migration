#!/usr/bin/python

import sys

def s3_storage_cost(store_size):

    # sizes in TBs
    size_tiers  = [1, 49, 450, 500, 4000, 5000]
    # prices in $/GB
    price_tiers = [0.125, 0.110, 0.095, 0.090, 0.080, 0.055]

    if len(size_tiers) != len(price_tiers):
        print 'Number of size tiers does not match number of price tiers.'
        sys.exit(-1)

    chunk_per_tier = [0.0] * len(size_tiers)
    cost_per_tier  = [0.0] * len(size_tiers)

    chunk_left = float(store_size)

    for t in range(len(size_tiers)):
        if t == len(size_tiers)-1:
            chunk_per_tier[t] = chunk_left
        else:
            chunk_per_tier[t] = size_tiers[t] if (chunk_left > size_tiers[t]) else max(0,chunk_left)

        # 1024x for TB to GB (TiB to GiB) conversion
        cost_per_tier[t]  = 1024 * chunk_per_tier[t] * price_tiers[t]
        chunk_left = max(0,chunk_left - size_tiers[t])

    store_cost = sum(cost_per_tier)
    return store_cost


def s3_transfer_cost(transfer_size):

    # sizes in TBs
    size_tiers  = [1/1024, 10, 40, 100, 350]
    # prices in $/GB
    price_tiers = [0, 0.120, 0.090, 0.070, 0.050]

    if len(size_tiers) != len(price_tiers):
        print 'Number of size tiers does not match number of price tiers.'
        sys.exit(-1)

    chunk_per_tier = [0.0] * len(size_tiers)
    cost_per_tier  = [0.0] * len(size_tiers)

    chunk_left = float(transfer_size)

    for t in range(len(size_tiers)):
        if t == len(size_tiers)-1:
            chunk_per_tier[t] = chunk_left
        else:
            chunk_per_tier[t] = size_tiers[t] if (chunk_left > size_tiers[t]) else max(0,chunk_left)

        # 1024x for TB to GB (TiB to GiB) conversion
        cost_per_tier[t]  = 1024 * chunk_per_tier[t] * price_tiers[t]
        chunk_left = max(0,chunk_left - size_tiers[t])

    transfer_cost = sum(cost_per_tier)
    return transfer_cost


def main(args):
    store_cost = s3_storage_cost(args[1])
    print '$' + str(store_cost)


if __name__ == "__main__":
    main(sys.argv)
